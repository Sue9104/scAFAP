#' @title Adaptive PolyA Site Detection
#' @description This is the main function for the PAS detection workflow.
#' @param gene.oi A character string for the gene of interest (e.g., "MYC").
#' @param inbam The file path to the input BAM file.
#' @param outdir The output directory for saving all results.
#' @param sample A character string specifying the sample name.
#' @param anno_dir The directory containing pre-processed annotation files.
#' @param plot A logical value indicating whether to generate diagnostic plots.
#' @param run_mode A character string specifying the analysis mode (e.g., 'fast').
#' @param scaler An integer for scaling parameters.
#' @return An invisible integer `1` on successful completion, or `NULL` if any step fails.
#' @export
#' @importFrom stringr str_glue
#' @importFrom logger log_info
#' @import GenomicRanges
adaptive_detect_pas <- function(gene.oi, inbam, outdir, sample,
                                anno_dir,
                                plot = T,
                                run_mode = 'fast', scaler = 2) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(magrittr)
    library(purrr)
    library(GenomicRanges)
  })

  # --- Step 1: Setup and load required data ---
  prefix <- stringr::str_glue("{outdir}/tmps/{gene.oi}.{sample}")
  dir.create(stringr::str_glue("{outdir}/tmps"), showWarnings = F, recursive = T)

  # Load gene information and other necessary environment data.
  env_data <- setup_and_load_data(gene.oi, outdir, sample, anno_dir)
  if (is.null(env_data)) return(NULL)
  gene_info <- env_data$gene_info

  # --- Step 2: Calculate coverage and junctions ---
  reads <- extract_bam_coverage(gene_info, inbam, prefix, env_data$samtools)
  if (reads < 200) return(NULL)

  junc_res <- analyze_splicing_junctions(gene_info, prefix, env_data$junc_anno, plot = T)
  if (is.null(junc_res)) return(NULL)
  newCoords <- junc_res$newCoord
  oldCoords <- as.integer(names(newCoords))

  # --- Step 3: Perform peak calling ---
  data <- perform_peak_calling(prefix, newCoords)
  cov_data <- data$cov_data
  raw_peaks <-data$raw_peaks
  if (is.null(data$raw_peaks)) return(NULL)

  # --- Step 4: Annotate peaks with poly-A reads ---
  peak_data <- annotate_and_refine_peaks(gene_info, data, prefix, newCoords, oldCoords, plot = T)
  peaks <- peak_data$peaks
  print(peak_data$plot)
  cov_data %>% ggplot(aes(x=start)) +
    geom_point(aes(y=counts)) + geom_point(aes(y=splines), color = 'blue') +
    geom_vline(xintercept = raw_peaks$summit, color = 'red')

  if (is.null(peaks)) return(NULL)

  # --- Step 5: Run MCMC analysis (if needed) and generate final output ---
  do_mcmc <- run_mcmc(
    data$cov_data, peaks,
    prefix, env_data$py3, env_data$mcmc_py
  )
  finalize_results(gene_info, peaks, do_mcmc, data$cov_data, oldCoords, prefix)
  return(1)
}


#' @title Setup Environment and Load Data
#' @description Prepares the environment and loads necessary data for the analysis,
#'   including paths to external tools like `samtools` and `junction_annotation.py`.
#' @param gene.oi The gene of interest.
#' @param outdir The output directory.
#' @param sample The sample name.
#' @param anno_dir The directory containing annotation files.
#' @return A list containing `gene_info` and paths to external tools, or `NULL` if
#'   the annotation file is not found.
#' @importFrom stringr str_glue
#' @importFrom logger log_warn
setup_and_load_data <- function(gene.oi, outdir, sample, anno_dir) {
  anno_file <- stringr::str_glue('{anno_dir}/genes/{gene.oi}.gene_info.rds')
  if (!file.exists(anno_file)) {
    logger::log_warn(gene.oi, " not find annotation !!!")
    return(NULL)
  }
  gene_info <- readRDS(anno_file)

  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/miniconda/envs/scpa/bin", sep = ":"))
  samtools <- system('which samtools', intern = T)
  junc_anno <- system('which junction_annotation.py', intern = T)

  py3 <- '~/miniconda/envs/pymc/bin/python'
  mcmc_py <- system.file("python", "mcmc.py", package = "scAFAP")
  # mcmc_py <- '~/pipelines/scAFAP/inst/python/mcmc.py'
  env_data <- list(gene_info = gene_info,
                   samtools = samtools,
                   junc_anno = junc_anno,
                   py3 = py3, mcmc_py = mcmc_py)
  return(env_data)
}


# ---
#' @title Extract Coverage from BAM
#' @description Extracts reads for a specific gene from a BAM file and generates coverage data.
#' @param gene_info A list containing information about the gene of interest.
#' @param inbam The file path to the input BAM file.
#' @param prefix The output file path prefix for intermediate files.
#' @param samtools The system path to the `samtools` executable.
#' @return An invisible integer `1` on successful completion.
#' @importFrom stringr str_glue str_remove
#' @importFrom logger log_info log_debug
extract_bam_coverage <- function(gene_info, inbam, prefix, samtools) {
  logger::log_info(stringr::str_glue("Extract {gene_info$name} reads from bam"))
  tagS <- ifelse(gene_info$strand == "+", "-F 20", "-f 16 -F 4")
  region <- stringr::str_remove(as.character(gene_info$region), ':[+-]$')
  gene_bam <- stringr::str_glue("{prefix}.bam")
  cmd <- stringr::str_glue(
    '{samtools} view -h -@ 2 {inbam} {region} {tagS} -o {gene_bam} -d xf:25 ',
    ' && {samtools} index {gene_bam}'
  )
  logger::log_debug(cmd)
  system(cmd, intern = T)
  reads <- bam_to_cov(gene_bam, prefix, gene_info$mispriming)
  return(reads)
}


# ---
#' @title Analyze Splicing Junctions
#' @description Runs an external tool (`RSeQC` or similar) to analyze splicing junctions,
#'   processes the output, and generates a diagnostic plot.
#' @param gene_info A list containing information about the gene of interest.
#' @param prefix The output file path prefix for intermediate and final files.
#' @param junc_anno The system path to the `junction_annotation.py` executable.
#' @param plot A logical value indicating whether to generate a PDF plot.
#' @return A list containing the splicing analysis results, or `NULL` if no junctions
#'   are found.
#' @importFrom stringr str_glue str_remove
#' @importFrom logger log_info log_debug
analyze_splicing_junctions <- function(gene_info, prefix, junc_anno, plot = F) {
  logger::log_info(stringr::str_glue("Calculate junctions of {gene_info$name}"))
  gene_bam <- stringr::str_glue("{prefix}.bam")

  # Construct the command to run the external junction annotation script.
  cmd <- stringr::str_glue("{junc_anno} -i {gene_bam} -r {gene_info$model_file} -o {prefix}")
  logger::log_debug(cmd)
  system(cmd, intern = F)

  junc.file <- stringr::str_glue("{prefix}.junction.xls")
  read_cov.file <- stringr::str_glue("{prefix}.cov.reads.csv")

  # Skip junctions if spliced > unspliced.
  junc_res <- jump_splicing_region(
    gene_info$name, gene_info$region, gene_info$utr3, gene_info$exon,
    junc.file, read_cov.file,
    minJuncCount = 20
  )

  # Generate and save a PDF plot if requested.
  if (plot) {
    pdf(stringr::str_glue("{prefix}.junction.pdf"), 6, 8)
    print(junc_res$plot)
    dev.off()
  }
  return(junc_res)
}

# ---
#' @title Perform Peak Calling
#' @description Reads coverage data, smooths it using splines, and calls peaks from the smoothed signal.
#' @param prefix The file path prefix for the input coverage data (`.cov3.all.csv`).
#' @param newCoords A named integer vector for coordinate mapping.
#' @return A list containing the raw peaks (`raw_peaks`), the coverage data (`cov_data`), and candidate regions (`candidates`), or `NULL` if no peaks are found.
#' @importFrom stringr str_glue
#' @import GenomicRanges
#' @importFrom dplyr if_else left_join mutate select arrange
#' @importFrom logger log_debug
perform_peak_calling <- function(prefix, newCoords) {
  cov3.df <- read.csv(stringr::str_glue("{prefix}.cov3.all.csv"))
  cov3.df <- df_to_new_coords(cov3.df, newCoords)
  reads <- sum(cov3.df$counts)
  if (reads < 100) return(NULL)

  cov3.gr <- GenomicRanges::makeGRangesFromDataFrame(cov3.df, keep.extra.columns = T)
  # zoom-in and smooth
  min_counts <- stats::quantile(cov3.df$counts[cov3.df$counts >= 3], 0.5)
  candidates.gr <- GenomicRanges::reduce(cov3.gr[cov3.gr$counts >= min_counts], min.gapwidth = 50)
  candidates.gr <- candidates.gr[GenomicRanges::width(candidates.gr) >= 70]
  if (length(candidates.gr) == 0) return(NULL)
  candidates.gr <- GenomicRanges::reduce(
    GenomicRanges::resize(candidates.gr, fix = 'center', width = GenomicRanges::width(candidates.gr) + 600)
  )
  GenomicRanges::start(candidates.gr) <-
    dplyr::if_else(GenomicRanges::start(candidates.gr) < 1, 1, GenomicRanges::start(candidates.gr))

  cov_data <- cov3.df[cov3.gr %within% candidates.gr, ] %>% arrange(start)
  zoom_in <- 1
  candidates_reads <- sum(cov_data$counts)
  use_cv <- sum(GenomicRanges::width(candidates.gr)) > 5000
  logger::log_debug('Smoothing, cv = ', use_cv,
                    ' | width = ', sum(GenomicRanges::width(candidates.gr)),
                    ' | reads = ', candidates_reads)
  if (candidates_reads < 5000) {
    zoom_in <- ceiling(5000 / 2 / candidates_reads)
    cov_data <- flank_cov(GenomicRanges::makeGRangesFromDataFrame(cov_data, keep.extra.columns = T), ex = zoom_in)
    zoom_in <- zoom_in * 2
  }
  cov_data$splines <- run_splines(cov_data$start, cov_data$counts, cv = use_cv)
  cov3.df <-
    left_join(cov3.df, cov_data %>% dplyr::select(start, splines), by = "start") %>%
    mutate(splines = splines / zoom_in)

  # find peaks
  raw_peaks <- find_flexible_peaks(cov_data$start, cov_data$splines,
                                   min_peak_height = max(3, zoom_in),
                                   min_peak_width = 100)
  if (nrow(raw_peaks) == 0) return(NULL)
  raw_peaks <- dplyr::mutate(raw_peaks, seqnames = cov3.df$seqnames[1], strand = cov3.df$strand[1], .before = 1)
  print(raw_peaks)

  data <- list(raw_peaks = raw_peaks, cov_data = cov3.df, candidates = candidates.gr)
  return(data)
}


# ---
#' @title Annotate and Refine Peaks
#' @description Adds poly-A reads to peaks and refines their positions.
#' @param gene_info A list containing information about the gene of interest.
#' @param data A list returned by `perform_peak_calling`, containing `raw_peaks`, `cov_data`, and `candidates`.
#' @param prefix The output file path prefix.
#' @param newCoords A named integer vector for coordinate mapping.
#' @param oldCoords A named integer vector for coordinate mapping.
#' @param plot A logical value indicating whether to generate a PDF plot.
#' @return A list containing the refined peaks (`peaks`) and a diagnostic plot (`plot`), or `NULL` if no valid peaks are found.
#' @importFrom stringr str_glue
#' @importFrom purrr map_vec
#' @importFrom tibble tibble
#' @import GenomicRanges
#' @importFrom logger log_info log_debug log_warn
#' @import ggplot2
annotate_and_refine_peaks <- function(gene_info, data, prefix, newCoords, oldCoords, plot = F) {
  raw_peaks <- data$raw_peaks
  cov_data <- data$cov_data

  # evalute sigma
  sigma <- purrr::map_vec(raw_peaks$summit, function(x) {
    ids <- which((cov_data$start >= x - 125) & (cov_data$start <= x + 125))
    unlist(lapply(ids, function(i) rep(cov_data[i, "start"], cov_data[i, "counts"]))) %>% sd()
  }) %>% round() %>% min()
  sigma <- max(30, min(sigma, 100))
  logger::log_info("predicted sigma: ", sigma)

  # read polyA reads
  cov3_pa_file <- stringr::str_glue("{prefix}.cov3.pa.csv")
  pa_reads <- 0
  preads <- integer(0)
  if (file.exists(cov3_pa_file)) {
    cov3pa.df <- df_to_new_coords(read.csv(cov3_pa_file), newCoords)
    candidates_pa <- GenomicRanges::makeGRangesFromDataFrame(cov3pa.df) %within% data$candidates
    pa_reads <- sum(cov3pa.df[candidates_pa, 'counts'])
    preads <-
      cov3pa.df[candidates_pa, ] %>% dplyr::arrange(end) %>%
      group_by(g = cut_width(start, width = 50, boundary = 0)) %>%
      arrange(desc(counts)) %>% dplyr::slice_head(n = 1) %>% ungroup() %>%
      dplyr::group_by(g = cumsum(c(0, diff(end) >= 50))) %>%
      dplyr::arrange(dplyr::desc(counts)) %>% dplyr::slice_head(n = 1) %>%
      dplyr::pull(end)
    logger::log_debug("PAS: {paste(preads, collapse = ', ')}")
  }
  q <- round(pa_reads / sum(cov_data$counts), 5)
  logger::log_info(stringr::str_glue("{basename(prefix)} polyA reads ratio: {q} ({pa_reads}/{sum(cov_data$counts)})"))

  # add preads to peaks
  q <- ifelse(q == 0, 0.005, q)
  d <- round(qnorm(1 - q) * sigma) - 20
  d <- ifelse(gene_info$strand == "+", d, -d)
  find_summit <- function(pread) {
    final_summit <- NA_integer_
    valids <- (sign(d)*(pread - raw_peaks$summit) >= 50) & (sign(d)*(pread - raw_peaks$summit) <= 200)
    if (sum(valids) > 0){
      if (sum(valids) == 1){
        summit <- raw_peaks$summit[which(valids)]
      } else if (sum(valids) > 1) {
        ids <- which(valids)
        dists <- pread - d - raw_peaks$summit[ids]
        if (any(abs(dists) < 100)){
          summit <- raw_peaks$summit[ ids[which.min(abs(dists))] ]
        }
      }
      if (preads[which.min(abs(summit + d - preads))] == pread){
        final_summit <- summit
      }
    }
    return(final_summit)
  }

  pos2splines <- cov_data %>% dplyr::pull(splines, name = start)
  pos2counts <- cov_data %>% dplyr::pull(counts, name = start)
  preads.df <- tibble::tibble(
    pread = preads,
    counts = pos2counts[as.character(preads)],
    in_peak = sapply(preads - d, function(x) any((raw_peaks$start <= x) & (raw_peaks$end >= x))),
    splines = pos2splines[as.character(preads - d)],
    summit = sapply(preads, find_summit)
  ) %>% dplyr::filter(pread != 0)

  min_pa_counts <- stats::quantile(cov_data[cov_data$counts > 1, 'splines'], 0.5, na.rm = T)
  preads.df <-
    preads.df %>%
    dplyr::filter( (!is.na(summit)) | (in_peak & (counts >= 3)) | ((splines > min_pa_counts) & (counts >= 3)) ) %>%
    dplyr::mutate(summit = dplyr::if_else(is.na(summit), pread - d, summit)) %>%
    dplyr::group_by(summit) %>%
    dplyr::arrange(summit, dplyr::desc(splines)) %>%
    slice_head(n=1) %>% ungroup()

  peaks <-
    dplyr::full_join(raw_peaks, dplyr::select(preads.df, pread, summit), by = "summit") %>%
    dplyr::mutate(pa = dplyr::if_else(is.na(pread), summit + d, pread),
                  is_peak = dplyr::if_else(is.na(height), F, T),
                  seqnames = dplyr::coalesce(seqnames, raw_peaks$seqnames[1]),
                  strand = dplyr::coalesce(strand, gene_info$strand),
                  start = dplyr::coalesce(start, summit - 3 * sigma),
                  end = dplyr::coalesce(end, summit + 3 * sigma),
                  width = dplyr::coalesce(width, 6 * sigma)) %>%
    dplyr::mutate(start = dplyr::if_else((gene_info$strand == "-") & (!is.na(pread)), pread - 20, start),
                  end = dplyr::if_else((gene_info$strand == "+") & (!is.na(pread)), pread + 20, end)) %>%
    dplyr::mutate(start = dplyr::if_else(start <= 0, 1, start),
                  end = dplyr::if_else(end > max(newCoords), max(newCoords), end),
                  width = end - start + 1,
                  pa = dplyr::case_when(
                    pa <= 0 ~ 1,
                    pa >= max(newCoords) ~ max(newCoords),
                    .default = pa
                  ))

  # gene belongs
  pa.gr <- GenomicRanges::GRanges(seqnames = gene_info$chrom, strand = gene_info$strand,
                                  ranges = IRanges::IRanges(start = oldCoords[peaks$pa]))
  valids <- which(! pa.gr %within% gene_info$others)
  if (length(valids) == 0) {
    logger::log_warn("No valid peaks: ", gene.oi)
  } else {
    peaks <- peaks[valids, ]
    peaks$region <- 'genic'
    peaks$region[pa.gr %within% gene_info$exon] <- 'exonic'
    peaks$region[pa.gr %within% gene_info$intron] <- 'intronic'
    peaks$region[pa.gr %within% gene_info$utr3] <- 'utr'
  }
  peaks <-
    peaks %>% dplyr::arrange(summit) %>%
    dplyr::mutate(gene = gene_info$name, id = 1:n(), .before = 1) %>%
    dplyr::mutate(d1 = d, d2 = pread - summit, polyA_ratio = q, sigma = sigma,
                  s1 = c(1000, diff(pa)), s2 = c(diff(pa), 1000),
                  do_mcmc = if_else(s1 < 350 | s2 < 350, T, F)) %>%
    dplyr::mutate(end = if_else((strand == '-') & (width > 350) & (!do_mcmc), start + 349, end),
                  start = if_else((strand == '+') & (width > 350) & (!do_mcmc), end - 349, start),
                  width = end - start + 1)
  peaks$total <- apply(peaks, 1, function(x) {
    ids <- (cov_data$start >= x[["start"]]) & (cov_data$start <= x[["end"]])
    sum(cov_data[ids, 'counts'])
  })
  write.csv(peaks, file = stringr::str_glue("{prefix}.peaks.csv"), quote = F, row.names = F)
  print(peaks)

  # Plot
  if (plot) {
    pdf(stringr::str_glue("{prefix}.peak.pdf"), 9, 3)
    preads_info <- stringr::str_glue("preads: {paste(preads, collapse = ' ')}")
    ymax <- stats::quantile(cov_data$counts, 0.999)
    p <-
      cov_data %>%
      ggplot2::ggplot(ggplot2::aes(x = start)) +
      ggplot2::geom_point(ggplot2::aes(y = counts), size = .2) +
      ggplot2::geom_line(ggplot2::aes(y = splines), color = 'red', linewidth = 1, na.rm = T) +
      ggplot2::geom_vline(xintercept = preads, color = 'blue', linetype = "dashed", linewidth = .6) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(title = stringr::str_glue("{gene_info$name}:{gene_info$strand}"),
                    caption = preads_info) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::coord_cartesian(ylim = c(NA, ymax), xlim = range(peaks$start - 50, peaks$end + 50))
    if (nrow(peaks) > 0) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(peaks$start, peaks$end), color = 'green') +
        ggplot2::geom_vline(xintercept = peaks$summit, color = 'orange')
    }
    print(p)
    dev.off()
  }

  return(list("peaks" = peaks, "plot" = p))
  # return(peaks)
}


# ---
#' @title Run MCMC Analysis
#' @description Prepares and runs the MCMC python script for multi-peak genes.
#' @param cov_data A data frame of coverage data.
#' @param peaks A data frame of refined peaks.
#' @param prefix The output file path prefix.
#' @param py3 The system path to the Python 3 executable.
#' @param mcmc_py The file path to the MCMC Python script.
#' @param run_mode A character string specifying the analysis mode (e.g., 'fast' or 'precise').
#' @param scaler An integer for scaling parameters.
#' @return A logical vector indicating which peaks underwent MCMC analysis.
#' @importFrom stringr str_glue
#' @import GenomicRanges
#' @importFrom logger log_info
run_mcmc <- function(cov_data, peaks, prefix, py3, mcmc_py,
                              run_mode = 'fast', scaler = 2) {
  do_mcmc <- (peaks$s1 < 350 | peaks$s2 < 350)
  if (sum(do_mcmc) > 0) {
    cov3.gr <- GenomicRanges::makeGRangesFromDataFrame(cov_data, keep.extra.columns = T)
    peak.gr <- GenomicRanges::makeGRangesFromDataFrame(peaks[do_mcmc,])
    hits <- GenomicRanges::findOverlaps(cov3.gr, peak.gr, maxgap = 50)@from
    olps.gr <- cov3.gr[unique(hits)]
    total_reads <- sum(olps.gr$counts)
    peak_width <- sum(GenomicRanges::width(GenomicRanges::reduce(olps.gr, min.gapwidth = 50)))
    bp <- ifelse(sum(olps.gr$counts) > 2000, 0, 3)
    mcmc_data <- flank_cov(olps.gr, ex = bp)
    if (scaler != 0) {
      s <- max(peak_width * scaler, 200) / total_reads
      mcmc_data <- mcmc_data %>% dplyr::mutate(counts = round(counts * s)) %>% dplyr::filter(counts > 0)
    }
    write.table(mcmc_data, stringr::str_glue("{prefix}.mcmc_input.csv"),
                sep = ",", quote = F, row.names = F, col.names = F)

    cmd <- stringr::str_glue(
      "{py3} {mcmc_py} ",
      "{prefix}.mcmc_input.csv {sum(do_mcmc)} ",
      "--sigma {peaks$sigma[1]} ",
      "--init {paste(peaks[do_mcmc,'summit'], collapse = ' ')} ",
      "--prefix {prefix}.mcmc \n"
    )
    cat(cmd, "\n", file = stringr::str_glue("{dirname(prefix)}/mcmc_cmd.sh"), append = TRUE)
    if (run_mode == 'precise') cmd <- stringr::str_glue("{cmd} --infer")
    logger::log_info(cmd)
    # system(cmd, intern = F)
  }
  return(do_mcmc)
}


# ---
#' @title Finalize and Save Results
#' @description Plots and saves the final results of the analysis.
#' @param gene_info A list containing information about the gene of interest.
#' @param peaks A data frame of refined peaks.
#' @param do_mcmc A logical vector indicating which peaks underwent MCMC analysis.
#' @param cov_data A data frame of coverage data.
#' @param oldCoords A named integer vector for coordinate mapping.
#' @param prefix The output file path prefix.
#' @return Returns `NULL` invisibly.
#' @importFrom stringr str_glue str_remove_all
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble tibble column_to_rownames
#' @importFrom purrr set_names
#' @importFrom logger log_error
finalize_results <- function(gene_info, peaks, do_mcmc, cov_data, oldCoords, prefix) {
  gene.oi <- gene_info$name

  # Write results
  final.df <- tibble::tibble(
    gene_id = gene_info$name,
    strand = gene_info$strand, seqnames = gene_info$seqnames,
    id = peaks$id, region = peaks$region, summit = peaks$summit,
    pread = peaks$pread, pread_ratio = peaks$polyA_ratio,
    mcmc = do_mcmc
  )
  # prob for each pos
  id2mid <- cumsum(do_mcmc) - 1
  mid2id <- seq_along(do_mcmc)
  names(mid2id) <- as.character(id2mid)
  probs.df <- tibble()
  if (sum(do_mcmc) > 0){
    mcmc_file <- stringr::str_glue("{prefix}.mcmc.model.csv")
    if (!file.exists(mcmc_file)) {
      logger::log_error('MCMC is not finished !!!')
      return(NULL)
    }
    mcmc_data <- read.csv(mcmc_file)
    post_mean <- purrr::set_names(mcmc_data$mean, nm = mcmc_data$X)

    probs.mcmc <-
      read.table(stringr::str_glue("{prefix}.mcmc.probs.txt"), header = T) %>%
      dplyr::distinct() %>%
      pivot_longer(cols=2:sum(do_mcmc), names_to = 'id', values_to = 'prob') %>%
      mutate(id = str_remove_all(id, 'Cluster_|_Prob'), id = mid2id[id])
    colnames(probs.mcmc)[1] <- "start"
    probs.df <- bind_rows(probs.df, probs.mcmc)
    mcmc_reads <- sum(cov_data[cov_data$start %in% probs.mcmc$start, 'counts'])
  }
  if (sum(do_mcmc) != nrow(peaks)){
    probs.not_mcmc <-
      apply(peaks[!do_mcmc, ], 1, function(x){
        tibble(start = seq(x[["start"]], x[["end"]]), id = x[["id"]], prob = 1)
      }) %>%
      purrr::list_rbind()
    probs.df <- bind_rows(probs.df, probs.not_mcmc)
  }

  for (i in seq_along(do_mcmc)){
    if (do_mcmc[i]){
      mid <- id2mid[i]
      final.df[i, 'pa_weights']  <- post_mean[stringr::str_glue('w[{mid}]')]
      final.df[i, 'pa_sigma']  <- post_mean[stringr::str_glue('sigma[{mid}]')]
      final.df[i, 'pa_mu']  <- post_mean[stringr::str_glue('mu[{mid}]')]
      final.df[i, 'pa_counts']  <- peaks[i, 'pa_weights'] * mcmc_reads
    } else {
      final.df[i, 'pa_weights']  <- 1
      final.df[i, 'pa_sigma']  <- peaks[i, 'sigma']
      final.df[i, 'pa_mu']  <- peaks[i, 'summit']
      final.df[i, 'pa_counts']  <- peaks[i, 'total']
    }
  }

  # Finalize columns and coordinates
  final.df <- final.df %>%
    dplyr::mutate(pa_weights = round(pa_counts/sum(pa_counts), 3)) %>%
    dplyr::mutate(
      d = dplyr::if_else(is.na(pread),
                         stats::qnorm(1 - 0.01) * pa_sigma - 20,
                         stats::qnorm(1 - pread_ratio) * pa_sigma - 20),
      d = round(d),
      pa = dplyr::if_else(strand == '+', pa_mu + d, pa_mu - d),
      pa = dplyr::if_else(is.na(pread), pa, pread)
    ) %>%
    dplyr::mutate(
      pa = oldCoords[pa],
      summit = oldCoords[summit],
      pa_mu = oldCoords[pa_mu]
    )
  print(final.df)
  write.csv(final.df, stringr::str_glue("{prefix}.pas.stats.csv"), quote = F, row.names = F)

  # Cell counts
  id2pas <- final.df %>% mutate(pas = str_glue('{gene.oi}:{pa}')) %>% pull(pas, name = id)
  probs.df <- probs.df %>% mutate(pas = id2pas[id], start = oldCoords[start])
  cov3_cb.df <-
    read.csv(stringr::str_glue("{prefix}.cov3.cb.csv")) %>%
    filter(start %in% probs.df$start) %>%
    left_join(probs.df, by = 'start', relationship = "many-to-many")
  matrix <-
    cov3_cb.df %>%
    group_by(CB, pas) %>% summarise(counts = sum(counts * prob)) %>% ungroup() %>%
    pivot_wider(id_cols = "pas", names_from = "CB", values_from = "counts", values_fill = 0) %>%
    tibble::column_to_rownames('pas')
  write.csv(matrix, stringr::str_glue("{prefix}.pas.matrix.csv"), quote = F, row.names = T)
}
