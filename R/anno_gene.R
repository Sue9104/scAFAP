library(dplyr)

#' @title Identify Mispriming Sites
#' @description This function identifies regions in the genome that contain
#'   repetitive sequences, which can cause mispriming, and saves the results.
#' @param fa_file The file path to the genome FASTA file.
#' @param anno_dir Output directory.
#' @return Invisibly returns `1`.
#' @importFrom Biostrings readDNAStringSet matchPattern
#' @importFrom stringr str_split_i
#' @import GenomicRanges
find_mispriming_sites <- function(fa_file, anno_dir) {
  suppressPackageStartupMessages({library(GenomicRanges)})
  message("Identifying mispriming sites...")
  # Load and process DNA sequences
  dna_seq <- Biostrings::readDNAStringSet(fa_file)
  names(dna_seq) <- stringr::str_split_i(names(dna_seq), ' ', 1)
  dna_seq <- dna_seq[!grepl('^G|^ML|^J|^K', names(dna_seq))]

  # Helper function to find repetitive stretches
  get_stretch_gr <- function(genome, base = "A", stretchLen = 6, mismatch = 0) {
    stretchA.gr <- c()
    for (seqname in names(genome)) {
      seq.chr = genome[[seqname]]
      stretchA <- Biostrings::matchPattern(
        pattern = paste(rep(base, stretchLen), collapse = ""),
        subject = seq.chr,
        max.mismatch = mismatch
      )
      if (length(stretchA) > 0) {
        stretchA.gr <- c(stretchA.gr, GenomicRanges::GRanges(seqnames = seqname, ranges = stretchA@ranges))
      }
    }
    stretchA.gr <- do.call("c", stretchA.gr)
    GenomicRanges::strand(stretchA.gr) <- ifelse(base == 'A', '+', '-')
    return(stretchA.gr)
  }

  misprimingA <- get_stretch_gr(dna_seq, base = 'A', stretchLen = 6, mismatch = 0)
  misprimingT <- get_stretch_gr(dna_seq, base = 'T', stretchLen = 6, mismatch = 0)
  misprimings <- c(misprimingA, misprimingT)
  misprimings <- GenomicRanges::resize(misprimings, fix = 'center', width = GenomicRanges::width(misprimings) + 10)
  misprimings <- GenomicRanges::reduce(misprimings)
  saveRDS(misprimings, file = stringr::str_glue('{anno_dir}/misprimings.rds'))
  return(1)
}

#' @title Process Genomic Annotations
#' @description This function processes the GTF file to extract gene models
#'   and identifies overlapping genes. It saves intermediate results for
#'   parallel processing.
#' @param gtf_file The file path to the GTF file.
#' @param anno_dir Output directory.
#' @param cores The number of parallel cores to use for processing.
#' @return A list containing processed genomic annotations, including `pc_genes`,
#'   `olp_genes`, `txid2genename`, and `model_file`.
#' @importFrom plyranges read_gff2
#' @importFrom dplyr filter pull
#' @importFrom tibble as_tibble
#' @import GenomicRanges
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom stringr str_glue
#' @importFrom logger log_debug
process_genomic_annotations <- function(gtf_file, anno_dir, cores) {
  message("Processing genomic regions and gene models...")
  # Load and process GTF file
  gencode.gr <- plyranges::read_gff2(gtf_file)

  txid2genename <- gencode.gr %>%
    dplyr::filter(type == "transcript") %>%
    tibble::as_tibble() %>%
    dplyr::pull(gene_name, name = transcript_id)

  all_genes.gr <- gencode.gr %>%
    dplyr::filter(type == "gene") %>%
    dplyr::filter(
      if ('gene_biotype' %in% colnames(GenomicRanges::mcols(.))) {
        gene_biotype == "protein_coding"
      } else {
        gene_type == "protein_coding"
      })
  pc_genes <- all_genes.gr$gene_name
  write.table(dplyr::tibble(gene = pc_genes),
              file = stringr::str_glue("{anno_dir}/genes.list"),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

  olps <- GenomicRanges::findOverlaps(all_genes.gr, drop.self = TRUE, drop.redundant = FALSE)
  olp_genes <- split(pc_genes[olps@to], pc_genes[olps@from])

  genes.gr <- gencode.gr[gencode.gr$gene_name %in% pc_genes]
  genes.gr <- GenomicRanges::split(genes.gr, ~gene_name)

  # Helper function to get UTRs
  get_utr <- function(gr_type, gr) {
    utr.gr <- GenomicRanges::GRanges()
    if ('three_prime_utr' %in% gr$type) {
      utr.gr <- GenomicRanges::reduce(gr[gr$type == gr_type])
    } else {
      if (gr_type == 'five_prime_utr') {
        utr.gr <- GenomicRanges::reduce(gr[(gr$type == "UTR") & (gr$exon_number == 1)])
      } else {
        utr.gr <- GenomicRanges::reduce(gr[(gr$type == "UTR") & (gr$exon_number != 1)])
      }
    }
    return(utr.gr)
  }
  tmp_dir <- stringr::str_glue('{anno_dir}/tmp')
  dir.create(tmp_dir, showWarnings = F, recursive = T)

  elements <- BiocParallel::bplapply(
    genes.gr[pc_genes],
    BPPARAM = BiocParallel::SnowParam(workers = cores, tasks = 500),
    function(gr, tmp_dir) {
      gene.oi <- unique(gr$gene_name)
      region <- range(gr)
      gene_strand <- as.character(GenomicRanges::strand(region)[1])
      gene_chrom <- as.character(GenomicRanges::seqnames(region)[1])
      pas <- GenomicRanges::reduce(
        GenomicRanges::resize(gr[gr$type %in% c('transcript')], fix = "center", width = 2 * 300)
      )
      utr3 <- get_utr('three_prime_utr', gr)
      utr5 <- get_utr('five_prime_utr', gr)

      exon <- GenomicRanges::reduce(gr[gr$type %in% c('exon')])
      intron <- GenomicRanges::setdiff(region, c(utr3, utr5, exon))
      res <- list(
        'name' = gene.oi, 'strand' = gene_strand, 'chrom' = gene_chrom,
        'region' = region, 'pas' = pas,
        'utr3' = utr3, 'utr5' = utr5,
        'exon' = exon, 'intron' = intron
      )
      saveRDS(res, file = stringr::str_glue('{tmp_dir}/{gene.oi}.rds'))
      return(1)
    },
    tmp_dir = tmp_dir
  )

  # Generate gene model bed file
  model_file <- stringr::str_glue('{anno_dir}/gene_model.bed')
  cmd <- stringr::str_glue(
    'gtfToGenePred {gtf_file} {anno_dir}/genePred.tsv && ',
    'genePredToBed {anno_dir}/genePred.tsv {model_file} && ',
    'rm {anno_dir}/genePred.tsv'
  )
  logger::log_debug(cmd)
  system(cmd, intern = FALSE)

  anno_res <- list(
    pc_genes = pc_genes,
    olp_genes = olp_genes,
    txid2genename = txid2genename,
    pc_genes = pc_genes,
    model_file = model_file
  )
  return(anno_res)
}

#' @title Save All Output Files
#' @description This function saves the processed genomic information to various
#'   output files, including gene model BED files and gene-specific RData files,
#'   by leveraging parallel processing.
#' @param anno_res A list of genomic annotations from `process_genomic_annotations`.
#' @param anno_dir Output directory.
#' @param cores The number of parallel cores for saving individual gene files.
#' @return Invisibly returns `NULL`.
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom stringr str_glue
#' @import GenomicRanges
#' @importFrom dplyr tibble
save_output_files <- function(anno_res, anno_dir, cores) {
  suppressPackageStartupMessages({library(GenomicRanges)})
  message("Saving final output files...")

  anno_subdir <- stringr::str_glue('{anno_dir}/genes')
  dir.create(anno_subdir, showWarnings = FALSE, recursive = TRUE)

  # Parallel saving of individual gene files
  BiocParallel::bplapply(
    anno_res$pc_genes,
    BPPARAM = BiocParallel::SnowParam(workers = cores, tasks = 500),
    function(gene.oi, olp_genes, txid2genename, model_file, anno_subdir) {
      suppressPackageStartupMessages({
        library(GenomicRanges)
      })
      # message(gene.oi)
      gene_struc <- readRDS(stringr::str_glue('{anno_dir}/tmp/{gene.oi}.rds'))

      misprimings <- readRDS(stringr::str_glue('{anno_dir}/misprimings.rds'))
      mispriming <- misprimings[misprimings %over% gene_struc$region]
      rm(misprimings)
      olps <- olp_genes[[gene.oi]]
      others <- GenomicRanges::GRanges()
      if (!is.null(olps)) {
        others <-
          lapply(olps, function(g) {
            gene_struc <- readRDS(stringr::str_glue('{anno_dir}/tmp/{g}.rds'))
            c(gene_struc$pas, gene_struc$utr)
          }) %>%
          GenomicRanges::GRangesList() %>%
          unlist() %>%
          GenomicRanges::reduce()
      }

      df <- read.table(model_file)
      txs <- names(txid2genename)[txid2genename == gene.oi]
      gene_model <- stringr::str_glue('{anno_subdir}/{gene.oi}.gene_model.bed')
      write.table(df[df$V4 %in% txs, ], file = gene_model,
                  sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

      gene_info <- c(
        gene_struc,
        list(
          "olps" = olps,
          "others" = others,
          "mispriming" = mispriming,
          'model_file' = gene_model
        )
      )
      saveRDS(gene_info,
              file = stringr::str_glue('{anno_subdir}/{gene.oi}.gene_info.rds'))
      return(NULL)
    },
    olp_genes = anno_res$olp_genes,
    txid2genename = anno_res$txid2genename,
    model_file = anno_res$model_file,
    anno_subdir = anno_subdir
  )
  message("Output files saved.")
}



#' @title Prepare Annotation Files
#' @description This is the main orchestrator function that runs a series of
#'   genomic analysis steps based on user-specified options. It handles
#'   mispriming site identification, genomic region processing, and saving
#'   of all output files.
#' @param fa_file The file path to the genome FASTA file.
#' @param gtf_file The file path to the GTF file.
#' @param anno_dir The output directory for annotation files.
#' @param cores The number of parallel cores to use.
#' @return Invisibly returns `NULL`.
#' @export
prepare_anno <- function(fa_file, gtf_file, anno_dir,
                         cores = 60
) {
  anno_dir <- normalizePath(anno_dir)
  dir.create(anno_dir, showWarnings = F, recursive = T)

  # --- Step 1: Analyze mispriming ---
  find_mispriming_sites(fa_file, anno_dir)

  # --- Step 2: Process genomic regions ---
  anno_res <- process_genomic_annotations(gtf_file, anno_dir, cores)

  # --- Step 3: Save final output ---
  save_output_files(anno_res, anno_dir, cores)

}
