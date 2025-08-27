#' Calculate Gene Coverage from a BAM File
#'
#' This function uses the `pandepth` tool to calculate gene sequencing coverage from
#' an input BAM file. It then reads the coverage statistics, filters for genes
#' that meet predefined thresholds (sites >= 100 and total_depth >= 200),
#' and returns the names of the genes that satisfy these criteria.
#'
#' @param inbam A string specifying the path to the input BAM file.
#' @param anno_dir A string specifying the path to the annotation directory
#'   containing the `coding_genes.bed` file.
#' @param outdir A string specifying the path to the output directory where
#'   temporary files will be stored.
#' @return A character vector containing the names of the genes that meet the
#'   coverage filtering criteria.
#' @export
#' @importFrom stringr str_glue
bam_to_genes <- function(inbam, anno_dir, outdir, sample, cpm = 5){
  # calculate cov
  genes.bed <- stringr::str_glue("{anno_dir}/coding_genes.bed")
  genes.df <- read.table(genes.bed)
  id2genes <- setNames(genes.df$V6,
                       stringr::str_glue('{genes.df$V1}_{genes.df$V2}_{genes.df$V3}'))

  cmd <- stringr::str_glue('pandepth -t 8 -d 5 -i {inbam} -b  {genes.bed} ',
                          '-o {outdir}/{sample}.gene_cov')
  logger::log_debug(cmd)
  system(cmd, intern = F)

  covs.df <- read.table(stringr::str_glue('{outdir}/{sample}.gene_cov.bed.stat.gz'),
                        col.names = c("chrom", "start", "end", "id", "width",
                                   "sites", "total_depth", "percent", "gene_depth"),
                        sep = '\t')
  covs.df$cpm_value <- round(covs.df$total_depth / sum(covs.df$total_depth) * 10^6, 2)
  ids <- covs.df[covs.df$cpm_value > cpm,"id"]
  genes <- as.vector(id2genes[ids])
  logger::log_debug(sample, " available genes: ", length(genes))
  return(genes)
}



#' @title Convert Rle to a Single-Base data.table
#' @description This function takes an Rle object, converts it to a GenomicRanges
#'   object, and then expands it into a data.table with one row for each base pair.
#'   It efficiently filters out regions with zero coverage.
#' @param cov.rle An Rle object containing coverage data.
#' @return A data.table with the following columns: 'seqnames', 'start', 'end',
#'   'strand', and 'counts'. Each row represents a single base pair.
#' @importFrom data.table as.data.table data.table
rle_to_dt <- function(cov.rle){
  dt <- data.table::as.data.table(as(cov.rle, "GRanges"))
  dt <- dt[dt$score > 0,]

  # Expand the data.table so that each row represents a single base pair.
  dt <-
    dt[ , data.table::data.table(
      seqnames = seqnames,
      start = seq(from = start, to = end),
      end = seq(from = start, to = end),
      strand = '*',
      counts = score
    ), by = 1:nrow(dt)]
  dt <- dt[,-1]

  return(dt)
}

#' @title Convert BAM to Coverage and Reads Data
#' @description This function processes a BAM file, removes reads from mispriming sites,
#'   and calculates read coverage and 3' stop-site counts, saving the output to CSV files.
#' @param gene_bam A character string specifying the path to the gene-level BAM file.
#' @param prefix A character string specifying the prefix for output file names.
#' @param mispriming A GenomicRanges object containing the locations of mispriming sites to be filtered.
#' @return The function invisibly returns `NULL` and saves several CSV files:
#'   \itemize{
#'     \item `.reads.csv`: Processed reads data.
#'     \item `.cov.reads.csv`: Coverage of all processed reads.
#'     \item `.cov3.all.csv`: Counts of 3' stop-sites for all reads.
#'     \item `.cov3.pa.csv`: Counts of 3' stop-sites for reads with a "pa" tag.
#'     \item `.cov3.cb.csv`: Counts of 3' stop-sites per cell barcode.
#'   }
#' @export
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse
#' @import GenomicRanges
#' @importFrom S4Vectors decode
#' @importFrom stringr str_glue
bam_to_cov <- function(gene_bam, prefix, mispriming) {
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
  })
  # Step 1: Read the BAM file and specify tags
  param <- Rsamtools::ScanBamParam(
    what = c("rname", "pos", "qwidth"),
    tag = c("CB", "pa")
  )
  gal <- GenomicAlignments::readGAlignments(gene_bam, param = param)
  GenomeInfoDb::seqlevels(gal) <- GenomeInfoDb::seqlevelsInUse(gal)
  if (length(gal) == 0) return(0)

  # Step 2: Remove reads overlapping with mispriming sites
  reads3_prime_end <- GenomicRanges::resize(GenomicRanges::granges(gal), fix = 'end', width = 1)
  gal <- gal[!reads3_prime_end %over% mispriming]
  reads <- length(gal)
  if (length(gal) == 0) return(0)

  # Step 3: Convert reads to a data.table for efficient manipulation
  tags <- S4Vectors::mcols(gal)
  dt <- data.table::as.data.table(GenomicRanges::granges(gal))

  # Add tag information to the data.table
  dt[, `:=`(pa = tags$pa, CB = tags$CB)]

  # Write processed reads to file
  # reads_to_save <- dt[, .(seqnames, start, end, strand, CB)]
  # data.table::fwrite(reads_to_save, file = stringr::str_glue("{prefix}.reads.csv"), row.names = FALSE)

  # Step 4: Calculate and write total coverage (optimized)
  # Directly calculate coverage on the GAlignments object for efficiency
  reads_granges <- GenomicRanges::granges(gal)
  reads_rle <- GenomicRanges::coverage(reads_granges)
  reads.dt <- data.table::as.data.table(as(reads_rle, "GRanges"))
  reads.dt <- reads.dt[score > 0,]
  reads.dt$strand <- dt$strand[1]
  data.table::setnames(reads.dt, "score", "counts")
  data.table::fwrite(reads.dt, file = stringr::str_glue("{prefix}.cov.reads.csv"), row.names = FALSE)

  # Step 5: Calculate and write 3' stop-site counts (optimized with data.table)
  # Calculate 3' stops
  dt[, stops := data.table::fifelse(strand == "+", end, start)]

  # Counts for all reads
  all_stops <- dt[, .N, by = .(seqnames, start = stops, end = stops, strand)]
  data.table::setnames(all_stops, "N", "counts")
  # skip if high cov sites < 50
  if (sum(all_stops$counts >=3) < 50) return(0)
  data.table::fwrite(all_stops, file = stringr::str_glue("{prefix}.cov3.all.csv"), row.names = FALSE)

  # Counts for reads with "pa" tag
  pa_dt <- dt[!is.na(pa)]
  if (nrow(pa_dt) > 0){
    pa_stops <- pa_dt[, .N, by = .(seqnames, start = stops, end = stops, strand)]
    data.table::setnames(pa_stops, "N", "counts")
    data.table::fwrite(pa_stops, file = stringr::str_glue("{prefix}.cov3.pa.csv"), row.names = FALSE)
  }

  # Counts per cell barcode
  cb_stops <- dt[, .N, by = .(seqnames, start = stops, end = stops, strand, CB)]
  data.table::setnames(cb_stops, "N", "counts")
  data.table::fwrite(cb_stops, file = stringr::str_glue("{prefix}.cov3.cb.csv"), row.names = FALSE)

  return(reads)
}


#' @title Convert Data Frame Coordinates to New Coordinate System
#'
#' @param df A data frame with 'start' and 'end' columns.
#' @param newCoords A named vector， names are raw (characters) and values are new coordinates (as numeric).
#'
#' @return A data frame with 'start' and 'end' columns remapped to the new coordinate system, and 'width' recalculated.
#' @importFrom dplyr mutate filter
df_to_new_coords <- function(df, newCoords) {
  df.new <- df %>%
    dplyr::mutate(
      start = newCoords[as.character(start)],   # Map using newCoords
      end = newCoords[as.character(end)]
    ) %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::mutate(width = end - start + 1)
  return(df.new)
}

#' @title Calculate Flanking Coverage
#' @description This function extends a GenomicRanges object by a specified width,
#'   calculates the coverage of the expanded regions, and then converts the coverage
#'   into a data.table with one row per base pair.
#' @param gr A GenomicRanges object with a metadata column named 'counts' to be
#'   used as weights for coverage calculation.
#' @param ex An integer specifying the number of base pairs to extend on each side.
#'   The total width of the flanking region will be `2 * ex`.
#' @return A data.table containing the single-base-pair coverage data.
#' @importFrom GenomicRanges resize width coverage
flank_cov <- function(gr, ex = 5) {

  # Step 1: Extend the GenomicRanges to include flanking regions, centered on the original ranges.
  gr_extended <- GenomicRanges::resize(gr, fix = 'center', width = GenomicRanges::width(gr) + 2 * ex)

  # Step 2: Calculate the coverage of the extended ranges, using 'counts' as weights.
  cov.rle <- GenomicRanges::coverage(gr_extended, weight = "counts")

  # Step 3: Convert the Rle object into a single-base data.table.
  rle_to_dt(cov.rle)
}

#' @title Run Smoothing Splines
#' @description This function fits a smoothing spline to x and y coordinates,
#'   returning the smoothed y-values. It can either use cross-validation to
#'   select the smoothing parameter or use a user-defined parameter.
#' @param x A numeric vector representing the x-coordinates for the spline fit.
#' @param y A numeric vector representing the y-coordinates for the spline fit.
#' @param spar A numeric value specifying the smoothing parameter. If `cv` is
#'   FALSE, this value will be used. If `NULL`, a default of 0.7 will be applied.
#' @param cv A logical value. If `TRUE`, cross-validation is used to automatically
#'   select the smoothing parameter.
#' @return A numeric vector representing the smoothed y-values from the spline fit.
#' @export
#' @importFrom stats smooth.spline
#' @importFrom dplyr coalesce
run_splines <- function(x, y, spar = NULL, cv = FALSE) {
  if (cv){
    fit <- stats::smooth.spline(x, y, cv = TRUE)
  } else {
    spar <- dplyr::coalesce(spar, 0.7)
    fit <- stats::smooth.spline(x, y, spar = spar)
  }
  final <- predict(fit, x = x)$y
  return(final)
}
