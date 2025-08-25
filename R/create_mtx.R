#' @title Merge and Convert APA Data to 10x Format
#' @description This function merges gene-level APA quantification files and converts the final expression matrix into the standard 10x Genomics format.
#' @param genedir Directory containing the gene-level output files.
#' @param sample Sample name, used to find and name files.
#' @param overwrite Logical value, whether to overwrite existing 10x output files. Defaults to TRUE.
#' @return This function invisibly returns NULL and saves the 10x format files to the specified output directory.
#' @export
#' @importFrom data.table fread rbindlist
#' @importFrom DropletUtils write10xCounts
#' @importFrom Matrix Matrix
#' @importFrom stringr str_glue str_split_i
merge_create_mtx <- function(outdir, sample, overwrite = TRUE) {
  suppressPackageStartupMessages({library(dplyr)})
  genes <- read.table(stringr::str_glue('{argv$outdir}/{argv$sample}.todo'))$V1

  # Retrieve all matrix files
  matrix_files <- Sys.glob(str_glue("{outdir}/tmps/{genes}.{sample}.pas.matrix.csv"))
  raw_table <- data.table::rbindlist(lapply(matrix_files, data.table::fread), fill = TRUE)
  raw_table[is.na(raw_table)] <- 0

  # Extract PAS names and cell barcodes
  pas_names <- raw_table[[1]]
  cell_names <- names(raw_table)[-1]

  # Convert data to a sparse matrix, which is more memory efficient than a dense matrix
  sparse_matrix <- Matrix::Matrix(as.matrix(raw_table[, -1]), sparse = TRUE)

  # Set row and column names
  rownames(sparse_matrix) <- pas_names
  colnames(sparse_matrix) <- cell_names

  # Create output directory
  output_dir <- stringr::str_glue("{outdir}/matrix/{sample}")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save the sparse matrix in 10x format
  DropletUtils::write10xCounts(
    path = output_dir,
    overwrite = overwrite,
    x = sparse_matrix,
    barcodes = colnames(sparse_matrix),
    gene.id = stringr::str_split_i(pas_names, ':', 1),
    gene.symbol = pas_names,
    version = "3"
  )

  files <- Sys.glob(stringr::str_glue("{outdir}/tmps/{genes}.{sample}.pas.stats.csv"))
  pas.df <-
    purrr::map_dfr(split(files, ceiling(seq_along(files) / 20)),
                   function(x) {
                     readr::read_csv(x, col_types = readr::cols(.default = "c"))
                   })
  write.csv(pas.df, stringr::str_glue("{outdir}/{sample}.pas.stats.csv"), quote = F, row.names = F)
}
