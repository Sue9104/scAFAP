library(argparser)
p <- arg_parser("Description")
p <- add_argument(p, "annodir", help="annotation directory", type="character")
p <- add_argument(p, "inbam", help="input bam", type="character")
p <- add_argument(p, "sample", help="sample name", type="character")
p <- add_argument(p, "outdir", help="output directory", type="character")
p <- add_argument(p, "--cores", default = 40, help="parallel cores", type = "integer")
argv <- parse_args(p)


# anno_dir = "~/projects/seq3/20250514/scPAPA/hg38"
# inbam = "~/projects/seq3/20250514/1_matrix/NC-1/outs/possorted_genome_bam.bam"
# sample = "NC-1"
# outdir = "~/projects/seq3/20250514/scPAPA/human"
# cores = 4
# argv <- list("annodir" = anno_dir, "inbam" = inbam, "sample" = sample, "outdir" = outdir, "cores" = cores)

# # test single gene
# adaptive_detect_pas('JAK1', inbam, outdir, sample, anno_dir)

library(scAFAP)
genes <- read.table(stringr::str_glue('{argv$annodir}/genes.list'), header = F)$V1
BiocParallel::bplapply(
  genes,
  BPPARAM = BiocParallel::SnowParam(workers = argv$cores, tasks = 500, stop.on.error = F),
  adaptive_detect_pas,
  inbam = argv$inbam,
  outdir = argv$outdir,
  sample = argv$sample,
  anno_dir = argv$annodir
)


merge_create_mtx(argv$outdir, argv$sample)
