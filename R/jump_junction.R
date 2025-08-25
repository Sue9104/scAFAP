#' @title Prepare Junction Data for Bezier Plotting
#'
#' @description Calculates midpoints and heights for plotting junctions as
#'   Bezier curves, typically for sashimi plots. The height is based on the
#'   log2 of junction counts.
#'
#' @param junctions A data frame containing junction information wit
#'   'start', 'end', and 'counts' columns.
#'
#' @return A tibble with 'x', 'y', 'group', and 'counts' columns,
#'   suitable for `ggplot2::geom_bezier`.
#' @importFrom tibble tibble
#' @importFrom dplyr mutate rowwise do bind_rows n
junc_to_bezier <- function(junctions){
  bezier_data <- tibble(x=1, y =1, group=1, counts=0)
  if (nrow(junctions) > 0){
    bezier_data <-
      junctions %>%
      dplyr::mutate(
        mid = (start + end) / 2,
        height = log2(counts + 1)
      ) %>%
      dplyr::mutate(group = dplyr::row_number()) %>% # Assign a unique group for each junction
      dplyr::rowwise() %>% # Process row by row
      dplyr::do({
        data.frame(
          x = c(.$start, .$mid, .$end),     # X-axis: start → midpoint → end
          y = c(1, .$height, 1),             # Y-axis: bottom → arc height → bottom
          group = .$group,
          counts = .$counts
        )
      }) %>%
      dplyr::ungroup() %>% # Remove rowwise grouping after do()
      dplyr::bind_rows()
  }
  return(bezier_data)
}

#' @title Generate a Sashimi Plot for Junctions
#'
#' @description Creates a ggplot2-based sashimi plot visualizing gene junctions,
#'   exons, and UTRs. Junctions are plotted as Bezier curves, with counts
#'   labeled.
#'
#' @param gene.name Character string, the name of the gene of interest.
#' @param gene.strand Character string, the strand of the gene ("+" or "-").
#' @param junctions A data frame containing junction information with
#'   'start', 'end', and 'counts' columns.
#' @param exon.gr A GRanges object representing exon regions.
#' @param utr.gr A GRanges object representing UTR regions.
#' @param xtick Numeric, interval for x-axis ticks.
#'
#' @return A ggplot object representing the sashimi plot.
#' @importFrom ggplot2 ggplot geom_label scale_size_area geom_rect
#'   coord_cartesian scale_x_continuous theme_minimal labs theme element_blank
#'   element_text
#' @importFrom tibble as_tibble
#' @importFrom ggforce geom_bezier
#' @importFrom logger log_warn
junc_sashimi_plot <- function(gene.name, junctions, exon.gr, utr.gr, xtick = 1000) {
  if (nrow(junctions) == 0) {
    logger::log_warn("No available junctions for sashimi plot.")
    return(ggplot2::ggplot())
  }
  gene.strand <- as_tibble(exon.gr)$strand[1]
  bezier_data <- junc_to_bezier(junctions)
  label_data <- dplyr::filter(bezier_data, .data$y != 1)

  p.junc <-
    ggplot2::ggplot() +
    ggforce::geom_bezier(
      data = bezier_data,
      ggplot2::aes(x = x, y = y, group = group, linewidth = counts),
      color = "brown2"
    ) +
    ggplot2::geom_label(
      data = label_data,
      ggplot2::aes(x = x, y = y / 2, label = counts),
      size = 4
    ) +
    # Use scale_size_area to map line width proportionally to counts
    ggplot2::scale_size_area(max_size = 1.5, guide = "none") + # Ensure `linewidth` is mapped properly
    ggplot2::geom_rect(
      data = tibble::as_tibble(exon.gr),
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -0.4, ymax = 1),
      fill = "gray"
    ) +
    ggplot2::geom_rect(
      data = tibble::as_tibble(utr.gr),
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
      fill = "royalblue"
    ) +
    ggplot2::coord_cartesian(xlim = c(min(bezier_data$x) , max(bezier_data$x) )) +
    ggplot2::scale_x_continuous(
      breaks = seq(
        floor(min(bezier_data$x) / xtick) * xtick, # Align breaks to xtick multiples
        ceiling(max(bezier_data$x) / xtick) * xtick,
        by = xtick
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = paste(gene.name, gene.strand),
      x = "Genomic Position",
      y = ""
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = .5)
    )
  return(p.junc)
}


#' @title Plot Coverage with Junction Arcs
#'
#' @description Generates a ggplot2 plot showing read coverage as points/smooth line,
#'   along with junction arcs. Exons, UTRs, and highlight regions can also be displayed.
#'
#' @param cov.df A data frame with coverage data, including 'start', 'end', and 'counts'.
#' @param junctions A data frame with junction data, including 'start', 'end', and 'counts'.
#' @param gene.gr A GRanges object for the gene, used for plot boundaries and title.
#' @param exon.gr A GRanges object representing exon regions.
#' @param utr.gr A GRanges object representing UTR regions.
#' @param highlight.gr A GRanges object for regions to highlight (optional).
#' @param add_to_title Character string to append to the plot title (optional).
#'
#' @return A ggplot object visualizing coverage and junctions.
#' @importFrom ggplot2 ggplot aes facet_wrap geom_point geom_smooth geom_rect
#'   geom_label coord_cartesian labs theme_minimal element_text element_rect
#'   scale_size_area
#' @importFrom dplyr filter mutate if_else
#' @importFrom stringr str_glue
#' @importFrom tibble as_tibble
#' @importFrom ggforce geom_bezier
plot_cov_with_junc <- function(cov.df, junctions, gene.gr, exon.gr, utr.gr,
                               highlight.gr = GenomicRanges::GRanges(),
                               add_to_title = "") {
  if (nrow(cov.df) == 0){
    logger::log_debug("No coverage data available for coverage plot.")
    return(ggplot2::ggplot())
  }
  if (!("reads" %in% names(cov.df))) {
    cov.df <- dplyr::mutate(cov.df, reads = "all")
  }

  # Define plot boundary based on gene.gr
  if ( length(highlight.gr) > 0){
    pmin <- max(GenomicRanges::start(highlight.gr) - 100, 0)
    pmax <- GenomicRanges::end(highlight.gr) + 100
  } else {
    pmin <- GenomicRanges::start(gene.gr)
    pmax <- GenomicRanges::end(gene.gr)
  }

  # Filter and preprocess coverage within gene boundary
  cov.df <-
    dplyr::filter(cov.df, start >= pmin, end <= pmax) %>%
    dplyr::mutate(mid = round((start + end) / 2))
  if (nrow(cov.df) == 0){
    logger::log_debug("No coverage data within gene boundaries for plot.")
    return(ggplot2::ggplot())
  }

  # Create bezier arcs for junctions
  juncBezier <- junc_to_bezier(junctions)

  # Adjust junction arc heights to fit coverage plot scale
  zoom <- 1 # Default zoom if no junctions
  if(nrow(juncBezier) > 0 && max(juncBezier$y) > 0) {
    zoom <- max(cov.df$counts) / max(juncBezier$y)
  }
  zoom2 <- max(cov.df$counts) / 20 # Constant for exon/UTR height

  if(nrow(juncBezier) > 0) {
    juncBezier <- juncBezier %>%
      dplyr::mutate(y = dplyr::if_else(.data$y == 1, -1, -zoom * .data$y), reads = "all") # Invert and scale
  }

  juncLabels <- dplyr::filter(juncBezier, .data$y != -1) %>%
    dplyr::mutate(y = .data$y / 2)

  # Prepare exon/UTR data for plotting beneath coverage
  exon_plot_data <- tibble::as_tibble(exon.gr)
  utr_plot_data <- tibble::as_tibble(utr.gr)

  # Plot
  p <- ggplot2::ggplot(cov.df, ggplot2::aes(x = mid, y = counts, color = reads)) +
    ggplot2::facet_wrap(~reads, ncol = 1, scales = "free_y") +
    ggplot2::geom_point(size = 0.5) +
    # ggplot2::geom_smooth(span = 0.2,  color = "gray50", se = FALSE) +
    ggplot2::geom_rect(data = tibble::as_tibble(highlight.gr),
                       ggplot2::aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = "grey", alpha = 0.4) +
    ggforce::geom_bezier(data = juncBezier,
                         ggplot2::aes(x = x, y = y, group = group, linewidth = counts),
                         color = "brown2") +
    ggplot2::scale_size_area(max_size = 3, guide = "none") +
    ggplot2::geom_label(data = juncLabels,
                        ggplot2::aes(x = x, y = y, label = counts), size = 3) +
    ggplot2::geom_rect(data = exon_plot_data,
                       ggplot2::aes(xmin = start, xmax = end, ymin = -zoom2, ymax = zoom2),
                       inherit.aes = FALSE,
                       color = "darkviolet", fill = "darkviolet") +
    ggplot2::geom_rect(data = utr_plot_data,
                       ggplot2::aes(xmin = start, xmax = end, ymin = -zoom2 * 0.8, ymax = zoom2 * 0.8),
                       inherit.aes = FALSE,
                       color = "palegreen4", fill = "palegreen4") +
    ggplot2::coord_cartesian(xlim = c(pmin, pmax)) +
    ggplot2::labs(x = "Genomic Position", y = "Counts", title = add_to_title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.grid = ggplot2::element_blank())

  return(p)
}

#' @title Analyze Gene for Alternative Splicing (Skipping Junctions)
#'
#' @description This function identifies "skipped" regions in a gene based on
#'   junction and coverage data, starting the analysis from the 3' end. It then
#'   transforms genomic coordinates to a new continuous space and generates
#'   sashimi and coverage plots both before and after the coordinate transformation.
#'
#' @param gene.gr A GRanges object representing gene regions
#' @param exon.gr A GRanges object representing exon regions within the gene.
#' @param utr.gr A GRanges object representing UTR regions within the gene.
#' @param junc.file RSeQC junction file
#' @param read_cov.file Coverage file for reads
#' @param minJuncCount Numeric, minimum read count for a junction to be considered.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{jumpRegion}: GRanges object of identified jumped regions.
#'     \item \code{newCoord}: Named numeric vector mapping original genomic
#'       coordinates to the new continuous coordinate space.
#'   }
#'
#' @details
#' The "jumping" junction identification logic is as follows:
#' \enumerate{
#'   \item Iterate from the 3' end of the gene towards the 5' end.
#'   \item For each pivot point (junction start/end depending on strand),
#'     sum the counts of all junctions originating/ending at that point.
#'   \item Compare this sum to the "unspliced" coverage at that exact point.
#'     Unspliced coverage is defined as total coverage at the point minus
#'     the sum of junction counts.
#'   \item If the total junction counts are greater than the unspliced counts,
#'     a "jumping" event is indicated.
#'   \item If multiple overlapping junctions satisfy the condition, the one
#'     with the highest count is chosen.
#'   \item Once a jumping region is identified, subsequent pivot points
#'     within that region are skipped to avoid redundant processing.
#' }
#'
#' @import GenomicRanges
#' @importFrom purrr map set_names
#' @importFrom stringr str_glue
#' @importFrom logger log_info log_warn
#' @importFrom cowplot plot_grid
#' @importFrom tibble as_tibble
#' @export
jump_splicing_region <- function(gene_name, gene.gr, utr.gr, exon.gr,
                                 junc.file, read_cov.file,
                                 minJuncCount = 20) {

  # Step 0: Input Validation and Preprocessing ------------------------------------------
  {
    # incase pseudogene have two regions
    gene.gr <- gene.gr[1]
    gene_strand <- as.character(GenomicRanges::strand(gene.gr)[1])
    chrom <- as.character(GenomicRanges::seqnames(gene.gr)[1])
    logger::log_info(stringr::str_glue("Starting analysis for gene: {gene_name}"))
    res <- list(
      "jumpRegion" = GRanges(),
      "newCoord" = set_names(1:width(gene.gr),
                             nm = as.character(as_tibble(gene.gr)$start:as_tibble(gene.gr)$end))
    )

    # read coverages and junctions
    {
      cov.df <- read.csv(read_cov.file)
      cov_rle <-
        makeGRangesFromDataFrame(cov.df, keep.extra.columns = T) %>%
        coverage(weight = 'counts') %>%
        unlist()
      if (!file.exists(junc.file) | file.size(junc.file) < 70){
        logger::log_warn("NoneExisted or NoJunctions for ", junc.file)
        return(res)
      }
      junctions <-
        readr::read_delim(junc.file,
                          col_names = c("seqnames", "start", "end", "counts", rep(NULL, Inf)),
                          delim = "\t", col_select = 1:4, skip = 1, show_col_types = F) %>%
        mutate(start = start + 1, width = end - start + 1,
               strand = gene_strand) %>%
        filter(counts >= minJuncCount)
      junc.gr <- makeGRangesFromDataFrame(junctions, keep.extra.columns = T)
      if (length(junc.gr) == 0) return(res)
      junc_rle <- coverage(junc.gr, weight = "counts") %>% unlist()

    }
  }

  # Step 1: Identify Skipped Regions from 3' to 5' ----------------------------------------------
  {
    jumped_regions <- GenomicRanges::GRanges()
    if (gene_strand == "+") {
      pivot_points <- sort(unique(GenomicRanges::end(junc.gr)), decreasing = TRUE)
    } else {
      pivot_points <- sort(unique(GenomicRanges::start(junc.gr)), decreasing = FALSE)
    }

    skipped_pivot_coords <- GenomicRanges::GRanges()
    for (i in 1:length(pivot_points)){
      piv_point <- pivot_points[i]
      # Skip if this pivot point is within an already identified jumped region
      piv_point.gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(piv_point, piv_point), gene_strand)
      if (countOverlaps(piv_point.gr, skipped_pivot_coords) >0) next
      spliced <- decode(junc_rle[piv_point])
      total <- decode(cov_rle[piv_point])
      unspliced <- total - spliced
      logger::log_debug("pivot: ", piv_point, " spliced: ", spliced, " unspliced: ", unspliced)

      # If junction reads are dominant, this indicates a potential jumping event.
      if (spliced > unspliced) {
        if (gene_strand == "+") {
          relevant_juncs <- junc.gr[GenomicRanges::end(junc.gr) == piv_point]
        } else { # gene_strand == "-"
          relevant_juncs <- junc.gr[GenomicRanges::start(junc.gr) == piv_point]
        }
        # Condition 2: If multiple relevant junctions exist, select the one with the highest count.
        selected_junc <- relevant_juncs[which.max(GenomicRanges::mcols(relevant_juncs)$counts)]
        skipped_pivot_coords <- c(skipped_pivot_coords, selected_junc)
      }
    }
    jumped_regions <- unlist(GenomicRanges::GRangesList(skipped_pivot_coords))
  }


  # Step 2: Coordinate Transformation ---------------------------------
  {
    retained_regions <- GenomicRanges::setdiff(gene.gr, jumped_regions)
    retained_regions <- sort(retained_regions)
    coords <-
      purrr::map(seq_along(retained_regions),
                 ~ GenomicRanges::start(retained_regions[.x]):GenomicRanges::end(retained_regions[.x])) %>%
      unlist()
    coord2new <- purrr::set_names(seq_along(coords), nm = as.character(coords))
  }


  # Step 3: Generate Plots (Before and After Transformation)  --------------------------------------------------------
  {
    # coverage
    {
      # Prepare coverage data for plotting (add 'reads' column for facet)
      coverage_data_before_transform <- dplyr::mutate(cov.df, reads = "all")
      p.cov_junc.before <- plot_cov_with_junc(
        coverage_data_before_transform,
        junctions,
        gene.gr,
        exon.gr,
        utr.gr,
        add_to_title = str_glue("{gene_name} {gene_strand} (Raw)")
      )
      x <- ggplot_build(p.cov_junc.before)$data[[1]]$x
      junc_data <- tibble(pos = x, spliced = if_else(pos <= length(junc_rle),pos, 1))
      junc_data$spliced <- decode(junc_rle[junc_data$spliced])
      p.cov_junc.before <-
        p.cov_junc.before +
        geom_point(data = junc_data, aes(x=pos,y=spliced), inherit.aes = F, size = 0.1, color = "blue")
    }
    # Sashimi Plots
    {
      xtick_sashimi <- floor(width(gene.gr) / 10)
      p.junc.before <- junc_sashimi_plot(gene_name, junctions, exon.gr, utr.gr, xtick = xtick_sashimi)
      p.junc.after <- junc_sashimi_plot(gene_name, as_tibble(jumped_regions), exon.gr, utr.gr, xtick = xtick_sashimi)
      p.cov_junc.before.junc <- p.cov_junc.before + coord_cartesian(xlim = range(ggplot_build(p.junc.before)$data[[1]]$x))
      plot_junc <- cowplot::plot_grid(p.cov_junc.before.junc, p.junc.before, p.junc.after,
                                      ncol = 1, labels = c("Cov", "All", "Skip"), align = "v")
      # print(plot_junc)
    }
  }


  # Step 4: Return Results  ----------------------
  res <- list(
    "jumpRegion" = jumped_regions,
    "newCoord" = coord2new,
    "plot" = plot_junc
  )
  return(res)
}
