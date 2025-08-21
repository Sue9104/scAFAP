#' @title Find Flexible Peaks in Smoothed Data
#' @description This function identifies flexible peaks in a smoothed data signal based on
#'   local trends, height, and width. It handles noise and merges adjacent regions.
#' @param pos A numeric vector of genomic positions.
#' @param values A numeric vector of smoothed signal values (e.g., from a spline fit).
#' @param min_cov The minimum smoothed value threshold for a region to be considered.
#' @param nups An integer specifying the number of consecutive upward trends to define a peak start.
#' @param ndowns An integer specifying the number of consecutive downward trends to define a peak end.
#' @param tolerance An integer allowing for some non-trend values within the `nups` or `ndowns` windows.
#' @param min_peak_height The minimum height for a peak to be considered.
#' @param min_peak_diff The minimum difference between the peak summit and its local valley.
#' @param min_peak_width The minimum width of a peak in base pairs.
#' @return A data.frame (tibble) with the identified peaks, including their height, summit,
#'   start, end, and width.
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile
#' @importFrom zoo rollapply
#' @importFrom tibble tibble
find_flexible_peaks <- function(pos, values,
                                min_cov = 3,
                                nups = 10, ndowns = 10, tolerance = 1,
                                min_peak_height = 3, min_peak_diff = 2,
                                min_peak_width = 100
) {

  # --- 1. Parameter validation ---
  if (!is.numeric(pos) || !is.numeric(values) || length(pos) != length(values)) {
    stop("`pos` and `values` must be same length")
  }
  empty_peaks <-
    data.frame(height = numeric(0), summit = numeric(0),
               start = numeric(0), end = numeric(0), width = numeric(0))

  # --- 2. Filter and merge adjacent regions ---
  min_cov_used <- stats::quantile(values, 0.1, na.rm = TRUE)
  if (stats::quantile(values, 0.99, na.rm = TRUE) > 20) {
    min_cov_used <- max(min_cov_used, min_cov)
  }
  is_valid <- values >= min_cov_used
  if (!any(is_valid)) return(empty_peaks)
  valid_rle <- rle(is_valid)
  points <- cumsum(c(1, valid_rle$lengths))
  valid_regions <- lapply(which(valid_rle$values), function(i){
    c(points[i], points[i] + valid_rle$lengths[i] - 1)
  })
  # merge regions within 50bp
  if (length(valid_regions) > 1) {
    merged_regions <- list()
    current_region <- valid_regions[[1]]
    for (j in 2:length(valid_regions)) {
      next_region <- valid_regions[[j]]
      if ((pos[next_region[1]] - pos[current_region[2]] - 1) < 50) {
        current_region[2] <- next_region[2]
      } else {
        merged_regions[[length(merged_regions) + 1]] <- current_region
        current_region <- next_region
      }
    }
    merged_regions[[length(merged_regions) + 1]] <- current_region
    valid_regions <- merged_regions
  }
  if (length(valid_regions) == 0) return(empty_peaks)

  # --- 3. diff map to +1, -1 and use rollapply to find trends ---
  min_peak_diff <- max(min_peak_diff, stats::quantile(values, 0.25) - min(values))
  peaks <- empty_peaks
  for (region_id in seq_along(valid_regions)) {
    region_start <- valid_regions[[region_id]][1]
    region_end <- valid_regions[[region_id]][2]
    # make sure the region is long enough to hold a peak pattern
    if ((region_end - region_start + 1) < (max(nups, ndowns) + 1)) next

    segment_values <- values[region_start:region_end]
    trend_signs <- sign(diff(segment_values))
    # --- vectorized search for upward trends ---
    up_windows <- zoo::rollapply(trend_signs, width = nups, FUN = function(x) sum(x == 1), align = "left", fill = 0)
    up_rle <- rle(up_windows >= (nups - tolerance))
    ups <- cumsum(c(1,up_rle$lengths))[which(up_rle$values)] + region_start - 1 + round(nups/2)

    # --- vectorized search for downward trends ---
    down_windows <- zoo::rollapply(trend_signs, width = ndowns, FUN = function(x) sum(x == -1), align = "left", fill = 0)
    down_rle <- rle( down_windows >= (ndowns - tolerance) )
    downs <- cumsum(down_rle$lengths)[which(down_rle$values)] + region_start - 1 + round(ndowns/2)

    # --- combine upward and downward trends to form peaks ---
    if (length(ups) > 0 && length(downs) > 0) {
      for (id in seq_along(ups)){
        s <- ups[id]
        # filter by peak width
        is_e <- pos[downs] - pos[s] > min_peak_width
        if (! any(is_e)) next
        e <- downs[min(which(is_e))]
        # is the nearest to e
        if ( (id < length(ups)) & (ups[(id+1)] < e) ) next
        ids <- which((pos >= pos[s]) & (pos <= pos[e]) )
        if (length(ids) < 50) next
        max_height <- max(values[ids])
        min_height <- min(values[ids])
        summit_id <- which(values == max_height)
        # shores <- max(values[c(s,e)])
        if ( (max_height > min_peak_height) & (max_height - min_height >= min_peak_diff)){
          peak <- tibble::tibble(height = round(max_height, 3),
                                 summit = pos[summit_id],
                                 start = pos[s], end = pos[e],
                                 width = end - start + 1)
          peaks <- dplyr::bind_rows(peaks, peak)
        }
      }
    }

    # Continuous upward trend
    if ((region_id == length(valid_regions)) & all(down_rle$values == F) & (length(ups) == 1)) {
      up_width <- up_rle$lengths[up_rle$values == T]
      up_ratio <- up_width / sum(up_rle$lengths)
      if ((up_width >= min_peak_width) & (up_ratio >= 0.9)){
        s <- ups
        e <- ups + up_width
        min_height <- values[s]
        max_height <- values[e]
        summit_id <- max(e - 100, s)
        if ( (max_height > min_peak_height) & (max_height - min_height >= min_peak_diff)){
          peak <- tibble::tibble(height = round(max_height, 3),
                                 summit = pos[summit_id],
                                 start = pos[s], end = pos[e],
                                 width = end - start + 1)
          peaks <- dplyr::bind_rows(peaks, peak)
        }
      }
    }
    # Continuous downward trend
    if ((region_id == 1) & all(up_rle$values == F) & (length(downs) == 1)) {
      down_width <- down_rle$lengths[down_rle$values == T]
      down_ratio <- down_width / sum(down_rle$lengths)
      if ((down_width >= min_peak_width) & (down_ratio >= 0.9)){
        s <- downs
        e <- downs + down_width
        min_height <- values[e]
        max_height <- values[s]
        summit_id <- min(s + 100, e)
        if ( (max_height > min_peak_height) & (max_height - min_height >= min_peak_diff)){
          peak <- tibble::tibble(height = round(max_height, 3),
                                 summit = pos[summit_id],
                                 start = pos[s], end = pos[e],
                                 width = end - start + 1)
          peaks <- dplyr::bind_rows(peaks, peak)
        }
      }
    }
  }
  return(peaks)
}
