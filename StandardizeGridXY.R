# ---- Downsample analyte spectrogram to the library grid ----
#
# functions:
#
# Uses bin-averaging (anti-alias) when exp grid is much finer than lib grid.
resample_to_lib_grid <- function(tt_exp, I_exp, grid_lib) {
  step_lib <- median(diff(grid_lib))
  step_exp <- median(diff(tt_exp))
  # If much finer, bin-average; else linear interpolation
  if (step_exp < 0.5 * step_lib) {
    # bin edges centered on lib grid
    edges <- c(grid_lib - step_lib/2, tail(grid_lib, 1) + step_lib/2)
    bin <- findInterval(tt_exp, edges, left.open = TRUE)  # 1..length(grid_lib)
    # fast mean per bin
    out <- rep(NA_real_, length(grid_lib))
    # avoid tapply order issues on empty bins
    ok <- bin >= 1 & bin <= length(grid_lib)
    if (any(ok)) {
      m <- tapply(I_exp[ok], bin[ok], mean, na.rm = TRUE)
      out[as.integer(names(m))] <- as.numeric(m)
    }
    # fill empty bins by linear interpolation over lib grid
    if (anyNA(out)) {
      out <- approx(x = grid_lib[!is.na(out)], y = out[!is.na(out)],
                    xout = grid_lib, rule = 2)$y
    }
    out
  } else {
    approx(tt_exp, I_exp, xout = grid_lib, rule = 2)$y
  }
}

# ---- Normalize analyte spectrogram to match library's scale (I_max = 10000) ----
normalize_to_Imax <- function(y, target_max = 10000) {
  m <- max(y, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(y)
  y * (target_max / m)
}

# ---- (Optional) tiny 2θ zero-shift correction after resampling ----
best_shift_correlation <- function(y_exp_on_grid, y_ref_on_grid, grid,
                                   max_shift = 0.20, step = 0.002) {
  shifts <- seq(-max_shift, max_shift, by = step)
  best <- -Inf; best_s <- 0; best_y <- y_exp_on_grid
  for (s in shifts) {
    y_s <- approx(x = grid + s, y = y_exp_on_grid, xout = grid, rule = 2)$y
    r <- suppressWarnings(cor(y_s, y_ref_on_grid))
    if (is.finite(r) && r > best) { best <- r; best_s <- s; best_y <- y_s }
  }
  list(cor = best, shift = best_s, y = best_y)
}
#
# run lines:
#
