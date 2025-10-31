#
#####===== Matching techniques for preprocessed, experimental .xy spectra =====####               
#                                                                                
                                                                                  # for use with an audited, standardized .xy reference library
###--- Default parameter settings ---###
#
#
# matching parameters (choice 2):      current settings:       suggested values:
       DEFAULT_LIB_STEP                <- 0.02                 # 0.02             # deg 2θ (your library step)
       DEFAULT_MAX_SHIFT               <- 0.20                 # 0.20             # deg 2θ search window (±)
       DEFAULT_SHIFT_STEP              <- 0.02                 # 0.02             # deg 2θ shift increment
       DEFAULT_TOP_N                   <- 10                   # 10               # how many matches to list
       DEFAULT_OVERLAY_K               <- 3L                   # 3L               # how many overlays to plot/save
       DEFAULT_AMORPH                  <- FALSE                # FALSE            # library 'amorphous' transformation option
       DEFAULT_AMORPH_FWHM             <- 0.20                 # 0.20             # moderate environmental lattice strain (for amorphous transform)
       out_dir                         <- NULL
       out_csv                         <- NULL

       DEFAULT_SCAN <- list(                                                      # global scan default parameters
         max_shift                     = 1.0,                  # 1.0               
         shift_step                    = 0.2,                  # 0.1
         refine                        = TRUE,                 # TRUE
         refine_win                    = 0.06,                 # 0.06
         refine_step                   = 0.02                  # 0.01
       )
       DEFAULT_CFG <- list(                                                       # peak weighting/scoring parameters
         tol                           = 0.75,                 # 0.08             # lateral tolerance (deg 2θ)
         heavy_thr                     = 0.30,                 # 0.50             # peak intensity threshhold for 'heavy' peak weight
         heavy_penalty                 = 1.0,                  # 1.0              # penalty for absent heavy peaks
         weight_gamma                  = 1.2,                  # 1.0              # exponent to modify library peaks
         min_matches                   = 0L,                   # 0L               # optional gate
         require_heavy_match           = FALSE                 # FALSE            # TRUE = "oh shit" lever (if too many/distracting matches)
       )
       DEFAULT_DET_EXP <- list(
         rel_min                       = 0.15,                 # 0.20             # minimum peak height
         min_spacing                   = 0.25,                 # 0.01             # allowed peak spacing
         prominence_rel                = 0.08,                 # 0.02             # accept peaks over this percentage of the max
         prom_window                   = 0.25,                 # 0.08             # half-window for local baseline
         smooth_k                      = 3,                    # 3                # MA window for detection only (odd)
         w_gamma                       = 0.85,                 # 0.85             # compress dynamic range of weights
         min_keep                      = 7                     # 8                # minimum amount of peaks to return (override rules)
       )
       options(warnPartialMatchDollar  = TRUE)
#
#
###--- top menu ---###
#
#
run_matching_menu <- function() {
              choice <- utils::menu(c(
                          "Peak-only (fast) — use indexed peaks ≥20% of max; shift ±1.0° step 0.1°",
                          "Full-profile (correlation) — clamped scaling; optional amorphization"),
                          title = "Choose matching method")
                        if (choice == 0) { 
                           message("Cancelled.")
                           return(invisible(NULL)) 
                        }
                        message("Select prepped experimental .xy") 
                        exp_xy <- try(file.choose(), silent = TRUE)                # Pick experimental file (prepped .xy) *after* method choice
                        if (inherits(exp_xy, "try-error") 
                           || !nzchar(exp_xy) || !file.exists(exp_xy)) {
                           message("No file selected.")
                           return(invisible(NULL))
                        }
             lib_dir <- if (exists("get_or_choose_lib_dir")) {                     # Pick/confirm library directory
                           get_or_choose_lib_dir(default = dirname(exp_xy))
                        } else { if (.Platform$OS.type == "windows") {
                                    if (is.na(d)) d <- ""
                                 } else {d <- readline(sprintf(
                                     "Enter library directory path [%s]: ", 
                                     dirname(exp_xy)))
                                     if (!nzchar(d)) d <- dirname(exp_xy)
                                 }
                                 if (!dir.exists(d)) { message("No valid library directory.")
                                    return(invisible(NULL))
                                 }
                                 normalizePath(d, winslash = "/")
                        }
                        message("Library directory: ", normalizePath(lib_dir, winslash = "/", 
                                mustWork = FALSE))
                        if (choice == 1) {                                         # Peak-only matcher (stick method)
                           compare_peaks_to_library(exp_xy_path = exp_xy,
                                                    lib_dir = lib_dir)
                        } else if (choice == 2) {                                  # Full-profile pattern matcher (3 correlation
                    use_amorph <- tolower(trimws(readline(
                                    "Turn ON amorphization for library? [y/N]: "))) 
                                    %in% c("y","yes")
                   amorph_fwhm <- 0
                                  if (use_amorph) {
                               ai <- suppressWarnings(as.numeric(readline(
                                     "Gaussian FWHM (deg 2θ) [default 0.20]: ")))
                                  if (!is.finite(ai)) ai <- 0.20 
                   amorph_fwhm <- ai
                                  }
                                  compare_xy_to_library(exp_xy_path = exp_xy, lib_dir = lib_dir)
                        }
                      invisible(TRUE)
                    }
#
#
###--- Helper functions ---###
#
#
resample_linear <- function(tt, I, grid) {
                     approx(tt, I, xout = grid, rule = 2, ties = mean)$y
                   }
detect_step <- function(tt) median(diff(tt), na.rm = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a                                 # tiny `%||%` helper
validate_list <- function(x, defaults) {
                   if (is.null(defaults) || !is.list(defaults))
                   stop("defaults must be a list", call. = FALSE)
                   if (is.null(x)) return(defaults)
                   if (!is.list(x)) x <- as.list(x)                               # allow named atomic vectors or lists
                   if (is.null(names(x)) || any(names(x) == "")) {                # only keep named entries
                 x <- x[!is.null(names(x)) & names(x) != ""]
                   }
                   modifyList(defaults, x)
                 }
# Plot to PNG (pattern match)
overlay_plot <- function(grid, x, y_shifted_scaled, title, diff = TRUE, file = NULL) {
                   if (!is.null(file)) png(file, width = 1400, height = 900, res = 150)
             op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
                   layout(matrix(c(1, if (diff) 2 else 1), nrow = if (diff) 2 else 1, byrow = TRUE),
                          heights = if (diff) c(2,1) else 1)
                   par(mar = c(4,4,2,1))
             yl <- range(x, y_shifted_scaled, finite = TRUE)
                   plot(grid, x, type = "l", col = 1, ylim = yl,
                        xlab = expression(2*theta*" (deg)"), ylab = "Intensity", main = title)
                   lines(grid, y_shifted_scaled, col = 2)
                   legend("topright", bty = "n", lty = 1, col = c(1,2), 
                          legend = c("experimental", "library (shifted, scaled)"))
                   if (diff) {
                      par(mar = c(4,4,1,1))
                      plot(grid, x - y_shifted_scaled, type = "l",
                           xlab = expression(2*theta*" (deg)"), ylab = "Difference")
                      abline(h = 0, lty = 3)
                   }
                   if (!is.null(file)) dev.off()
                }
# Plot to PNG (peak match)
stick_plot <- function(exp_peaks, lib_peaks_shifted, title) {                
           op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
                 par(mar = c(4,4,2,1))
           xr <- range(c(exp_peaks$pos, lib_peaks_shifted$pos))
                 plot(xr, c(0, 1.05), type = "n",
                      xlab = expression(2*theta*" (deg)"), 
                      ylab = "Relative peak height", main = title)
                 segments(exp_peaks$pos, 0, exp_peaks$pos, exp_peaks$w, col = 1, lwd = 2)
                 segments(lib_peaks_shifted$pos, 0, lib_peaks_shifted$pos, 
                          lib_peaks_shifted$w, col = 2, lwd = 2)
                 legend("topright", bty = "n", lty = 1, col = c(1,2), 
                        legend = c("experimental peaks", "library peaks (shifted)"))
              }
# Atomic CSV and PNG writer
make_unique_path <- function(path) {
               ext  <- if (grepl("\\.[^.]+$", path)) sub(".*(\\.[^.]+)$", "\\1", path) else ""
               base <- if (nzchar(ext)) sub(paste0(ext, "$"), "", path) else path
                       for (i in 1:999) {
                           cand <- paste0(base, "-new-", sprintf("%02d", i), ext)
                           if (!file.exists(cand)) return(cand)
                       }
                       stop("Could not create a unique filename for: ", path)
                    }
safe_write_csv <- function(filename, df, ...) {
                     dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
              tmp <- paste0(filename, ".tmp")
               ok <- tryCatch({utils::write.csv(df, tmp, row.names = FALSE, ...)
                               TRUE}, 
                               error = function(e) { message("CSV write failed: ", 
                               conditionMessage(e)); FALSE 
                               }
                     )
                     if (!ok) return(invisible(NULL))                              # Try to finalize atomically; if locked, fall back
                     if (!file.rename(tmp, filename)) {
                 alt <- make_unique_path(filename)
                        if (!file.rename(tmp, alt)) {message(
                           "CSV finalize failed (file may be locked). Temp kept at: ", 
                           tmp) return(invisible(tmp))
                        } else { message("Target CSV busy; saved as: ", alt)
                            return(invisible(alt))
                        }
                     } else { if (file.size(filename) <= 0) message(                  # sanity check
                                 "Warning: CSV saved but size is 0 bytes: ", filename)
                              return(invisible(filename))
                     }
                  }
safe_write_png <- function(filename, plotter, width = 1400, height = 900, res = 150, bg = "white") {
                    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
             tmp <- paste0(filename, ".tmp.png")                            # explicit .png to satisfy some viewers
              ok <- tryCatch( {grDevices::png(tmp, width = width, height = height, 
                                              res = res, bg = bg)
                                  on.exit(try(grDevices::dev.off(), silent = TRUE), 
                                          add = TRUE) plotter() TRUE}, 
                                  error = function(e) { message("PNG write failed: ", 
                                             conditionMessage(e))
                                             FALSE 
                                          } 
                            )
                    if (!ok) return(invisible(NULL))
                    if (!file.rename(tmp, filename)) {
                alt <- make_unique_path(filename)
                       if (!file.rename(tmp, alt)) {message(
                          "PNG finalize failed (file may be locked). Temp kept at: ", tmp)
                          return(invisible(tmp))
                       } else {
                           message("Target PNG busy; saved as: ", alt)
                           return(invisible(alt))
                       }
                    } else { if (file.size(filename) <= 0) message(
                                "Warning: PNG saved but size is 0 bytes: ", filename)
                                return(invisible(filename))
                    }
                  }
# If a usable global lib_dir exists, prompt to reuse or rewrite
get_or_choose_lib_dir <- function(default = NULL) {
                    prev <- tryCatch(get("lib_dir", envir = .GlobalEnv), error = function(e) NULL)
                            if (is.character(prev) && length(prev) == 1L && nzchar(prev) && dir.exists(prev)) {
                        ans <- tolower(trimws(readline(sprintf(
                                 "Reuse existing library directory?\n  %s\n[Y/n]: ",
                                 normalizePath(prev, winslash = "/")))))
                               if (ans %in% c("", "y", "yes")) {
                                  message("Using library directory: ", normalizePath(prev, winslash = "/"))
                                  return(prev)
                               }
                            }
                      sel <- ""
                            if (requireNamespace("rstudioapi", quietly = TRUE) 
                               && rstudioapi::isAvailable()) {
                        sel <- tryCatch(rstudioapi::selectDirectory(path = 
                                  default %||% getwd()), error = function(e) "")
                            } else if (.Platform$OS.type == "windows") {
                               sel <- utils::choose.dir(default = default %||% getwd())
                                      if (is.na(sel)) sel <- ""
                            } else {                                               # Fallback: prompt for a path
                       prompt <- sprintf("Enter library directory path [%s]: ",
                                         default %||% getwd())
                          sel <- readline(prompt)
                                 if (!nzchar(sel)) sel <- default %||% getwd()
                           }
                           if (!nzchar(sel) || !dir.exists(sel)) 
                              stop("No valid library directory chosen.")
                    sel <- normalizePath(sel, winslash = "/")
                           assign("lib_dir", sel, envir = .GlobalEnv)             # persist for next run
                           message("Using library directory: ", sel)  
                           sel
                         }
# global best pattern shift
best_shift_metrics <- function(x, y_tt, y_I, grid, max_shift, shift_step,
                               enforce_nonneg = TRUE, eps = 1e-9) {
            same_len  <- length(y_tt) == length(grid)                             # x: experimental on 'grid'; y: library (tt = y_tt, I = y_I)
            same_ends <- isTRUE(all.equal(c(y_tt[1], y_tt[length(y_tt)]), 
                                c(grid[1], grid[length(grid)]), tolerance = 1e-6))
            same_step <- isTRUE(all.equal(detect_step(y_tt), 
                                          detect_step(grid), tolerance = 1e-9))
            fast_mode <- same_len && same_ends && same_step
                 step <- detect_step(grid)
               shifts <- seq(-max_shift, max_shift, by = shift_step)
            roll_by_k <- function(v, k) { 
                       n <- length(v)
                            if (k > 0) c(rep(NA_real_, k), v[1:(n - k)]) 
                            else if (k < 0) c(v[(1 - k):n], rep(NA_real_, -k)) 
                            else v
                         }
          scaled_rmse <- function(xs, ys) {                                       # compute scale with non-negative-residual
                    a_ls <- sum(xs * ys) / sum(ys * ys)
                       a <- max(0, a_ls)                                          # LS scale (clamped at 0)
                            if (enforce_nonneg) {
                               pos <- ys > eps & is.finite(xs)
                               if (any(pos)) {                                    # cap so x - a*y >= 0 for all points with y > eps
                                  a_cap <- min(xs[pos] / ys[pos], na.rm = TRUE)
                                  if (is.finite(a_cap)) a <- min(a, a_cap)
                               }
                            }
                      rm <- sqrt(mean((xs - a * ys)^2))
                            list(a = a, rm = rm)
                         }
                 best <- list(r= -Inf,shift= 0,cos= NA_real_, rmse = Inf, a = 1)
                         for (s in shifts) {
                 y_shift <- if (fast_mode && abs(s/step-round(s/step)) < 1e-6) {
                          k <- as.integer(round(s / step)); roll_by_k(y_I, k)
                            } else { approx(y_tt + s, y_I, xout = grid, rule = 2, 
                                            ties = mean)$y
                            }
                      ok <- is.finite(x) & is.finite(y_shift)
                            if (!any(ok)) next
                      xs <- pmax(x[ok], 0)                                        # guard against tiny negatives after baseline
                      ys <- pmax(y_shift[ok], 0)
                            if (sd(xs) == 0 || sd(ys) == 0) next
                      r  <- suppressWarnings(cor(xs, ys))
                     cos <- sum(xs * ys) / sqrt(sum(xs^2) * sum(ys^2))
                     sr  <- scaled_rmse(xs, ys)
                            if (is.finite(r) && r > best$r) {
                       best <- list(r = r, shift = s, cos = cos, rmse = sr$rm, a = sr$a)
                            }
                         }
                         best
                      }
# Area-preserving Gaussian “amorphization” steps
.gauss_kernel <- function(fwhm, step_deg, n_sigma = 3) {
                    if (!is.finite(fwhm) || fwhm <= 0) return(1)                  # degenerate without smoothing
           sigma <- fwhm / (2 * sqrt(2 * log(2)))
           half  <- max(1L, ceiling(n_sigma * sigma / step_deg))
           xk    <- (-half:half) * step_deg
           k     <- exp(-0.5 * (xk / sigma)^2)
                    k / sum(k)                                                    # unit area
                 }
amorphize_gauss <- function(I, step_deg, fwhm) {
                      if (!is.finite(fwhm) || fwhm <= 0) return(as.numeric(I))
                 k <- .gauss_kernel(fwhm, step_deg)
               pad <- length(k) %/% 2
                 n <- length(I)
              ypad <- c(rev(I[2:(pad+1)]), I, rev(I[(n-pad):(n-1)]))              # pad grid space to avoid edge losses and center
             yconv <- stats::filter(ypad, k, sides = 2, circular = FALSE)
                      as.numeric(yconv[(pad+1):(pad+n)])
                   }

# Peak extraction (library)
extract_peaks <- function(tt, I, det = DEFAULT_DET_LIB) {
             det <- validate_list(det, DEFAULT_DET_LIB)
                    with(det, {
                  I <- pmax(I, 0)
                  n <- length(I)
                       if (!n || !any(is.finite(I))) {
                          return(data.frame(pos = numeric(0), w = numeric(0), 
                                            is_imax = logical(0))) 
                       }
          Imax_val  <- max(I, na.rm = TRUE) 
                       if (!is.finite(Imax_val) || Imax_val <= 0) Imax_val <- 1
          imax_idx  <- which.max(I)
                d2  <- diff(sign(diff(I)))                                         # enforce local maxima
                idx <- which(c(FALSE, d2 == -2, FALSE))
                       if (!(imax_idx %in% idx)) idx <- sort(unique(c(idx, imax_idx)))
               rel  <- I[idx] / Imax_val                                           # detect relative heights
               keep <- which(rel >= rel_min) 
                       if (!length(keep)) keep <- which.max(rel)                   # always keep tallest peak
               idx  <- idx[keep]
                rel <- rel[keep]
                pos <- tt[idx]
                  w <- rel                          
                       if (is.finite(w_gamma) && w_gamma != 1) {                   # dynamic-range compression
                     w <- w^W_gamma
                     w <- w / max(w, na.rm = TRUE)
                       }
              ord   <- order(-w)                                                   # enforce minimum spacing
              sel   <- logical(length(idx))
              taken <- rep(FALSE, length(idx))
                       for (k in ord) {
                          if (taken[k]) next
                sel[k] <- TRUE
                taken  <- taken | (abs(pos - pos[k]) < min_spacing)
                       }
            idx_sel <- idx[sel]
            pos_sel <- pos[sel]
            w_sel   <- w[sel]
        is_imax_sel <- idx_sel == imax_idx                                         # mark peak at global Imax
                       data.frame(pos = pos_sel, w = w_sel, is_imax = is_imax_sel, row.names = NULL)
                    })
                }
# library indexing
build_library_index <- function(lib_dir, rel_min = 0.20, min_spacing = 0.06) {
  files <- list.files(lib_dir, pattern = "\\.xy$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No .xy files in: ", lib_dir)
  idx <- lapply(files, function(f) {
    dt <- read_xy(f)
    pk <- extract_peaks(dt$tt, dt$I)
    list(file = basename(f), pos = pk$pos, w = pk$w, is_imax = pk$is_imax)
  })
  structure(idx, class = "lib_peaks_index", dir = lib_dir, built = Sys.time())
}
save_lib_index <- function(index, lib_dir) {
  out_dir <- file.path(lib_dir, "Results"); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  saveRDS(index, file.path(out_dir, "lib_peaks_index.rds"))
}
load_lib_index <- function(lib_dir) {
  path <- file.path(lib_dir, "Results", "lib_peaks_index.rds")
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}
# light smoothing for detection only (does NOT change intensities used for weights)
.smooth_ma <- function(y, k = 5) {
  k <- max(3, as.integer(k)); if (k %% 2 == 0) k <- k + 1
  pad <- (k - 1)/2
  ypad <- c(rep(y[1], pad), y, rep(y[length(y)], pad))
  as.numeric(stats::filter(ypad, rep(1/k, k), sides = 2))[ (pad+1):(length(ypad)-pad) ]
}

# local maxima on a smoothed trace
.local_max <- function(z) {
  n <- length(z)
  which( (c(-Inf, z[-n]) < z) & (z >= c(z[-1], -Inf)) )                         # preserve right edge of plateaus
}

# identify experimental peaks (assume data is noisy)
extract_peaks_exp <- function(tt, I, det = DEFAULT_DET_EXP){
  det = DEFAULT_DET_EXP
  with(det, {
  I <- pmax(I, 0)
  n <- length(I)
  if (!n || !any(is.finite(I))) {
    return(data.frame(pos = numeric(0), w = numeric(0), is_imax = logical(0)))
  }
  Imax_val <- max(I, na.rm = TRUE); if (!is.finite(Imax_val) || Imax_val <= 0) Imax_val <- 1
  imax_idx <- which.max(I)
  Js <- .smooth_ma(I, k = smooth_k)                                             # find peak candidates on smoothed series
  cand <- .local_max(Js)
  if (!length(cand)) cand <- imax_idx
  rel <- I[cand] / Imax_val
  step <- max(1L, round( (prom_window) / max(1e-9, median(diff(tt), na.rm = TRUE)) ))
  get_local_min <- function(i) min( I[max(1, i - step) : min(n, i + step)] , na.rm = TRUE)
  base_loc <- vapply(cand, get_local_min, numeric(1))
  prom <- I[cand] - base_loc
  prom_rel <- prom / Imax_val
  keep <- which( (rel >= rel_min) | (prom_rel >= prominence_rel) )
  if (!length(keep)) keep <- which.max(rel)                                     # accept peaks if high OR prominent
  idx  <- cand[keep]; rel <- rel[keep]                                          # always keep strongest peak
  ord   <- order(-rel)                                                          # enforce minimum spacing
  sel   <- logical(length(idx))
  taken <- rep(FALSE, length(idx))
  pos0  <- tt[idx]
  for (k in ord) {
    if (taken[k]) next
    sel[k] <- TRUE
    taken  <- taken | (abs(pos0 - pos0[k]) < min_spacing)
  }
  idx <- idx[sel]; rel <- rel[sel]                                              # relax thresholds if necessary for results
  if (length(cand) > length(idx) && length(idx) < min_keep) {
    extra <- setdiff(cand, idx)
    rel_extra <- I[extra] / Imax_val
    ord_extra <- order(-rel_extra)
    for (j in ord_extra) {
      i <- extra[j]
      if (any(abs(tt[idx] - tt[i]) < min_spacing)) next
      idx <- c(idx, i); rel <- c(rel, I[i]/Imax_val)
      if (length(idx) >= min_keep) break
    }
  }
  w <- rel                                                                      # weight compression and normalization
  if (is.finite(w_gamma) && w_gamma != 1) {
    w <- w^w_gamma
    w <- w / max(w, na.rm = TRUE)
  }
  if (!(imax_idx %in% idx)) {                                                   # ensure global Imax is present and marked
    idx <- c(idx, imax_idx); w <- c(w, 1.0)
  }
  ord_final <- order(tt[idx])
  pos <- tt[idx][ord_final]
  w   <- w[ord_final]
  is_imax <- (idx[ord_final] == imax_idx)
  data.frame(pos = pos, w = w, is_imax = is_imax, row.names = NULL)
  })
}

# Peak matching
match_score_shift <- function(exp_peaks, lib_peaks, shift,
                              tol = NULL, heavy_thr = NULL, heavy_penalty = NULL,
                              weight_gamma = 1.0, return_pairs = FALSE,
                              cfg = NULL) {
  eff <- if (is.null(cfg)) DEFAULT_CFG else validate_list(cfg, DEFAULT_CFG)
  if (!is.null(tol))           eff$tol           <- suppressWarnings(as.numeric(tol))
  if (!is.null(heavy_thr))     eff$heavy_thr     <- suppressWarnings(as.numeric(heavy_thr))
  if (!is.null(heavy_penalty)) eff$heavy_penalty <- suppressWarnings(as.numeric(heavy_penalty))
  if (!is.null(weight_gamma))  eff$weight_gamma  <- suppressWarnings(as.numeric(weight_gamma))
  if (!is.finite(eff$weight_gamma)) eff$weight_gamma <- 1.0
  v <- suppressWarnings(as.numeric(lib_peaks$w %||% numeric(0)))
  v[!is.finite(v)] <- 0
  if (abs(eff$weight_gamma - 1) > 1e-9) v <- v ^ eff$weight_gamma
  if (!nrow(exp_peaks) || !length(v)) {
    heavy_mask <- (lib_peaks$w %||% numeric(0)) >= eff$heavy_thr
    total_lib_w <- sum(v); heavy_total_w <- sum(v[heavy_mask])
    out <- list(score=0, matches=0L,
                matched_lib_w=0, total_lib_w=total_lib_w,
                heavy_missed_w=0, heavy_total_w=heavy_total_w,
                pairs_idx=list(e=integer(0), l=integer(0)))
    return(out)
  }
  ord_e <- order(exp_peaks$pos)                                                 # map back to original row indices
  ep    <- exp_peaks[ord_e, , drop = FALSE]
  ord_l <- order(lib_peaks$pos)
  lp <- data.frame(
    pos = lib_peaks$pos[ord_l] + shift,
    w   = v[ord_l],
    is_heavy = (lib_peaks$w[ord_l] >= eff$heavy_thr))
  i <- j <- 1L
  matched_l <- logical(nrow(lp))
  e_idx_ord <- integer(0); l_idx_ord <- integer(0)
  while (i <= nrow(ep) && j <= nrow(lp)) {
    d <- lp$pos[j] - ep$pos[i]
    if (abs(d) <= eff$tol) {
      if (!matched_l[j]) { matched_l[j] <- TRUE; e_idx_ord <- c(e_idx_ord, i); l_idx_ord <- c(l_idx_ord, j) }
      i <- i + 1L; j <- j + 1L
    } else if (d < -eff$tol) j <- j + 1L else i <- i + 1L
  }
    matched_lib_w <- sum(lp$w[matched_l])
    total_lib_w   <- sum(lp$w)
    heavy_total_w  <- sum(lp$w[lp$is_heavy])
    heavy_missed_w <- if (heavy_total_w > 0) sum(lp$w[lp$is_heavy & !matched_l]) else 0
    coverage <- if (total_lib_w > 0) matched_lib_w / total_lib_w else 0
    penalty  <- if (heavy_total_w > 0) eff$heavy_penalty * (heavy_missed_w / heavy_total_w) else 0
    score    <- max(0, coverage - penalty)
    list(score=score, matches=as.integer(sum(matched_l)),
                matched_lib_w=matched_lib_w, total_lib_w=total_lib_w,
                heavy_missed_w=heavy_missed_w, heavy_total_w=heavy_total_w,
                pairs_idx=if (isTRUE(return_pairs))
                list(e=ord_e[e_idx_ord], l=ord_l[l_idx_ord]) else NULL)
}

# Fit peak matches
best_shift_for_library <- function(exp_peaks, lib_peaks, scan = DEFAULT_SCAN, 
                                   cfg  = DEFAULT_CFG) {
                     scan <- validate_list(scan, DEFAULT_SCAN)
                     cfg  <- validate_list(cfg,  DEFAULT_CFG)
                     best <- list(score = -Inf, shift = NA_real_, matches = 0L,
                                  matched_lib_w = 0, total_lib_w = 0,
                                  heavy_missed_w = 0, heavy_total_w = 0,
                                  pairs_idx= list(e= integer(0), l= integer(0)))
                   shifts <- seq(-scan$max_shift,scan$max_shift,by=scan$shift_step) # coarse scan (use scan$*)
                             for (s in shifts) {
                          sc <- match_score_shift(exp_peaks, lib_peaks, shift = s, cfg = cfg, return_pairs = TRUE)
                                if (is.finite(sc$score) && sc$score > best$score) best <- c(sc, list(shift = s))
                                }
                             if (isTRUE(scan$refine) && is.finite(best$shift)){     # optional refinement (use scan$*)
                        fine <- seq(best$shift - scan$refine_win, best$shift +
                                    scan$refine_win, by = scan$refine_step)
                                for (s in fine) {
                                   sc <- match_score_shift(exp_peaks, lib_peaks, 
                                           shift = s, cfg = cfg, return_pairs = TRUE)
                                   if (is.finite(sc$score) && sc$score > best$score) 
                              best <- c(sc, list(shift = s))
                                }
                             }
                             best
                          }

#
## --- peak matching script ---
#
compare_peaks_to_library <- function(exp_xy_path, lib_dir,
                                     rel_min     = DEFAULT_DET_EXP$rel_min,
                                     min_spacing = DEFAULT_DET_EXP$min_spacing,
                                     max_shift   = DEFAULT_SCAN$max_shift,
                                     shift_step  = DEFAULT_SCAN$shift_step,
                                     tol         = DEFAULT_CFG$tol,
                                     overlay_k   = DEFAULT_OVERLAY_K) {
  scan    <- list(max_shift = max_shift, shift_step = shift_step)
  cfg     <- list(tol = tol)
  message("Library directory: ", normalizePath(lib_dir, winslash = "/", mustWork = FALSE))
  dt  <- read_xy(exp_xy_path)                                                   # determine experimental peaks
  exp_peaks <- extract_peaks_exp(dt$tt, dt$I)
  idx <- load_lib_index(lib_dir)                                                # index library
  if (is.null(idx)) {
    message("No library peak index found — building now (once).")
    idx <- build_library_index(lib_dir, rel_min = rel_min, min_spacing = min_spacing)
    save_lib_index(idx, lib_dir)
    message("Saved index: ", file.path(lib_dir, "Results", "lib_peaks_index.rds"))
  } else {
    message("Loaded library peak index from Results/lib_peaks_index.rds")
  }
  res <- lapply(idx, function(entry) {                                          # evaluate library
    if (is.null(entry$is_imax)) {
      ii <- if (length(entry$w)) which.max(entry$w) else integer(0)
      entry$is_imax <- rep(FALSE, length(entry$w))
      if (length(ii)) entry$is_imax[ii] <- TRUE
    }
     best <- best_shift_for_library(exp_peaks, entry, scan=scan, cfg=cfg)       
     exp_n   <- nrow(exp_peaks)
     lib_n   <- length(entry$pos)
     pairs_n <- best$matches
     exp_imax_row <- which(exp_peaks$is_imax)
     lib_imax_row <- which(entry$is_imax)
     exp_imax_matched <- length(exp_imax_row) > 0 && any(best$pairs_idx$e %in% exp_imax_row)
     lib_imax_matched <- length(lib_imax_row) > 0 && any(best$pairs_idx$l %in% lib_imax_row)
     if (best$matches > 0) {                                                    # Calculate peak offset residuals after shift
       e_idx <- best$pairs_idx$e
       l_idx <- best$pairs_idx$l
       exp_pos <- exp_peaks$pos[e_idx]
       lib_pos <- entry$pos[l_idx] + best$shift
       deltas  <- lib_pos - exp_pos
       delta_mean_abs   <- mean(abs(deltas))
       delta_median_abs <- stats::median(abs(deltas))
       delta_max_abs    <- max(abs(deltas))
       frac_small_delta <- mean(abs(deltas) <= (tol/2))
     } else {
       delta_mean_abs <- delta_median_abs <- delta_max_abs <- frac_small_delta <- NA_real_
     }
     data.frame(                                                                # Build result metadata and produce result files
       file             = entry$file,
       mineral          = sub("\\.xy$", "", entry$file, ignore.case = TRUE),
       score            = as.numeric(best$score),                               # library-only, penalized by heavy misses
       best_shift_deg   = as.numeric(best$shift),                               # single global shift for this file
       matches          = as.integer(best$matches),                             # count of matched library peaks
       matched_lib_w    = as.numeric(best$matched_lib_w),
       total_lib_w      = as.numeric(best$total_lib_w),
       heavy_missed_w   = as.numeric(best$heavy_missed_w),
       heavy_total_w    = as.numeric(best$heavy_total_w),
       exp_peaks_n      = as.integer(nrow(exp_peaks)),
       lib_peaks_n      = as.integer(length(entry$pos)),
       exp_imax_matched = exp_imax_matched,
       lib_imax_matched = lib_imax_matched,
       stringsAsFactors = FALSE
     )
  })
  tab <- do.call(rbind, res)
  tab <- tab[order(-tab$score, -tab$matched_lib_w, -tab$matches, abs(tab$best_shift_deg)), , drop = FALSE]
  rownames(tab) <- NULL
  base14 <- gsub("[^A-Za-z0-9_.-]", "_", substr(sub("\\.xy$", "", basename(exp_xy_path), ignore.case = TRUE), 1, 14))
  out_dir <- file.path(dirname(exp_xy_path), "Results"); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_csv <- file.path(out_dir, paste0(base14, ".csv"))
  utils::write.csv(tab, out_csv, row.names = FALSE)
  message("Saved match report: ", out_csv)
  topk <- utils::head(tab, overlay_k)                                           # plot top K experimental vs shifted lib peaks
  for (i in seq_len(nrow(topk))) {
    entry <- idx[[ which(vapply(idx, function(e) e$file, "") == topk$file[i]) ]]
    shifted <- data.frame(pos = entry$pos + topk$best_shift_deg[i], w = entry$w)
    png_path <- file.path(out_dir, sprintf("%s_vs_%s_top%02d.png", base14, topk$mineral[i], i))
    grDevices::png(png_path, width = 1400, height = 900, res = 150)
    stick_plot(exp_peaks, shifted,
               sprintf("#%d  %s  (Score=%.3f, shift=%.1f°, matches=%d)",
                       i, topk$mineral[i], topk$score[i], topk$best_shift_deg[i], topk$matches[i]))
    grDevices::dev.off()
  }
  invisible(tab)
}

#
## --- .xy pattern matching script ---
#
compare_xy_to_library <- function(exp_xy_path, lib_dir,
                                  lib_step        = DEFAULT_LIB_STEP,
                                  max_shift       = DEFAULT_MAX_SHIFT,
                                  shift_step      = DEFAULT_SHIFT_STEP,
                                  top_n           = DEFAULT_TOP_N,
                                  overlay_k       = DEFAULT_OVERLAY_K,
                                  amorphize       = DEFAULT_AMORPH,                   
                                  amorph_fwhm_deg = DEFAULT_AMORPH_FWHM){
  message("Library directory: ", normalizePath(lib_dir, winslash = "/", mustWork = FALSE))
  if (isTRUE(amorphize)) {
    amorph_fwhm_deg <- as.numeric(amorph_fwhm_deg)
    if (!is.finite(amorph_fwhm_deg) || amorph_fwhm_deg < 0) amorph_fwhm_deg <- 0
    amorph_fwhm_deg <- min(amorph_fwhm_deg, 0.25)
  }
  if (isTRUE(amorphize) && amorph_fwhm_deg > 0) {
    message(sprintf("Amorphization: ON  (Gaussian FWHM = %.3f° 2θ) for all library patterns.", amorph_fwhm_deg))
  } else {
    message("Amorphization: OFF")
  }
  # load experimental xy file and library files
  exp_dt  <- read_xy(exp_xy_path)                                               # expects columns: tt, I
  grid    <- exp_dt$tt                                                          # canonical grid = experimental tt
  x       <- exp_dt$I                                                           # experimental intensity
  lib_files <- list.files(lib_dir, pattern = "\\.xy$", full.names = TRUE, ignore.case = TRUE)
  if (!length(lib_files)) stop("No .xy files found in: ", lib_dir)
  # prepare results and output file
  res <- vector("list", length(lib_files))
  base_dir  <- dirname(exp_xy_path)
  out_dir   <- file.path(base_dir, "Results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base_name <- sub("\\.xy$", "", basename(exp_xy_path), ignore.case = TRUE)
  base14    <- gsub("[^A-Za-z0-9_.-]", "_", substr(base_name, 1, 14))
  out_csv <- file.path(out_dir, paste0(base14, ".csv"))
  # iterate library evaluation and build result table
  for (i in seq_along(lib_files)) {
    f <- lib_files[i]
    li <- read_xy(f)
    yI <- li$I                                                                  # amorphization step
    if (isTRUE(amorphize) && amorph_fwhm_deg > 0) {
      yI <- amorphize_gauss(yI, detect_step(li$tt), amorph_fwhm_deg)
    }
    best <- best_shift_metrics(x, li$tt, yI, grid, max_shift, shift_step, enforce_nonneg = TRUE, eps = 1e-8)
    res[[i]] <- data.frame(
      file = basename(f),
      mineral = sub("\\.xy$", "", basename(f), ignore.case = TRUE),
      r = best$r,
      cosine = best$cos,
      rmse = best$rmse,
      best_shift_deg = best$shift,
      scale_a = best$a,
      amorphize = isTRUE(amorphize) && amorph_fwhm_deg > 0,
      amorph_fwhm_deg = if (isTRUE(amorphize) && amorph_fwhm_deg > 0) amorph_fwhm_deg else 0,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, res)
  tab <- tab[order(-tab$r, tab$rmse), ]
  rownames(tab) <- NULL
  if (is.null(out_csv)) {                                                       # write CSV
    base <- sub("\\.xy$", "", exp_xy_path, ignore.case = TRUE)
    out_csv <- paste0(base, "_match_report.csv")
  }
  safe_write_csv(out_csv, tab)
  if (is.null(out_dir)) out_dir <- dirname(exp_xy_path)                         # make plots and show results for top K
  overlay_k <- as.integer(overlay_k)
  topk <- utils::head(tab, overlay_k)
  for (i in seq_len(nrow(topk))) {
    f  <- file.path(lib_dir, topk$file[i])
    li <- read_xy(f)
    yI <- li$I                                                                  # reapply amorphization if enabled
    if (isTRUE(amorphize) && amorph_fwhm_deg > 0) {
      yI <- amorphize_gauss(yI, detect_step(li$tt), amorph_fwhm_deg)
    }
    y_shift <- approx(li$tt + topk$best_shift_deg[i], yI, xout = grid, rule = 2, ties = mean)$y
    y_plot  <- topk$scale_a[i] * y_shift
    title_tag <- sprintf(
      "#%d  %s  (Score=%.3f, shift=%.1f°, matches=%d, libW=%.2f/%.2f, heavyMiss=%.2f/%.2f)",
      i, topk$mineral[i], topk$score[i], topk$best_shift_deg[i], topk$matches[i],
      topk$matched_lib_w[i], topk$total_lib_w[i],
      topk$heavy_missed_w[i], topk$heavy_total_w[i],
      if (isTRUE(amorphize) && amorph_fwhm_deg > 0) sprintf(", amorph=%.3f°", amorph_fwhm_deg) else ""
    )
    png_path <- file.path(out_dir, sprintf("%s_vs_%s_top%02d.png", base14, topk$mineral[i], i))
    safe_write_png(png_path, function() {
      overlay_plot(grid, x, y_plot, title = title_tag, diff = TRUE, file = NULL)
    })
  }
  message("Saved match report: ", out_csv)
  message("Saved overlay PNGs for top ", nrow(topk), " in: ", out_dir)
  print(utils::head(tab, top_n))
  message("Saved match report: ", out_csv)
  message("Saved overlay PNGs for top ", nrow(topk), " in: ", out_dir)
  invisible(tab)
}

#
## --- auto-run (top menu) ---
#
if (interactive()) {
  run_matching_menu()
}
