#
####===== Preprocessing techniques for raw diffratogram spectra (.xy) =====####
#
#
###--- Default parameter settings ---###
#
#
# preprocessing parameters:          current settings:              suggested values:
      DEFAULT_REL_I       <- c(0.21, 1.00, 0.46, 0.13, 0.08)        # c(1.00, 0.46, 0.21, 0.13, 0.08) for 5 line aluminum
#
###--- Dependency check ---###
#
#
required_pkgs <- c("data.table","baseline","prospectr","limSolve","signal","waveslim")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing)) {
   stop("Required packages missing: ", paste(missing, collapse=", "), "\nInstall with: install.packages(c(", 
        paste(sprintf('\"%s\"', missing), collapse=", "), "))", call. = FALSE)
}
library(data.table)                                                               # data.table attached, others called directly with ::
#
#
###--- Helper functions ---###
#
#
.pick_file <- function() {
                if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
                   p <- rstudioapi::selectFile(caption = "Choose experimental .xy file",
                                               path = getwd(), filter = "XY files (*.xy)|*.xy")
                   if (!is.null(p) && nzchar(p)) return(p)
                }
                message("Using base file chooser (RStudio API not available).")
                suppressWarnings(file.choose(new = FALSE))
              }
read_xy <- function(path) {
             dt <- data.table::fread(path, header = FALSE, col.names = c("tt","I"))
             dt[order(tt)]
           }
plot_before_after <- function(tt, y_before, y_after, cex_main_before = 0.85, 
                              cex_main_after = 0.85, main_before = "Raw", 
                              main_after="Baseline subtracted", wrap = 50) {
                  mb <- if (wrap) paste(strwrap(main_before, width = wrap), 
                                        collapse = "\n") 
                           else main_before
                  ma <- if (wrap) paste(strwrap(main_after,  width = wrap), 
                                        collapse = "\n") 
                           else main_after
                  op <- par(no.readonly = TRUE); on.exit(par(op))
                        par(mfrow = c(1,2), mar = c(4,4,4.5,1))
                        on.exit(par(op), add = TRUE)
                        plot(tt, y_before, type="l", xlab=expression(2*theta~"(deg)"), 
                             ylab="Intensity", main=mb, cex.main=cex_main_before)
                        plot(tt, y_after, type="l", xlab=expression(2*theta~"(deg)"), 
                             ylab="Intensity", main=ma, cex.main=cex_main_after)
                     }
write_xy_file <- function(path, x, y, digits_x = 5, digits_y = 6) {               # 2 column text, no headers
                   fmt <- paste0("%.", digits_x, "f %.", digits_y, "f")
                   con <- file(path, "wt"); on.exit(close(con), add = TRUE)
                   for (i in seq_along(x)) writeLines(sprintf(fmt, x[i], y[i]), con)
                       invisible(path)
                 }

# Output path = input path
make_out_path <- function(infile, suffix = "_prepped.xy") {
                   if (grepl("\\.xy$", infile, ignore.case = TRUE)) {
                   base <- sub("\\.xy$", "", infile, ignore.case = TRUE)
                   out  <- paste0(base, suffix)
                   } else {
                       base <- infile
                       out  <- paste0(base, suffix)
                   }
                   if (!file.exists(out)) return(out)
                   root <- sub("\\.xy$", "", out, ignore.case = TRUE)
                   for (i in 1:999) {
                       cand <- paste0(root, "-", sprintf("%02d", i), ".xy")
                       if (!file.exists(cand)) return(cand)
                   }
                   stop("Could not create a unique output filename.")
                 }

#
###--- Recipe tracking ---###
#

apply_recipe_to_file <- function(file_path, recipe, save_baseline = FALSE, suffix = "_prepped.xy") {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  dt <- read_xy(file_path)
  rawI <- dt$I
  base <- recipe$baseline                                                         # baseline from the recipe
  corrI <- do.call(baseline_subtract_vec, c(list(I = rawI, method = base$method), base$args))
  y <- corrI
  for (op in recipe$ops) {                                                        # replay steps
    nm <- op$name
    if (nm == "sg") {
      y <- sg_smooth(y, n = op$n, p = op$p)
    } else if (nm == "wavelet") {
      y <- wavelet_denoise(y, wf = op$wf, L = op$L, mode = op$mode,
                           alpha = op$alpha, level_min = op$level_min)
    } else if (nm == "norm_max1") {
      y <- normalize_max1(y)
    } else if (nm == "norm_area1") {
      step_est <- median(diff(dt$tt))
      y <- normalize_area1(y, step = step_est)
    } else if (nm == "scale_Imax") {
      y <- normalize_to_Imax(y, target_max = op$target_max)
    }
  }
  if (grepl("\\.xy$", file_path, ignore.case = TRUE)) {
    out_prepped  <- sub("\\.xy$", paste0(suffix), file_path, ignore.case = TRUE)
    out_baseline <- sub("\\.xy$", "_baseline.xy", file_path, ignore.case = TRUE)
  } else {
    out_prepped  <- paste0(file_path, suffix)
    out_baseline <- paste0(file_path, "_baseline.xy")
  }
  write_xy_file(out_prepped,  dt$tt, y)
  if (isTRUE(save_baseline)) write_xy_file(out_baseline, dt$tt, corrI)
  invisible(list(prepped = out_prepped, baseline = if (save_baseline) out_baseline else NULL))
}

# Apply the recorded recipe to a dataframe and return series
apply_recipe_to_data <- function(dt, recipe) {
  rawI <- dt$I
  base <- recipe$baseline
  corrI <- do.call(baseline_subtract_vec, c(list(I = rawI, method = base$method), base$args))
  y <- corrI
  for (op in recipe$ops) {
    nm <- op$name
    if (nm == "sg") {
      y <- sg_smooth(y, n = op$n, p = op$p)
    } else if (nm == "wavelet") {
      y <- wavelet_denoise(y, wf = op$wf, L = op$L, mode = op$mode,
                           alpha = op$alpha, level_min = op$level_min)
    } else if (nm == "norm_max1") {
      y <- normalize_max1(y)
    } else if (nm == "norm_area1") {
      step_est <- median(diff(dt$tt))
      y <- normalize_area1(y, step = step_est)
    } else if (nm == "scale_Imax") {
      y <- normalize_to_Imax(y, target_max = op$target_max)
    }
  }
  list(corrI = corrI, y = y)
}

#
###--- transform functions ---###
#

normalize_max1  <- function(y) { 
                     m <- max(y, na.rm=TRUE) 
                     if (!is.finite(m) || m==0) y else y/m 
                   }
normalize_area1 <- function(y) {
                     s <- sum(y, na.rm=TRUE)
                     if (!is.finite(s) || s==0) y else y/s 
                   }
moving_average <- function(y, k = 3) {
                    if (k <= 1) return(y)
                    filt <- rep(1/k, k)
                    z <- as.numeric(stats::filter(y, filt, sides = 2))
                    z[is.na(z)] <- 0
                    z
                  }

# log-SNIP with optional auto offset
snip_baseline <- function(y, iterations = 60, log = FALSE, offset = NULL) {
                   y <- as.numeric(y)
                   n <- length(y)
                   if (!is.finite(iterations) || iterations < 1) return(rep(0, n))
                   if (log) {                                                     # choose a stable positive offset if not supplied
                      if (is.null(offset)) {
                         yy <- y[is.finite(y) & y > 0]
                         offset <- if (length(yy)) max(1e-6, 0.01 * median(yy)) 
                         else 1e-6
                      }
                   z <- log(pmax(y, 0) + offset)
                   b <- z
                   for (k in 1:iterations) {
                       zl <- c(rep(b[1], k), b[1:(n - k)])
                       zr <- c(b[(k + 1):n], rep(b[n], k))
                       avg <- 0.5 * (zl + zr)
                       b <- pmin(b, avg)
                   }
                   bl <- pmax(exp(b) - offset, 0)                                 # back-transform baseline
                   } else {
                          b <- y
                          for (k in 1:iterations) {
                              yl <- c(rep(b[1], k), b[1:(n - k)])
                              yr <- c(b[(k + 1):n], rep(b[n], k))
                              avg <- 0.5 * (yl + yr)
                              b <- pmin(b, avg)
                          }
                          bl <- b
                   }
                   bl
                }

# scale working spectra
normalize_to_Imax <- function(y, target_max = 10000) {
                       m <- max(y, na.rm = TRUE)
                       if (!is.finite(m) || m == 0) return(y)
                       y * (target_max / m)
                     }

build_grid <- function(tt_min, tt_max, step_deg) {
                if (!is.finite(step_deg) || step_deg <= 0) 
                   stop("Library step must be > 0")
                seq(from = tt_min, to = tt_max, by = step_deg)
              }

# unified subtractor (supports SNIP with log/offset)
baseline_subtract_vec <- function(I, method = "rollingBall", ...) {
                           y <- as.numeric(I)
                           if (toupper(method) == "SNIP") {
                              dots <- list(...)
                              iters  <- if (!is.null(dots$iterations)) 
                                           as.integer(dots$iterations) 
                                        else 60L
                              uselog <- isTRUE(dots$log)
                              off    <- if (!is.null(dots$offset)) dots$offset 
                                        else NULL
                              bl <- snip_baseline(y, iterations = iters, 
                                                  log = uselog, offset = off)
                              corr <- y - bl
                           } else {
                               m <- matrix(y, nrow = 1)
                               blfit <- baseline::baseline(m, method = method, ...)
                               corr <- as.numeric(baseline::getCorrected(blfit))
                           }
                           corr[corr < 0] <- 0
                           corr
                        }

# Savitzky–Golay smoothing (poly order p, window n must be odd)
sg_smooth <- function(y, n = 5, p = 2) {
               if (n %% 2 == 0) n <- n + 1                                        # make sure n is odd
               if (n <= p) n <- p + 3 + (p %% 2)                                  # ensure window > poly order
                  signal::sgolayfilt(y, p = p, n = n)
             }

# Wavelet denoising (configurable):
     # wf:        waveslim wavelet name ("la8" = Daubechies db4; "la10" ~ db5; "la12" ~ db6)
     # L:         number of decomposition levels (higher = more aggressive, typically 3–6)
     # mode:      "soft" (default) or "hard" threshold
     # alpha:     threshold multiplier (1.0 default, >1.0 = more aggressive)
     # level_min: set thresholds from at this detail level (1 = finest details)
wavelet_denoise <- function(y, wf = "la8", L = 4, mode = c("soft","hard"), alpha = 1.0, level_min = 1) {
                     mode <- match.arg(mode)
                     n <- length(y)
                     pow2 <- 2^ceiling(log2(n))
                     pad_len <- pow2 - n
                     ypad <- if (pad_len > 0) c(y, rep(tail(y, 1), pad_len))      # padded margins for y
                             else y
                     wf <- resolve_wavelet_name(wf)                               # ensure wf call is valid (check aliases)
                     L <- clamp_levels(L, length(ypad))                           # clamp L to allowable range 
                     d <- waveslim::dwt(ypad, wf = wf, n.levels = L,              # discrete wavelet transform default settings
                                        boundary = "periodic")
                     sigma_hat <- function(w) stats::mad(w, center = 0) / 0.6745
                     shrink <- function(w, thr, mode) 
                                 if (mode == "soft") sign(w)*pmax(abs(w)-thr,0) 
                                 else w*(abs(w)>thr)
                                 for (lvl in seq_len(L)) {
                                     if (lvl >= level_min) {
                                        w   <- d[[paste0("d", lvl)]]
                                        thr <- alpha * sqrt(2 * log(length(w))) * sigma_hat(w)
                                        d[[paste0("d", lvl)]] <- shrink(w, thr, mode)
                                     }
                                 }
                                 yrec <- waveslim::idwt(d)
                                 yrec[seq_len(n)]
                     }
# Preview 3 wavelet variants (Daubechies db4 = waveslim "la8")
wavelet_preview3 <- function(x, y, variants = list(list(L = 4, alpha = 1.0),
                                                   list(L = 5, alpha = 1.5),
                                                   list(L = 6, alpha = 2.0)),
                             wf = "la8", mode = "soft", level_min = 1) {
                      ys <- lapply(variants, function(v) wavelet_denoise(y,       # compute three candidates
                                   wf = wf, L = v$L, mode = mode, alpha = v$alpha, 
                                   level_min = level_min))                        # labels (coerce types explicitly to satisfy vapply)
                      Ls <- vapply(variants, function(v) 
                                   as.integer(v[["L"]]), integer(1))
                      alphas <- vapply(variants, function(v) 
                                       as.numeric(v[["alpha"]]), numeric(1))
                      labs <- sprintf("db4 soft  L=%d  \u03B1=%.1f", Ls, alphas)
                      op <- par(no.readonly = TRUE)
                      on.exit(par(op), add = TRUE)                                # 2x2 plot: raw + 3 variants
                      par(mfrow = c(2,2), mar = c(4,4,2.5,1))
                      plot(x, y, type = "l", main = "Raw (before)",
                           xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
                      plot(x, ys[[1]], type = "l", main = labs[1],
                           xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
                      plot(x, ys[[2]], type = "l", main = labs[2],
                           xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
                      plot(x, ys[[3]], type = "l", main = labs[3],
                           xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
                      cat("\nPreviewed variants:\n  [1]", labs[1], "\n  [2]", 
                          labs[2], "\n  [3]", labs[3], "\n  [0] Cancel\n")        # user choose which to apply
                      pick <- suppressWarnings(as.integer(readline(
                                               "Apply which variant? (0/1/2/3): ")))
                      if (!is.finite(pick) || pick == 0) return(NULL)
                      if (pick %in% 1:3) {
                         return(list(y = ys[[pick]], label = sprintf(
                                "Wavelet(db4,soft,L=%d,alpha=%.1f,level_min=%d)",
                                Ls[pick], alphas[pick], level_min)))
                      }
                   NULL
                   }
# Map aliases to waveslim codes and validate
resolve_wavelet_name <- function(wf_in) {
                          aliases <- c("db2"="d4","db4"="la8","db5"="la10","db6"="la12",
                                       "db7"="la14","db8"="la16","haar"="haar")
                          valid <- c("haar","d4","la8","la10","la12","la14","la16",
                                     "coif1","coif2","coif3","bl14","bl18","bl20",
                                     "mb4","mb8","fk4","fk6","fk8","fk14","fk18","fk22")
                          w <- tolower(trimws(wf_in))
                          w <- gsub("\\s+", "", w)
                          if (w %in% names(aliases)) w <- aliases[[w]]
                          if (!w %in% valid) {
                             message("Unrecognized wavelet '", wf_in, "'. Using 'la8' (Daubechies db4).")
                             w <- "la8"
                          }
                          w
                        }
# Clamp decomposition level to what the padded length allows
clamp_levels <- function(L, n) {
                  pow2 <- 2^ceiling(log2(n))                                      # padded length used in dwt 
                  maxL <- floor(log2(pow2))                                       # maximum allowable levels
                  Lc <- max(1L, min(as.integer(L), maxL)) 
                  if (Lc != L) message("Adjusted levels L from ", L, " to ", Lc, 
                                       " (signal length constraint).")
                  Lc
                }

# non-negative least squares scalar
.nnls_scalar <- function(x, y) {
                  den <- sum(x^2, na.rm = TRUE)
                  if (den <= 0) return(0)                                         # for alpha ≥ 0, y = alpha * x
                  max(0, sum(x*y, na.rm = TRUE) / den)
                }
# Convert to and from sigma (standard deviation) with FWHM measurements (in deg 2θ)
.fwhm_sigma <- function(fwhm) fwhm / (2*sqrt(2*log(2)))
.sigma_fwhm <- function(sigma) 2*sqrt(2*log(2)) * sigma

# Build a matrix of Gaussian or Lorentzian curves at specified positions based on FWHM and intensity measurements
.draw_curves <- function(tt, pos, w, fwhm, shape = c("gauss","lorentz")) {
                  stopifnot(is.numeric(tt), is.numeric(pos), is.numeric(w),
                  length(pos) == length(w), is.numeric(fwhm), fwhm > 0)
                  shape <- match.arg(shape)
                  if (shape == "gauss") {
                     sig <- .fwhm_sigma(fwhm)
                     vapply(seq_along(pos), function(i) {w[i] * exp(-0.5 * 
                            ((tt - pos[i]) / sig)^2)}, numeric(length(tt)))
                  } else {
                      g <- fwhm / 2
                      vapply(seq_along(pos), function(i) {w[i] * (g / ((tt - 
                             pos[i])^2 + g^2))}, numeric(length(tt)))
                  }
                }

#
#
###--- sample holder corrections (parameterized for aluminum) ---###
#
#

#  theoretical diffraction pattern of a conventional Al sample holder based on fcc lattice calculations (face centered cubic, i.e. hkl)
.smphldr_peaks <- function(lambda = 1.5406, a = 4.0495) {                         # defaults: lambda = Cu Kα₁, a = lattice parameter for Al
                    hkls <- rbind(c(1,1,1),c(2,0,0),c(2,2,0),c(3,1,1),c(2,2,2))   # 5 Al peaks (@ 2θ = 38.5, 44.7, 65.1, 78.2, 82.4)
                    relI <- DEFAULT_REL_I                                         # relative intensities for Al peaks
                    w    <- relI / max(relI)                                      # normalized peak weights
                    dA   <- a / sqrt(rowSums(hkls^2))                             # d-spacing (in Angstroms)
                    th2  <- 2 * asin(pmin(1, lambda / (2*dA))) * 180/pi           # rad conversion to deg 2θ
                    list(pos = th2, w = w, params = list(lambda = lambda, a = a))
                  }
# pattern template for a conventional aluminum sample holder 
.smphldr_tmpl <- function(tt=NULL, fwhm=0.30, lambda=1.5406, a=4.0495, pos=NULL,
                          shift_deg=0, shape = c("gauss","lorentz"), w = NULL) {
           shape <- match.arg(shape)
                    if (is.null(tt)) {
                       if (exists("dt", inherits = TRUE) && is.list(dt) && 
                          !is.null(dt$tt)) { tt <- dt$tt
                       } else { stop ("Grid not found. Provide a 2θ grid (tt).",
                                      call.=FALSE)
                       }
                    }
                    if (is.null(pos) || is.null(w)) {
               spec <- .smphldr_peaks(lambda = lambda, a = a)
               pos  <- spec$pos
               w    <- spec$w
                    }
        pos_draw <- pos + shift_deg
               M <- .draw_curves(tt, pos_draw, w, fwhm, shape)
             tpl <- rowSums(M) 
           denom <- sqrt(sum(tpl^2, na.rm = TRUE))                                 # L2 normalization (RSE) so alpha is a projection
                    if (is.finite(denom) && denom > 0) tpl <- tpl / denom
                    list(file = "Al_holder", pos = pos, w = w, template = tpl,
                         params = list(lambda = lambda, a = a, fwhm = fwhm,
                                       shape = shape, shift_deg = shift_deg))
                 }
# prepare plot masks (non-destructive)
.smphldr_mask <- function(tt, half_window = 0.25, lambda = 1.5406, a = 4.0495){
                   centers <- .smphldr_peaks(lambda = lambda, a = a)$th2
                   rng <- range(tt, na.rm = TRUE)
                   centers_in <- centers[centers >= rng[1] & centers <= rng[2]]
                   keep <- rep(TRUE, length(tt))
                   if (length(centers_in)) {
                      for (c0 in centers_in) {
                         keep <- keep & !(tt >= (c0 - half_window) & tt <= (c0 + half_window))
                      }
                   }
                   list(keep = keep, centers_in = centers_in,                      # only those visible in current plot range
                        centers_all = centers, half_window = half_window)          # full set is tracked for recipe
                 }
# Recipe label for plots
op_label <- function(op) {
              switch(op$name, "sg" = sprintf("SG(p=%d,n=%d)", op$p, op$n),
                     "wavelet" = sprintf("Wav(%s,%s,L=%d,α=%.1f,≥%d)", 
                                   op$wf, op$mode, op$L, op$alpha, op$level_min),
                     "norm_max1" = "norm_max1",
                     "norm_area1" = "norm_area1",
                     "scale_Imax" = sprintf("Imax→%s", op$target_max), op$name)
            }
recipe_final_label <- function(recipe) {if (!length(recipe$ops)) 
                                           return("Baseline only")
                        op_label(recipe$ops[[length(recipe$ops)]])
                      }
recipe_compact_label <- function(recipe) {if (!length(recipe$ops)) 
                                             return("Baseline only")
                          final <- op_label(recipe$ops[[length(recipe$ops)]])
                          n_prior <- length(recipe$ops) - 1L
                          if (n_prior > 0) sprintf("%s (+%d earlier)",final,n_prior) 
                             else final
                        }
#
#
### ---- Preprocessing script ---- ###
#
#

preprocess_xy <- function() {
       file_path <- .pick_file()
                    if (!nzchar(file_path) || !file.exists(file_path)) {
                       stop("No file selected or path invalid.", call. = FALSE)
                    }
 baseline_method <- "SNIP"
   baseline_args <- list(iterations = 60, log = TRUE)                             # default for log-SNIP
              dt <- read_xy(file_path)
            rawI <- dt$I
           corrI <- do.call(baseline_subtract_vec, c(list(I = rawI, 
                            method = baseline_method), baseline_args))
               y <- corrI                                                         # current spectrum (working vector) + history
   steps_applied <- c(paste0("baseline_", baseline_method))
          recipe <- list(baseline = list(method=baseline_method,args=baseline_args),     
                                         ops = list(),params = list())            # recipe and additional params like lib_step, etc. appended here for completed steps
                    plot_before_after(dt$tt, rawI, corrI)                         # initial preview
#                   
             ####--- Main Menu ---####                   
  repeat {                                                                        # loop menu options:
          choice <- utils::menu(c("End, keep baseline-subtract only)",                            # 1)                                                
                                  "Sample holder correction",                                     # 2)
                                  "Normalize area and intensity to 1",                            # 3)
                                  "Rescale by custom parameters",                                 # 4)
                                  "Denoise (SG / wavelet smoothing)",                             # 5)
                                  "Multiple steps (2–5 in sequence)",                             # 6)
                                  "Save processed spectrum to .xy",                               # 7)
                                  "Save and repeat steps on another file",                        # 8)
                                  "Abort preprocess"),                                            # 9)
                                  title = "\nSpectrum Processing Options:")
          if (choice == 0) break
          switch(choice, {                                                        # 1) stop and exit      
                 return(invisible(list(
                 file = file_path,
                 tt   = dt$tt,
                 raw  = rawI,
                 baseline_corrected = corrI,
                 processed = y,
                 steps = steps_applied,
                 baseline_method = baseline_method,
                 baseline_args   = baseline_args,
                 params  = attr(y, "preproc_params"),
                 lib_grid = if (!is.null(attr(y, "preproc_params")))
                               attr(y, "preproc_params")$lib_grid 
                            else NULL )))
          },
           
          {                                                                        # 2) sample holder subtract submenu
              method <- utils::menu(c("Subtract empty-holder scan",                                  
                                "Subtract aluminum sample holder from template",
                                "Mask only"), title =
                                "Sample Holder Signal Correction")
              if (method == 0) {                                                   # user cancelled
              } else if (method == 1) {                                            # Empty-holder subtraction
                  message("Choose empty holder scan (.xy)")
                  hold_path <- file.choose()
                  if (!nzchar(hold_path) || !file.exists(hold_path)) {
                     message("No valid selection")
                  } else {
           y_before <- y
                  h <- read_xy(hold_path)
                 hI <- baseline_subtract_vec(h$I,method="SNIP",iterations=60)
                       if (!isTRUE(all.equal(h$tt, dt$tt))) {                      # correct grid resolution
                    hI <- approx(h$tt, hI, xout = dt$tt, rule = 2)$y
                       }
              alpha <- .nnls_scalar(hI, y)
                  y <- pmax(0, y - alpha*hI)
      steps_applied <- c(steps_applied, sprintf("Al-holder subtract (α=%.3f)",alpha))
                       plot_before_after(dt$tt,y_before,y,main_before="Before:", 
                           main_after=sprintf("After: holder subtract (α=%.3f)", 
                                              alpha))
                       if (tolower(readline("Keep this change? [y/n]: ")) == "n") {
                     y <- y_before; steps_applied <- head(steps_applied, -1) 
                          message("Sample holder correction reverted")
                       }
                  }
              } else if (method == 2) {                                            # sample holder template subtraction
            y_before <- y                                                    
                tpl0 <- .smphldr_tmpl(tt = dt$tt, shape = "gauss")                 # build conventional Al sample holder template (defaults to Al)
           tpl_entry <- list(file = tpl0$file, pos = tpl0$pos, w = tpl0$w)
           exp_peaks <- extract_peaks_exp(dt$tt, y)                        
                scan <- list(max_shift = 3, shift_step = 0.2, refine = TRUE,       # shift allowances for iterative peak scanning
                             refine_win = 0.05, refine_step = 0.01)
                cfg  <- list(tol = 0.3, heavy_thr = 0.5, heavy_penalty = 1,        # peak scoring params
                             weight_gamma = 1.0)
                best <- best_shift_for_library(exp_peaks, tpl_entry, scan = scan,  # match experimental peaks to template, find best shift
                                               cfg = cfg)    
           tpl_shift <- .smphldr_tmpl(tt = dt$tt, w = tpl0$w, fwhm = 0.2,          # attempt to match template fwhm to experimental peaks
                                      shift_deg = best$shift, pos = tpl0$pos,      # render shifted template
                                      shape = tpl0$params$shape)
         tpl_shifted <- tpl_shift$template                                         # extract peak positions vector
               alpha <- .nnls_scalar(tpl_shifted, y)                               # alpha scalar projects y onto template (0 = no signal)
                   y <- pmax(0, y - alpha * tpl_shifted)
              recipe <- sprintf("Holder subtract: shift=%.3f°, FWHM=%.2f°, α=%.3f", 
                                best$shift, tpl0$params$fwhm, alpha / 100)
       steps_applied <- c(steps_applied, recipe)
         title_after <- paste(strwrap(recipe, width = 48), collapse = "\n")
                        plot_before_after(dt$tt, y_before, y, main_before = 
                                          "Before sample holder subtract:",
                                          main_after = title_after)
                        if (tolower(readline("Keep this change? [Y/n]: "))=="n") {
                      y <- y_before; steps_applied <- head(steps_applied, -1)
                           message("Reverted.")
                        }
              } else if (method == 3) {
                   m <- .smphldr_mask(dt$tt,half_window=0.5,lambda=1.5406,a=4.0495)
  attr(y, "mask_al") <- m
       steps_applied <- c(steps_applied, 
                          if (length(m$centers_in)) sprintf("Al-mask (±%.2f°) at %s",
                              m$half_window,paste(round(m$centers_in,2),collapse = ", "))
                          else sprintf("Al-mask (±%.2f°) [no Al peaks in range]", 
                              m$half_window))
                          plot_before_after(dt$tt, y, y, main_before =             # plot & annotate only visible masks 
                                            "Mask preview", main_after= 
                                            if (length(m$centers_in))
                                               sprintf("Masked %d windows at: %s", 
                                               length(m$centers_in), paste(round(m$centers_in, 2), 
                                               collapse = ", ")) 
                                            else "No Al peaks within this 2θ range")
                                            if (length(m$centers_in)) {
                                               abline(v = m$centers_in, lty = 3)
                                        usr <- par("usr")                          # mask labels
                                      y_top <- usr[4]
                                               text(x = m$centers_in, y = y_top, 
                                                    labels = round(m$centers_in, 2), 
                                                    pos = 3, cex = 0.8)
                                            }
                     }
          },
           
          {                                                                        # 3) Normalize: area=1, max=1
             y <- normalize_max1(y)
             steps_applied <- c(steps_applied, "norm_max1")
             step_est <- median(diff(dt$tt))
             y <- normalize_area1(y, step = step_est)
             steps_applied <- c(steps_applied, "norm_area1")
          },
          
          {                                                                        # 4) Rescale by parameters
             repeat {
               cat("\n--- Set scaling / grid / alignment parameters ---\n")
               def_target_max <- 10000
               def_lib_step   <- 0.02
               def_max_shift  <- 0.20
               def_shift_step <- 0.002
               target_max <- suppressWarnings(as.numeric(readline(
                             sprintf("Target I_max for library scale [default %.0f]: ", 
                                     def_target_max))))
               if (!is.finite(target_max)) target_max <- def_target_max
               lib_step <- suppressWarnings(as.numeric(readline(
                           sprintf("Library step angle, deg [default %.3f]: ", 
                                   def_lib_step))))
               if (!is.finite(lib_step)) lib_step <- def_lib_step
               max_shift <- suppressWarnings(as.numeric(readline(
                            sprintf("Max zero-shift to search, deg [default %.3f]: ", 
                                    def_max_shift))))
               if (!is.finite(max_shift)) max_shift <- def_max_shift
               shift_step <- suppressWarnings(as.numeric(readline(
                             sprintf("Zero-shift step size, deg [default %.3f]: ", 
                                     def_shift_step))))
               if (!is.finite(shift_step)) shift_step <- def_shift_step
               y <- normalize_to_Imax(y, target_max = target_max)
               recipe$ops <- c(recipe$ops, list(list(name = "scale_Imax", 
                                                     target_max = target_max)))
               recipe$params$lib_step   <- lib_step
               recipe$params$max_shift  <- max_shift
               recipe$params$shift_step <- shift_step
               tt_min <- min(dt$tt, na.rm = TRUE)
               tt_max <- max(dt$tt, na.rm = TRUE)
               lib_grid <- build_grid(tt_min, tt_max, lib_step)
               attr(y, "preproc_params") <- list(target_max = target_max,
                                                 lib_step = lib_step,
                                                 max_shift = max_shift,
                                                 shift_step = shift_step,
                                                 lib_grid = lib_grid)
               steps_applied <- c(steps_applied,
                                  sprintf("scale_to_Imax(%.0f)", target_max),
                                  sprintf("lib_step=%.3f", lib_step),
                                  sprintf("shift_win=±%.3f", max_shift),
                                  sprintf("shift_step=%.3f", shift_step))
               cat("\nScaling and parameter setup complete.\n")
               redo <- readline("Press Enter to return to menu, or type 'r' to redo scaling: ")
               if (tolower(redo) != "r") break
             }
         },
           
         {                                                                      # 5) Smoothing (SG / Wavelet / Preview)
             sm_choice <- utils::menu(
               c("Savitzky–Golay (p=2, n=5)",
                 "Savitzky–Golay (p=2, n=7)",
                 "Wavelet: Daubechies db4 (quick)",
                 "Wavelet: Custom settings",
                 "Wavelet: Preview 3 variants (db4)"),
               title = "Choose smoothing"
             )
             if (sm_choice == 0) next
             y_before <- y
             y_new <- NULL
             step_tag <- NULL
             label <- NULL
             meta <- NULL
             if (sm_choice == 1) {                                              # SG (p=2, n=5)
                 y_new <- sg_smooth(y, n = 5, p = 2)
                 step_tag <- "SG(p=2,n=5)"
                 label <- "Savitzky–Golay (p=2, n=5)"
                 meta <- list(name = "sg", p = 2, n = 5)
             } else if (sm_choice == 2) {                                       # SG (p=2, n=7)
                 y_new <- sg_smooth(y, n = 7, p = 2)
                 step_tag <- "SG(p=2,n=7)"
                 label <- "Savitzky–Golay (p=2, n=7)"
                 meta <- list(name = "sg", p = 2, n = 7)
             } else if (sm_choice == 3) {                                       # quick wavelet preset (db4 = waveslim 'la8')
                 y_new <- wavelet_denoise(y, wf = "la8", L = 4, mode = "soft", 
                                          alpha = 1.0, level_min = 1)
                 step_tag <- "Wavelet(db4,soft,L=4,alpha=1.0,level_min=1)"
                 label <- "Wavelet (Daubechies db4), soft, L=4, α=1.0"
                 meta <- list(name = "wavelet", wf = "la8", L = 4, mode = "soft", 
                              alpha = 1.0, level_min = 1)
             } else if (sm_choice == 4) {                                       # custom wavelet
                 cat("\n--- Wavelet (custom) ---\n")
                 wf_in <- readline("Wavelet [db4/la8/haar/coif2...] (default la8): ")
                 if (!nzchar(wf_in)) 
                    wf_in <- "la8"
                    L_in <- suppressWarnings(as.integer(readline("Levels L (3–6 typical) [default 4]: ")))
                 if (!is.finite(L_in)) 
                    L_in <- 4L
                    mode_in <- tolower(readline("Threshold mode [soft/hard] (default soft): "))
                 if (!nzchar(mode_in) || !mode_in %in% c("soft","hard")) 
                    mode_in <- "soft"
                    alpha_in  <- suppressWarnings(as.numeric(readline(
                                 "Threshold multiplier α (1.0 default; >1 stronger): ")))
                 if (!is.finite(alpha_in)) 
                    alpha_in <- 1.0
                    lvlmin_in <- suppressWarnings(as.integer(readline(
                                 "Start thresholding from detail level (1=fine) [default 1]: ")))
                 if (!is.finite(lvlmin_in)) 
                    lvlmin_in <- 1L
                    y_new <- wavelet_denoise(y, wf = wf_in, L = L_in, mode = mode_in,
                                             alpha = alpha_in, level_min = lvlmin_in)
                    step_tag <- sprintf("Wavelet(%s,%s,L=%d,alpha=%.2f,level_min=%d)",
                                        wf_in, mode_in, L_in, alpha_in, lvlmin_in)
                    label <- sprintf("Wavelet (%s, %s, L=%d, α=%.2f, level≥%d)",
                                     wf_in, mode_in, L_in, alpha_in, lvlmin_in)
                    meta <- list(name = "wavelet", wf = wf_in, L = L_in, mode = mode_in,
                                 alpha = alpha_in, level_min = lvlmin_in)
             } else if (sm_choice == 5) {                                       # preview three db4 variants
                 prev <- wavelet_preview3(dt$tt, y, variants = list(
                           list(L = 4, alpha = 1.0),
                           list(L = 5, alpha = 1.5),
                           list(L = 6, alpha = 2.0)),
                           wf = "la8", mode = "soft", level_min = 1)
                 if (is.null(prev)) next                                        # back to submenu
                 y_before <- y
                 y_new    <- prev$y
                 step_tag <- prev$label
                 label    <- prev$label
                 meta     <- if (!is.null(prev$meta)) prev$meta 
                             else list(name="wavelet", wf="la8", L = 
                                       if (!is.null(prev$L)) prev$L 
                                       else 4, mode="soft", alpha = 
                                       if (!is.null(prev$alpha)) prev$alpha 
                                       else 1.0, level_min = 1)
             }
             plot_before_after(dt$tt, y_before, y_new, main_before = "Before smoothing",
                               main_after  = paste("After:", label))
             keep <- tolower(trimws(readline("Keep these smoothing changes? [Y/n]: ")))
             if (keep %in% c("", "y", "yes")) {
                y <- y_new
                steps_applied <- c(steps_applied, step_tag)
                if (!is.null(meta)) 
                   recipe$ops <- c(recipe$ops, list(meta))                      # recipe tracking
                   cat("Changes kept.\n")
             } else {
                 cat("Changes discarded. Keeping previous series.\n")
             }
           },
         
           {                                                                    # 6) batch steps 2–5)
             repeat {
               sub_choice <- utils::menu(c("Normalize max=1","Normalize area=1",
                               "Rescale by factor","Smooth (moving average)","Done"),
                               title = "Batch steps")
               if (sub_choice %in% 1:4) do_step(as.character(sub_choice + 1))
               if (sub_choice == 5) break
               plot_before_after(dt$tt, corrI, y, main_after = paste("After:", 
                                 paste(steps_applied, collapse = ", ")))
             }
           },
           
           {                                                                    # 7) save current version to .xy
             y_out <- y                                                         # default name from current file
             default_out <- if (grepl("\\.xy$", file_path, ignore.case = TRUE)){
                               sub("\\.xy$", "_prepped.xy", file_path, ignore.case = TRUE)
                            } else { paste0(file_path, "_prepped.xy")
                            }
             cat("\nDefault output:\n", default_out, "\n", sep = "")
             out_path <- readline("Enter output path (or press Enter for default): ")
             if (!nzchar(out_path)) out_path <- default_out
             if (file.exists(out_path) && tolower(readline("File exists. Overwrite? [y/N]: "))!= "y") {
                cat("Save cancelled.\n")
             } else { write_xy_file(out_path, dt$tt, y_out) 
                 cat("Wrote: ", out_path, "\n", sep = "")
             }
           },
           
           {                                                                    # 8) save and apply same recipe to a new spectrum
             if (length(recipe$ops) == 0) {
                 cat("No steps recorded — baseline correction only\n")
                 next
             }
             cat("\nSelect a NEW file to process with the same steps...\n")
             new_file <- .pick_file()
             if (!nzchar(new_file) || !file.exists(new_file)) { 
                cat("No file selected.\n"); next 
             }
             dt_new <- read_xy(new_file)                                        # process the new file (baseline + recorded recipe)
             res_new <- apply_recipe_to_data(dt_new, recipe)
             lbl_after <- recipe_compact_label(recipe)                          # label plot 'After:' panel with the last processing step, else use recipe
             plot_before_after(dt_new$tt, res_new$corrI, res_new$y,
                               main_before = "Baseline-corrected (new file)",
                               main_after  = paste("After:", lbl_after))
             out_prepped <- make_out_path(new_file, suffix = "_prepped.xy")     # save current spectrum next to the new file, avoiding overwrite
             write_xy_file(out_prepped, dt_new$tt, res_new$y)
             cat("Processed and saved: ", out_prepped, "\n", sep = "")
             file_path <- new_file                                              # make the new file the current working dataset
             dt <- dt_new
             rawI <- dt$I
             corrI <- res_new$corrI
             y <- res_new$y
             steps_applied <- c(paste0("baseline_", recipe$baseline$method),
                                vapply(recipe$ops, op_label, character(1)))
             cat("Now working on: ", basename(file_path), "\n", sep = "")
           },
           
           {                                                                    # 9) Abort
             stop("Aborted by user.", call. = FALSE)
           }
    )                                                                           # end switch
    
  }                                                                             # end repeat 
  return(invisible(list(file = file_path, tt = dt$tt, raw = rawI,               # user aborted (choice == 0) or fell through
                        baseline_corrected = corrI, processed = y,
                        steps = steps_applied, baseline_method = baseline_method,
                        baseline_args = baseline_args, params = attr(y, "preproc_params"),
                        lib_grid = if (!is.null(attr(y, "preproc_params")))
                                      attr(y, "preproc_params")$lib_grid 
                                   else NULL )))
  }                                                                             # end preprocessing script

#
#
###--- autorun ---###
#
#

preprocess_xy()




