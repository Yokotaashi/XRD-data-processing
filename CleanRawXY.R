# ===== preprocess_xy_interactive.R =====
#
# ---- dependency check ----
#
required_pkgs <- c("data.table","baseline","prospectr","limSolve","signal","waveslim")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing)) {
  stop("Required packages missing: ", paste(missing, collapse=", "),
       "\nInstall with: install.packages(c(", 
       paste(sprintf('\"%s\"', missing), collapse=", "), "))",
       call. = FALSE)
}
# only data.table attached; others called directly with ::
library(data.table)
#
# functions
#
# file picker
.pick_file <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::selectFile(caption = "Choose experimental .xy file",
                                path = getwd(),
                                filter = "XY files (*.xy)|*.xy")
    if (!is.null(p) && nzchar(p)) return(p)
  }
  message("Using base file chooser (RStudio API not available).")
  suppressWarnings(file.choose(new = FALSE))
}
# file reader
read_xy <- function(path) {
  dt <- data.table::fread(path, header = FALSE, col.names = c("tt","I"))
  dt[order(tt)]
}
normalize_max1  <- function(y) { m <- max(y, na.rm=TRUE); if (!is.finite(m) || m==0) y else y/m }
normalize_area1 <- function(y) { s <- sum(y, na.rm=TRUE); if (!is.finite(s) || s==0) y else y/s }
moving_average <- function(y, k = 3) {
  if (k <= 1) return(y)
  filt <- rep(1/k, k)
  z <- as.numeric(stats::filter(y, filt, sides = 2))
  z[is.na(z)] <- 0
  z
}
plot_before_after <- function(tt, rawI, corrI, main_before="Raw", main_after="Baseline subtracted") {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1,2), mar = c(4,4,2,1))
  plot(tt, rawI, type="l", xlab=expression(2*theta~"(deg)"), ylab="Intensity", main=main_before)
  plot(tt, corrI, type="l", xlab=expression(2*theta~"(deg)"), ylab="Intensity", main=main_after)
}
# file writer
write_xy_file <- function(path, x, y, digits_x = 5, digits_y = 6) {
  # two-column txt, no headers
  fmt <- paste0("%.", digits_x, "f %.", digits_y, "f")
  con <- file(path, "wt"); on.exit(close(con), add = TRUE)
  for (i in seq_along(x)) writeLines(sprintf(fmt, x[i], y[i]), con)
  invisible(path)
}

# ---- log-SNIP with optional auto offset ----
snip_baseline <- function(y, iterations = 60, log = FALSE, offset = NULL) {
  y <- as.numeric(y)
  n <- length(y)
  if (!is.finite(iterations) || iterations < 1) return(rep(0, n))
  
  if (log) {
    # choose a stable positive offset if not supplied
    if (is.null(offset)) {
      yy <- y[is.finite(y) & y > 0]
      offset <- if (length(yy)) max(1e-6, 0.01 * median(yy)) else 1e-6
    }
    z <- log(pmax(y, 0) + offset)
    b <- z
    for (k in 1:iterations) {
      zl <- c(rep(b[1], k), b[1:(n - k)])
      zr <- c(b[(k + 1):n], rep(b[n], k))
      avg <- 0.5 * (zl + zr)
      b <- pmin(b, avg)
    }
    bl <- pmax(exp(b) - offset, 0)   # back-transform baseline
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

# scaling
normalize_to_Imax <- function(y, target_max = 10000) {
  m <- max(y, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(y)
  y * (target_max / m)
}

build_grid <- function(tt_min, tt_max, step_deg) {
  if (!is.finite(step_deg) || step_deg <= 0) stop("Library step must be > 0")
  seq(from = tt_min, to = tt_max, by = step_deg)
}

# ---- unified subtractor (supports SNIP with log/offset) ----
baseline_subtract_vec <- function(I, method = "rollingBall", ...) {
  y <- as.numeric(I)
  if (toupper(method) == "SNIP") {
    dots <- list(...)
    iters  <- if (!is.null(dots$iterations)) as.integer(dots$iterations) else 60L
    uselog <- isTRUE(dots$log)              # default FALSE unless set TRUE
    off    <- if (!is.null(dots$offset)) dots$offset else NULL
    bl <- snip_baseline(y, iterations = iters, log = uselog, offset = off)
    corr <- y - bl
  } else {
    m <- matrix(y, nrow = 1)
    blfit <- baseline::baseline(m, method = method, ...)
    corr <- as.numeric(baseline::getCorrected(blfit))
  }
  corr[corr < 0] <- 0
  corr
}

# --- Savitzky–Golay smoothing (poly order p, window n must be odd) ---
sg_smooth <- function(y, n = 5, p = 2) {
  if (n %% 2 == 0) n <- n + 1          # make sure n is odd
  if (n <= p) n <- p + 3 + (p %% 2)    # ensure window > poly order
  signal::sgolayfilt(y, p = p, n = n)
}

# --- Wavelet denoising (configurable) ---
# wf: waveslim wavelet name ("la8" = Daubechies db4; "la10" ~ db5; "la12" ~ db6)
# L:  number of decomposition levels (higher = more aggressive, typical 3–6)
# mode: "soft" (default) or "hard" thresholding
# alpha: threshold multiplier (1.0 default; >1.0 = more aggressive)
# level_min: start thresholding from this detail level (1 = finest details)
wavelet_denoise <- function(y, wf = "la8", L = 4, mode = c("soft","hard"),
                            alpha = 1.0, level_min = 1) {
  mode <- match.arg(mode)
  n <- length(y)
  pow2 <- 2^ceiling(log2(n))
  pad_len <- pow2 - n
  ypad <- if (pad_len > 0) c(y, rep(tail(y, 1), pad_len)) else y
  
  # ensure wf is valid (prevents wave.filter error)
  wf <- resolve_wavelet_name(wf)
  # clamp L to allowable range
  L <- clamp_levels(L, length(ypad))
  
  d <- waveslim::dwt(ypad, wf = wf, n.levels = L, boundary = "periodic")
  
  sigma_hat <- function(w) stats::mad(w, center = 0) / 0.6745
  shrink <- function(w, thr, mode) if (mode == "soft") sign(w)*pmax(abs(w)-thr,0) else w*(abs(w)>thr)
  
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

# ---- Preview 3 wavelet variants (Daubechies db4 = waveslim "la8") ----
wavelet_preview3 <- function(x, y,
                             variants = list(
                               list(L = 4, alpha = 1.0),
                               list(L = 5, alpha = 1.5),
                               list(L = 6, alpha = 2.0)
                             ),
                             wf = "la8", mode = "soft", level_min = 1) {
  # compute three candidates
  ys <- lapply(variants, function(v)
    wavelet_denoise(y, wf = wf, L = v$L, mode = mode,
                    alpha = v$alpha, level_min = level_min))
  
  # labels (coerce types explicitly to satisfy vapply)
  Ls     <- vapply(variants, function(v) as.integer(v[["L"]]),    integer(1))
  alphas <- vapply(variants, function(v) as.numeric(v[["alpha"]]), numeric(1))
  labs <- sprintf("db4 soft  L=%d  \u03B1=%.1f", Ls, alphas)
  
  # 2x2 plot: raw + 3 variants
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(2,2), mar = c(4,4,2.5,1))
  plot(x, y, type = "l", main = "Raw (before)",
       xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
  plot(x, ys[[1]], type = "l", main = labs[1],
       xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
  plot(x, ys[[2]], type = "l", main = labs[2],
       xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
  plot(x, ys[[3]], type = "l", main = labs[3],
       xlab = expression(2*theta*" (deg)"), ylab = "Intensity")
  
  # ask which one to apply
  cat("\nPreviewed variants:\n  [1]", labs[1], "\n  [2]", labs[2],
      "\n  [3]", labs[3], "\n  [0] Cancel\n")
  pick <- suppressWarnings(as.integer(readline("Apply which variant? (0/1/2/3): ")))
  if (!is.finite(pick) || pick == 0) return(NULL)
  if (pick %in% 1:3) {
    return(list(
      y = ys[[pick]],
      label = sprintf("Wavelet(db4,soft,L=%d,alpha=%.1f,level_min=%d)",
                      Ls[pick], alphas[pick], level_min)
    ))
  }
  NULL
}

# Map friendly names -> waveslim codes and validate
resolve_wavelet_name <- function(wf_in) {
  aliases <- c(
    "db2"="d4",     # Daubechies-2 (a.k.a. d4)
    "db4"="la8",    # Daubechies-4
    "db5"="la10",   # Daubechies-5
    "db6"="la12",   # Daubechies-6
    "db7"="la14",
    "db8"="la16",
    "haar"="haar"
  )
  valid <- c("haar","d4","la8","la10","la12","la14","la16",
             "coif1","coif2","coif3","bl14","bl18","bl20","mb4","mb8",
             "fk4","fk6","fk8","fk14","fk18","fk22")
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
  pow2 <- 2^ceiling(log2(n))          # padded length used in our dwt
  maxL <- floor(log2(pow2))           # maximum allowable levels
  Lc <- max(1L, min(as.integer(L), maxL))
  if (Lc != L) message("Adjusted levels L from ", L, " to ", Lc, " (signal length constraint).")
  Lc
}

#
# tracking preprocessing recipe
apply_recipe_to_file <- function(file_path, recipe, save_baseline = FALSE, suffix = "_prepped.xy") {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  dt <- read_xy(file_path)
  rawI <- dt$I
  
  # baseline from the recipe
  base <- recipe$baseline
  corrI <- do.call(baseline_subtract_vec, c(list(I = rawI, method = base$method), base$args))
  y <- corrI
  
  # replay steps
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
  
  # output paths
  if (grepl("\\.xy$", file_path, ignore.case = TRUE)) {
    out_prepped  <- sub("\\.xy$", paste0(suffix), file_path, ignore.case = TRUE)
    out_baseline <- sub("\\.xy$", "_baseline.xy", file_path, ignore.case = TRUE)
  } else {
    out_prepped  <- paste0(file_path, suffix)
    out_baseline <- paste0(file_path, "_baseline.xy")
  }
  
  # write
  write_xy_file(out_prepped,  dt$tt, y)
  if (isTRUE(save_baseline)) write_xy_file(out_baseline, dt$tt, corrI)
  
  invisible(list(prepped = out_prepped, baseline = if (save_baseline) out_baseline else NULL))
}

# Output path next to input file
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

# Recipe label for plots
op_label <- function(op) {
  switch(op$name,
         "sg"         = sprintf("SG(p=%d,n=%d)", op$p, op$n),
         "wavelet"    = sprintf("Wav(%s,%s,L=%d,α=%.1f,≥%d)", op$wf, op$mode, op$L, op$alpha, op$level_min),
         "norm_max1"  = "norm_max1",
         "norm_area1" = "norm_area1",
         "scale_Imax" = sprintf("Imax→%s", op$target_max),
         op$name
  )
}
# Short label for the final step only
recipe_final_label <- function(recipe) {
  if (!length(recipe$ops)) return("Baseline only")
  op_label(recipe$ops[[length(recipe$ops)]])
}

# Optional: compact summary like "Wavelet(db4…) (+2 earlier)"
recipe_compact_label <- function(recipe) {
  if (!length(recipe$ops)) return("Baseline only")
  final <- op_label(recipe$ops[[length(recipe$ops)]])
  n_prior <- length(recipe$ops) - 1L
  if (n_prior > 0) sprintf("%s  (+%d earlier)", final, n_prior) else final
}


# ---- preprocessing script----
#

preprocess_xy_interactive <- function() {
  file_path <- .pick_file()
  if (!nzchar(file_path) || !file.exists(file_path)) {
    stop("No file selected or path invalid.", call. = FALSE)
  }
  
  # defaults
  baseline_method <- "SNIP"
  baseline_args   <- list(iterations = 60, log = TRUE)  # <- if you’re using log-SNIP
  
  # load + baseline
  dt   <- read_xy(file_path)
  rawI <- dt$I
  corrI <- do.call(baseline_subtract_vec,
                   c(list(I = rawI, method = baseline_method), baseline_args))
  
  # working vector + history
  y <- corrI
  steps_applied <- c(paste0("baseline_", baseline_method))
  recipe <- list(
    baseline = list(method = baseline_method, args = baseline_args),
    ops = list(),          # we’ll append here as steps are KEPT
    params = list()        # optional bag for things like lib_step, etc.
  )
  
  # initial preview
  plot_before_after(dt$tt, rawI, corrI)
  
  # ===== OUTER MENU LOOP =====
  repeat {
    choice <- utils::menu(
      c("Stop (return baseline-subtracted only)",  # 1
        "Normalize max = 1",                       # 2
        "Normalize area = 1",                      # 3
        "Rescale by parameters",                   # 4
        "Smooth (SG / Wavelet)",                   # 5
        "Multiple steps (2–5 in sequence)",        # 6
        "Save processed to .xy",                   # 7
        "Apply same steps to another file (auto-save)",  # 8 
        "Abort"),                                  # 9
      title = "\nAction options:"
    )
    if (choice == 0) break
    
    switch(choice,
           { # 1) Stop -> return and exit
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
                 attr(y, "preproc_params")$lib_grid else NULL
             )))
           },
           
           { # 2) Normalize max=1
             y <- normalize_max1(y)
             steps_applied <- c(steps_applied, "norm_max1")
           },
           
           { # 3) Normalize area=1
             step_est <- median(diff(dt$tt))
             y <- normalize_area1(y, step = step_est)
             steps_applied <- c(steps_applied, "norm_area1")
           },
           
           { # 4) Rescale by parameters (your existing Option 4 block goes here unchanged)
             repeat {
               cat("\n--- Set scaling / grid / alignment parameters ---\n")
               def_target_max <- 10000
               def_lib_step   <- 0.02
               def_max_shift  <- 0.20
               def_shift_step <- 0.002
               
               target_max <- suppressWarnings(as.numeric(readline(
                 sprintf("Target I_max for library scale [default %.0f]: ", def_target_max))))
               if (!is.finite(target_max)) target_max <- def_target_max
               
               lib_step <- suppressWarnings(as.numeric(readline(
                 sprintf("Library step angle, deg [default %.3f]: ", def_lib_step))))
               if (!is.finite(lib_step)) lib_step <- def_lib_step
               
               max_shift <- suppressWarnings(as.numeric(readline(
                 sprintf("Max zero-shift to search, deg [default %.3f]: ", def_max_shift))))
               if (!is.finite(max_shift)) max_shift <- def_max_shift
               
               shift_step <- suppressWarnings(as.numeric(readline(
                 sprintf("Zero-shift step size, deg [default %.3f]: ", def_shift_step))))
               if (!is.finite(shift_step)) shift_step <- def_shift_step
               
               y <- normalize_to_Imax(y, target_max = target_max)
               
               recipe$ops <- c(recipe$ops, list(list(name = "scale_Imax", target_max = target_max)))
               recipe$params$lib_step  <- lib_step
               recipe$params$max_shift <- max_shift
               recipe$params$shift_step<- shift_step
               
               tt_min <- min(dt$tt, na.rm = TRUE)
               tt_max <- max(dt$tt, na.rm = TRUE)
               lib_grid <- build_grid(tt_min, tt_max, lib_step)
               
               attr(y, "preproc_params") <- list(
                 target_max  = target_max,
                 lib_step    = lib_step,
                 max_shift   = max_shift,
                 shift_step  = shift_step,
                 lib_grid    = lib_grid
               )
               
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
           
           { # 5) Smoothing (SG / Wavelet / Preview)
             sm_choice <- utils::menu(
               c("Savitzky–Golay (p=2, n=5)",
                 "Savitzky–Golay (p=2, n=7)",
                 "Wavelet: Daubechies db4 (quick)",
                 "Wavelet: Custom settings",
                 "Wavelet: Preview 3 variants (db4)"),
               title = "Choose smoothing"
             )
             if (sm_choice == 0) next  # back to main menu
             
             y_before <- y
             y_new <- NULL
             step_tag <- NULL
             label <- NULL
             meta <- NULL
             
             if (sm_choice == 1) {
               # SG (p=2, n=5)
               y_new <- sg_smooth(y, n = 5, p = 2)
               step_tag <- "SG(p=2,n=5)"
               label <- "Savitzky–Golay (p=2, n=5)"
               meta <- list(name = "sg", p = 2, n = 5)
               
             } else if (sm_choice == 2) {
               # SG (p=2, n=7)
               y_new <- sg_smooth(y, n = 7, p = 2)
               step_tag <- "SG(p=2,n=7)"
               label <- "Savitzky–Golay (p=2, n=7)"
               meta <- list(name = "sg", p = 2, n = 7)
               
             } else if (sm_choice == 3) {
               # Quick wavelet preset (db4 ≈ waveslim 'la8')
               y_new <- wavelet_denoise(y, wf = "la8", L = 4, mode = "soft", alpha = 1.0, level_min = 1)
               step_tag <- "Wavelet(db4,soft,L=4,alpha=1.0,level_min=1)"
               label <- "Wavelet (Daubechies db4), soft, L=4, α=1.0"
               meta <- list(name = "wavelet", wf = "la8", L = 4, mode = "soft", alpha = 1.0, level_min = 1)
               
             } else if (sm_choice == 4) {
               # Custom wavelet (wavelet_denoise() resolves/validates wf internally)
               cat("\n--- Wavelet (custom) ---\n")
               wf_in <- readline("Wavelet [db4/la8/haar/coif2...] (default la8): ")
               if (!nzchar(wf_in)) wf_in <- "la8"
               L_in <- suppressWarnings(as.integer(readline("Levels L (3–6 typical) [default 4]: ")))
               if (!is.finite(L_in)) L_in <- 4L
               mode_in <- tolower(readline("Threshold mode [soft/hard] (default soft): "))
               if (!nzchar(mode_in) || !mode_in %in% c("soft","hard")) mode_in <- "soft"
               alpha_in  <- suppressWarnings(as.numeric(readline("Threshold multiplier α (1.0 default; >1 stronger): ")))
               if (!is.finite(alpha_in)) alpha_in <- 1.0
               lvlmin_in <- suppressWarnings(as.integer(readline("Start thresholding from detail level (1=fine) [default 1]: ")))
               if (!is.finite(lvlmin_in)) lvlmin_in <- 1L
               
               y_new <- wavelet_denoise(y, wf = wf_in, L = L_in, mode = mode_in,
                                        alpha = alpha_in, level_min = lvlmin_in)
               step_tag <- sprintf("Wavelet(%s,%s,L=%d,alpha=%.2f,level_min=%d)",
                                   wf_in, mode_in, L_in, alpha_in, lvlmin_in)
               label <- sprintf("Wavelet (%s, %s, L=%d, α=%.2f, level≥%d)",
                                wf_in, mode_in, L_in, alpha_in, lvlmin_in)
               meta <- list(name = "wavelet", wf = wf_in, L = L_in, mode = mode_in,
                            alpha = alpha_in, level_min = lvlmin_in)
               
             } else if (sm_choice == 5) {
               # Preview three db4 variants
               prev <- wavelet_preview3(
                 dt$tt, y,
                 variants = list(
                   list(L = 4, alpha = 1.0),
                   list(L = 5, alpha = 1.5),
                   list(L = 6, alpha = 2.0)
                 ),
                 wf = "la8", mode = "soft", level_min = 1
               )
               if (is.null(prev)) next  # back to menu
               
               y_before <- y
               y_new    <- prev$y
               step_tag <- prev$label
               label    <- prev$label
               meta <- if (!is.null(prev$meta)) prev$meta else list(name="wavelet", wf="la8",
                                                                    L = if (!is.null(prev$L)) prev$L else 4,
                                                                    mode="soft",
                                                                    alpha = if (!is.null(prev$alpha)) prev$alpha else 1.0,
                                                                    level_min = 1)
             }
             plot_before_after(dt$tt, y_before, y_new,
                               main_before = "Before smoothing",
                               main_after  = paste("After:", label)
             )
               
             keep <- tolower(trimws(readline("Keep these smoothing changes? [Y/n]: ")))
             if (keep %in% c("", "y", "yes")) {
               y <- y_new
               steps_applied <- c(steps_applied, step_tag)
               if (!is.null(meta)) recipe$ops <- c(recipe$ops, list(meta))  # <-- the only place we log
               cat("Changes kept.\n")
             } else {
               cat("Changes discarded. Keeping previous series.\n")
             }
           },
           
           { # 6) Multiple steps (batch of 2–5)
             repeat {
               sub_choice <- utils::menu(
                 c("Normalize max=1","Normalize area=1","Rescale by factor","Smooth (moving average)","Done"),
                 title = "Batch steps"
               )
               if (sub_choice %in% 1:4) do_step(as.character(sub_choice + 1))
               if (sub_choice == 5) break
               plot_before_after(dt$tt, corrI, y,
                                 main_after = paste("After:", paste(steps_applied, collapse = ", ")))
             }
           },
           
           { # 7) Save processed to .xy (current series only)
             y_out <- y
             
             # default name from current file
             default_out <- if (grepl("\\.xy$", file_path, ignore.case = TRUE)) {
               sub("\\.xy$", "_prepped.xy", file_path, ignore.case = TRUE)
             } else {
               paste0(file_path, "_prepped.xy")
             }
             
             cat("\nDefault output:\n", default_out, "\n", sep = "")
             out_path <- readline("Enter output path (or press Enter for default): ")
             if (!nzchar(out_path)) out_path <- default_out
             
             if (file.exists(out_path) && tolower(readline("File exists. Overwrite? [y/N]: ")) != "y") {
               cat("Save cancelled.\n")
             } else {
               write_xy_file(out_path, dt$tt, y_out)
               cat("Wrote: ", out_path, "\n", sep = "")
             }
           },
           
           { # 8) Apply same steps to another file (auto-save using new file name)
             if (length(recipe$ops) == 0) {
               cat("No steps recorded beyond baseline — nothing to replay.\n")
               next
             }
             cat("\nSelect a NEW file to process with the recorded steps...\n")
             new_file <- .pick_file()
             if (!nzchar(new_file) || !file.exists(new_file)) { cat("No file selected.\n"); next }
             
             # Process the new file (baseline + recorded ops)
             dt_new <- read_xy(new_file)
             res_new <- apply_recipe_to_data(dt_new, recipe)
             
             # Label the 'After:' panel by the FINAL step, not the first one
             lbl_after <- recipe_compact_label(recipe)  # or use recipe_final_label(recipe)
             
             # Show Before/After on the NEW file
             plot_before_after(
               dt_new$tt,
               res_new$corrI,
               res_new$y,
               main_before = "Baseline-corrected (new file)",
               main_after  = paste("After:", lbl_after)
             )
             
             # Auto-save next to the NEW file, avoiding overwrite
             out_prepped <- make_out_path(new_file, suffix = "_prepped.xy")
             write_xy_file(out_prepped, dt_new$tt, res_new$y)
             cat("Processed and saved: ", out_prepped, "\n", sep = "")
             
             # Optional: also save baseline-only for the new file
             # out_base <- make_out_path(new_file, suffix = "_baseline.xy")
             # write_xy_file(out_base, dt_new$tt, res_new$corrI)
             # cat("Baseline-only saved: ", out_base, "\n", sep = "")
             
             # Optional: make the NEW file the current working dataset
             file_path <- new_file
             dt <- dt_new
             rawI <- dt$I
             corrI <- res_new$corrI
             y <- res_new$y
             steps_applied <- c(
               paste0("baseline_", recipe$baseline$method),
               vapply(recipe$ops, op_label, character(1))
             )
             cat("Now working on: ", basename(file_path), "\n", sep = "")
           },
           
            { # 9) Abort
             stop("Aborted by user.", call. = FALSE)
           }
    ) # end switch
    
  } # end repeat
    # If we get here, the user closed the menu (choice == 0) or fell through.
  return(invisible(list(
    file = file_path,
    tt = dt$tt,
    raw = rawI,
    baseline_corrected = corrI,
    processed = y,
    steps = steps_applied,
    baseline_method = baseline_method,
    baseline_args = baseline_args,
    params = attr(y, "preproc_params"),
    lib_grid = if (!is.null(attr(y, "preproc_params")))
      attr(y, "preproc_params")$lib_grid else NULL
  )))
} # end of preprocess

# autorun
preprocess_xy_interactive()