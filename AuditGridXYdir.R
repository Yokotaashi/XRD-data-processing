#
# functions:
#
read_xy <- function(path) {
  dt <- fread(path, header = FALSE, col.names = c("tt","I"))
  dt[order(tt)]
}

# Quick grid stats for a single file
.xy_grid_stats <- function(dt) {
  dtt <- diff(dt$tt)
  list(
    n = nrow(dt),
    tt_min = min(dt$tt), tt_max = max(dt$tt),
    step_median = as.numeric(median(dtt)),
    step_iqr = as.numeric(IQR(dtt)),
    uniform = isTRUE(all(abs(dtt - median(dtt)) < 1e-8))
  )
}

# 1) Audit the whole library
audit_xy_library <- function(lib_dir, pattern = "\\.xy$") {
  files <- list.files(lib_dir, pattern = pattern, full.names = TRUE)
  stopifnot(length(files) > 0)
  res <- rbindlist(lapply(files, function(f) {
    dt <- read_xy(f)
    gs <- .xy_grid_stats(dt)
    data.table(
      file = basename(f),
      n = gs$n,
      tt_min = gs$tt_min,
      tt_max = gs$tt_max,
      step_median = gs$step_median,
      step_iqr = gs$step_iqr,
      uniform = gs$uniform,
      I_max = max(dt$I, na.rm = TRUE),
      I_area = sum(pmax(dt$I,0), na.rm = TRUE) * gs$step_median
    )
  }))
  setorder(res, file)
  res[, scaled_max1 := abs(I_max - 1) < 1e-6]       # are refs already max-normalized?
  res[, scaled_area1 := abs(I_area - 1) < 1e-3]     # (only true if you chose area=1)
  res
}

# 2) Decide on a canonical grid (either adopt the first file's grid, or form a consensus)
canonical_grid <- function(lib_dir, method = c("first","consensus")) {
  method <- match.arg(method)
  files <- list.files(lib_dir, pattern = "\\.xy$", full.names = TRUE)
  stopifnot(length(files) > 0)
  dts <- lapply(files, read_xy)
  if (method == "first") {
    return(dts[[1]]$tt)
  }
# consensus: use median step and unified range across all files in this dir
  tt_min <- max(sapply(dts, function(x) min(x$tt)))
  tt_max <- min(sapply(dts, function(x) max(x$tt)))
  if (tt_min >= tt_max) stop("No overlapping 2θ range across library.")
  step_med <- median(unlist(lapply(dts, function(x) diff(x$tt))), na.rm = TRUE)
  step_med <- as.numeric(format(signif(step_med, 3), scientific = FALSE))
  seq(from = tt_min, to = tt_max, by = step_med)
}

# 3) Simple resampler + normalizers
resample_to_grid <- function(tt, I, grid) {
  stats::approxfun(tt, I, rule = 2)(grid)
}
normalize_max1  <- function(y) { m <- max(y, na.rm=TRUE); if (!is.finite(m)||m==0) y else y/m }
normalize_area1 <- function(y, step) { a <- sum(y, na.rm=TRUE)*step; if (!is.finite(a)||a==0) y else y/a }

#
## run lines:
#
lib_dir <- "XRD_library"  # <- folder (in default directory)
audit <- audit_xy_library(lib_dir)
View(audit)  # or print(audit)
