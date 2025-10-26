# pick two xy files to compare
f1 <- file.choose()
f2 <- file.choose()
d1 <- read_xy(f1)   # expects columns: tt, I
d2 <- read_xy(f2)
op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(1, 2), mar = c(4,4,2,1))

# use a common y-scale so intensities are comparable and plot
yl <- range(d1$I, d2$I, na.rm = TRUE)
plot(d1$tt, d1$I, type = "l", ylim = yl,
     xlab = expression(2*theta*" (deg)"), ylab = "Intensity",
     main = basename(f1))
plot(d2$tt, d2$I, type = "l", ylim = yl,
     xlab = expression(2*theta*" (deg)"), ylab = "Intensity",
     main = basename(f2))
