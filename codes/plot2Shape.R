plot2Shape <- function (shapeMatrix, background = NULL, colDots = rgb(0, 0, 
                                                                      1, 0.1), colDotsBg = rgb(0, 0, 0, 0.1), colLine = "steelblue", 
                        colLineBg = "gray50", cex = 0.5, lwd = 4, ylim, ...) 
{
  n <- nrow(shapeMatrix)
  mu <- colMeans(shapeMatrix, na.rm = TRUE)
  m <- length(mu)
  span <- round(m/2)
  if (is.null(background)) {
    if (missing(ylim)) 
      ylim <- range(mu, na.rm = TRUE)
    plot(mu, col = colDots, pch = 19, cex = cex, xaxt = "n", 
         xlab = "", #ylab = paste0("Mean value (n=", n, ")"), 
         ylim = ylim, ...)
    # axis(1, at = c(0, m), labels = c(-span, paste0("+", span)))
    # abline(v = span, lty = 2, col = "gray30")
    lines(lowess(mu, f = 1/10), col = colLine, lwd = lwd)
  }
  else {
    mu1 <- mu
    mu2 <- colMeans(background, na.rm = TRUE)
    if (missing(ylim)) 
      ylim <- range(mu1, mu2, na.rm = TRUE)
    plot(mu1, col = colDots, pch = 19, cex = cex, xaxt = "n", 
         xlab = "", #ylab = paste0("Mean value (n=", n, ")"), 
         ylim = ylim, ...)
    points(mu2, pch = 19, cex = cex, col = colDotsBg)
    # axis(1, at = c(0, m), labels = c(-span, paste0("+", span)))
    # abline(v = span, lty = 2, col = "gray30")
    lines(lowess(mu1, f = 1/10), col = colLine, lwd = lwd)
    lines(lowess(mu2, f = 1/10), col = colLineBg, lwd = lwd)
  }
}
