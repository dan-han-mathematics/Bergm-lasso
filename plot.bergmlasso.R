plot.bergmlasso <- function(x, ...) {
  
  seqq <- 4
  par(mfrow = c(min(4, x$dim), 3), 
      oma   = c(0, 0, 3, 0), 
      mar   = c(4, 3, 1.5, 1))
  
  for (i in 1:x$dim) {
    if (i %in% c(5, 9, 13)) {
      dev.new()
      par(mfrow = c(min(4, x$dim - (i - 1)), 3), 
          oma = c(0, 0, 3, 0), 
          mar = c(4, 3, 1.5, 1))
    }
    plot(density(x$Theta[, i]), 
         main = "", 
         axes = FALSE, 
         xlab = bquote(paste(theta[.(i)], " (", .(x$specs[i]), ")")),
         ylab = "", lwd = 2)
    axis(1)
    axis(2)
    traceplot(x$Theta[, i], type = "l", xlab = "Iterations", ylab = "")
    autocorr.plot(x$Theta[, i],lag.max=1000, auto.layout = FALSE, ...)
    if (x$dim > 4) seqq <- seq(4, x$dim, 4)
    if (i %in% union(x$dim, seqq)) {
      title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
    }
  }
}