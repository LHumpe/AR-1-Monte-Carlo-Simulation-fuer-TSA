if (!require('matlib')) install.packages('matlib'); library('matlib')

#' Calculation of arithmetic mean
Mean <- function(y) {
  sumOfy <- sum(y)
  lenOfy <- length(y)
  
  return(sumOfy / lenOfy)
}

#' Calculation of empirical variance
Variance <- function(y){
  meanOfy <- Mean(y)
  lenOfy <- length(y)
  ySquared <- c()
  yCentered <- y - meanOfy
  
  for (t in 1:lenOfy){
    ySquared[t] <- yCentered[t] ^ 2
  }
  return((sum(ySquared) / (lenOfy - 1)))
}

#' Calculation of empirical autocovariance
Autocovariance <- function(y, maxlag = 10) {
  meanOfy <- Mean(y)
  lenOfy <- length(y)
  yCentered <- y - meanOfX
  autocovariances <- c()
  
  for (t in 0:maxlag) {
    autocovariances[t+1] <- sum(yCentered * lag(yCentered, -t)) / (lenOfy)
    }

  return(autocovariances)
}

#' Calculation of the empirical autocorrelation
Autocorrelation <- function(y, maxlag = 10, plotting = TRUE) {
  meanOfy <- Mean(y)
  lenOfy <- length(y)
  yCentered <- y - meanOfy
  autocorrelations <- c()
  confidence = 1.96 / sqrt(lenOfy)
  
  for (t in 0:maxlag) {
    autocorrelations[t + 1] <- sum(yCentered * lag(yCentered, - t)) / (lenOfy)
  }
  
  autocorrelations <- autocorrelations / autocorrelations[1]
  
  if (plotting){
    PlotCorrelation(autocorrelations, confidence, maxlag, "ACF")
  }
  
  return(autocorrelations)
}

#' Calculation of partial autocorrelation using Yule-Walker Method
PartialAutocorrelation <- function(y, maxlag = 10, plotting = TRUE) {
  autocor <- Autocorrelation(y, maxlag = maxlag + 1, plotting=FALSE)
  autocorrelations = c(autocor[2])
  
  
  for (j in 2:maxlag){ 
    autocorrelations[k] <- (inv(toeplitz(autocor[1:j])) %*% autocor[2:(j+1)])[k]
  }
  
  if (plotting) {
    lenOfy <- length(y)
    confidence = 1.96 / sqrt(lenOfy)
    PlotCorrelation(partialAutocorrelations, confidence, maxlag, "PACF")
  }

  return(autocorrelations)
}

#' Wrapper function for (P)ACF plots
PlotCorrelation <- function(corrObject, confidence, maxlag, ylabel) {
  yLower = min(-confidence, corrObject) - 0.005
  yUpper = max(confidence, corrObject) + 0.005
  
  if (ylabel == "ACF") {
    xMin = 0
  } else {
    xMin = 1
  }

  plot(corrObject, 
       ylim = c(yLower, yUpper),
       xlim = c(xMin, maxlag),
       ylab = ylabel,
       xlab = "Lag",
       type = "n")
  
  abline(h = c(0, -confidence, confidence),
         col = c("black", "blue", "blue"), 
         lty = c(1, 2 ,2))
  
  segments(xMin:maxlag, 0, xMin:maxlag, corrObject)
}

#' Simulates an autoregressive process 
ar.sim <- function(n, mu, sigma, phi, padding = 100) {
  padding = padding
  n = n + padding
  
  epsilon <- rnorm(n, mean = mu, sd = sigma)
  y = c()
  y[1] <- epsilon[1]
  
  for (t in 2:n){
    y[t] <- y[(t - 1)] * phi + epsilon[t]
  }
  return(ts(y[(padding + 1):n]))
}

#' Monte-Carlo Simulation returning the first autocorrelation 
mc.sim_pacf = function(n, phi, mu, sigma){
  phi.est=vector()
  
  for (i in 1:10000) {
    y = ar.sim(n, mu, sigma, phi)
    phi.est[i] = Autocorrelation(y, maxlag = 2, plotting = FALSE)[2]
  }
  
  se = sqrt(Variance(phi.est))
  
  return(list(mean(phi.est), mean(phi.est) - phi, se))
}
