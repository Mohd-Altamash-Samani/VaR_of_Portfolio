library(quantmod)

start <- format(as.Date("2018-01-01"), "%Y-%m-%d")
end <- format(as.Date("2019-12-31"), "%Y-%m-%d")

# Load stock data
TCS <- getSymbols('TCS.NS', src = "yahoo", from = start, to = end, auto.assign = FALSE)
HDFCBANK <- getSymbols('HDFCBANK.NS', src = "yahoo", from = start, to = end, auto.assign = FALSE)
RELIANCE <- getSymbols('RELIANCE.NS', src = "yahoo", from = start, to = end, auto.assign = FALSE)
SUNPHARMA <- getSymbols('SUNPHARMA.NS', src = "yahoo", from = start, to = end, auto.assign = FALSE)
TATAMOTORS <- getSymbols('TATAMOTORS.NS', src = "yahoo", from = start, to = end, auto.assign = FALSE)

# Plot stock prices
par(mfrow = c(2, 3))
ts.plot(TCS$TCS.NS.Adjusted, main = "TCS Stock Price", ylab = "TCS Stock Price")
ts.plot(HDFCBANK$HDFCBANK.NS.Adjusted, main = "HDFCBANK Stock Price", ylab = "HDFCBANK Stock Price")
ts.plot(RELIANCE$RELIANCE.NS.Adjusted, main = "RELIANCE Stock Price", ylab = "RELIANCE Stock Price")
ts.plot(SUNPHARMA$SUNPHARMA.NS.Adjusted, main = "SUNPHARMA Stock Price", ylab = "SUNPHARMA Stock Price")
ts.plot(TATAMOTORS$TATAMOTORS.NS.Adjusted, main = "TATAMOTORS Stock Price", ylab = "TATAMOTORS Stock Price")
par(mfrow = c(1, 1))

library(quantmod)
library(QRM)

# Convert daily returns to vectors
a <- as.vector(dailyReturn(TCS))
b <- as.vector(dailyReturn(HDFCBANK))
c <- as.vector(dailyReturn(RELIANCE))
d <- as.vector(dailyReturn(SUNPHARMA))
e <- as.vector(dailyReturn(TATAMOTORS))

# Combine data into a matrix
data <- cbind(a, b, c, d, e)
colnames(data) <- c("TCS", "HDFCBANK", "RELIANCE", "SUNPHARMA", "TATAMOTORS")

# Initialize vectors for parameters
nu <- mu <- sigma <- numeric(ncol(data))
names <- c("TCS", "HDFCBANK", "RELIANCE", "SUNPHARMA", "TATAMOTORS")

# Fit t-distribution to each stock
par(mfrow=c(2,3))
for (i in 1:ncol(data)) {
  tfit <- fit.st(data[, i])  # Fit Student's t-distribution
  print(tfit)  # Check the structure
  
  # Ensure estimates exist and have correct length
  if (!is.null(tfit$par.ests) && length(tfit$par.ests) == 3) {
    nu[i] <- tfit$par.ests[1]
    mu[i] <- tfit$par.ests[2]
    sigma[i] <- tfit$par.ests[3]
    cat("✅ Fitting successful for", names[i], "\n")
  } else {
    cat("⚠️ Warning: fit.st() failed for", names[i], "\n")
  }
  
  hist(data[, i], nclass=20, probability=TRUE, ylim=c(0,35), 
       main=paste(names[i], "Returns Histogram"), 
       xlab=paste(names[i], "Returns"))
  
  curve(dnorm(x, mean=mean(data[, i]), sd=sd(data[, i])), col="red", lwd=3, add=TRUE)
  curve(dt((x - mu[i]) / sigma[i], df=nu[i]) / sigma[i], col="blue", lwd=3, add=TRUE)
}
par(mfrow=c(1,1))

# Transform data using empirical cumulative distribution function (ECDF)
datatrans <- data
for (i in 1:ncol(data)) {
  datatrans[, i] <- ecdf(data[, i])(data[, i])
}

par(mfrow=c(2,3))
for (i in 1:(ncol(datatrans)-1)) {
  for (j in (i+1):ncol(datatrans)) {
    plot(datatrans[,i], datatrans[,j], 
         xlab=paste(names[i],"Returns ECDF"),
         ylab=paste(names[j],"Returns ECDF"),
         main=paste("Scatterplot of ECDFs of", names[i], "and", names[j]))
  }
}
par(mfrow=c(1,1))

library(VineCopula)
for (i in 1:(ncol(datatrans)-1)) {
  for (j in (i+1):ncol(datatrans)) {
    summary(BiCopSelect(datatrans[,i], datatrans[,j], familyset=NA))
  }
}

library(copula)
rotateddata <- 1 - datatrans
fit <- fitCopula(gumbelCopula(dim=5), as.matrix(rotateddata), method="itau")
print(fit)

w <- c(0.20, 0.20, 0.20, 0.20, 0.20)

# Survival Gumbel and t-distributed marginals
n <- 1000000
Percentiles <- c()
for (i in 1:20) {
  sim <- rCopula(n, gumbelCopula(1.791, dim=5))
  sim <- 1 - sim
  
  sim2 <- as.matrix(qt(sim[,1], df=nu[1]) * sigma[1] + mu[1])
  for (i in 2:5) {
    sim2 <- cbind(sim2, qt(sim[,i], df=nu[i]) * sigma[i] + mu[i])
  }
  Percentiles <- c(Percentiles, sort(sim2 %*% w)[0.01 * n])
}

VaR <- mean(Percentiles)
print(VaR)

# Multivariate normal distribution
library(mvtnorm)
n <- 1000000
Percentiles <- c()
for (i in 1:20) {
  sim <- rmvnorm(n, mean=colMeans(data), sigma=var(data), method="chol")
  Percentiles <- c(Percentiles, sort(sim %*% w)[0.01 * n])
}

VaR <- mean(Percentiles)
print(VaR)

