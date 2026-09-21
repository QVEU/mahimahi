# Fits a 2-component normal mixture by a crude EM, then DELIBERATELY returns
# the components in an arbitrary order determined by the RNG -- mirroring the
# behaviour of real normalmixEM that the SCISSORS ordering fix guards against.
normalmixEM <- function(x, k = 2, ...) {
  q <- stats::quantile(x, c(0.25, 0.75))
  mu <- as.numeric(q); sigma <- rep(stats::sd(x) / 2, 2); lambda <- c(0.5, 0.5)
  for (iter in 1:50) {
    d <- sapply(1:2, function(j) lambda[j] * stats::dnorm(x, mu[j], sigma[j]))
    d[!is.finite(d)] <- 1e-300
    w <- d / pmax(rowSums(d), 1e-300)
    lambda <- pmax(colMeans(w), 1e-6)
    mu <- sapply(1:2, function(j) sum(w[, j] * x) / sum(w[, j]))
    sigma <- pmax(sapply(1:2, function(j) sqrt(sum(w[, j] * (x - mu[j])^2) / sum(w[, j]))), 1e-6)
  }
  # Arbitrary component order, as the real implementation gives.
  ord <- if (stats::runif(1) > 0.5) c(1, 2) else c(2, 1)
  structure(list(mu = mu[ord], sigma = sigma[ord], lambda = lambda[ord],
                 x = x, ft = "normalmixEM"), class = "mixEM")
}
