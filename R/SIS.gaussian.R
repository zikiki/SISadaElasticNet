function (X, Y, pred, scale = F)
{
  if (scale == T) {
    X <- scale(X)
  }
  p = dim(X)[2]
  IndicesSIS <- rep(0, p)
  beta <- rep(0, p)
  for (jj in 1:p) {
    beta[jj] <- abs(glm2(Y ~ X[, jj], family = gaussian)$coefficients[2])
  }
  IndicesSIS <- sort(beta, index = TRUE, decreasing = TRUE)$ix
  Xc <- X[, IndicesSIS[1:pred]]
  Xsis <- Xc
  new <- list(Xs = Xc, Index = IndicesSIS)
  return(new)
}
