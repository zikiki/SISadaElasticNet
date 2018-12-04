#' Sure Independence Screening
#' @param x Ultra high dimensional dependent variables matrix.
#' @param y Continuous response vector if \code{family} is \code{"gaussian"}.
#' @param pred number of relevant variable to select, and \code{pred} is smaller than n. A conservative \code{pred} is \code{(n-1)} and also can try \code{nlog(n)}.
#' @param scale If scale=TRUE, X will be scaled.
#' @return A smaller dimensional set of regressors for further modeling, \code{AdaElasticNet} model object,
#' @author Ziqi Zang
#'
#' @references
#' Jianqing Fan and Jinchi Lv (2008)
#' Sure Independence Screening for Ultrahigh Dimensional Feature Space (with discussion).
#' \emph{Journal of Royal Statistical Society} B, 70, 849-911.
#'
#' @export SIS.gaussian
#'
#' @examples
#'
#' x=matrix(rnorm(100*20),100,200)
#' y=rnorm(100)
#' Xsis=SIS.gaussian(x,y,pred = 20)

SIS.gaussian = function (X, Y, pred, scale = F)
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
  return(Xsis)
}

