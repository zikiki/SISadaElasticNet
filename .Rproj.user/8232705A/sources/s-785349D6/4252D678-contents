#' Tune Hyperparameters for Adaptive Elastic Net
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"}.
#' @param family Model family, only \code{"gaussian"}.
#' @param init Type of the penalty used in the initial
#' estimation step. In the paper, they use \code{"enet"}.
#' @param alphas Vector of candidate \code{alpha}s to use in
#' \code{\link[glmnet]{cv.glmnet}}.
#' @param tune Parameter tuning method for each estimation step.
#' Possible options are \code{"cv"}, \code{"ebic"}, \code{"bic"},
#' and \code{"aic"}. Default is \code{"cv"}.
#' @param nfolds Fold numbers of cross-validation when \code{tune = "cv"}.
#' @param rule Lambda selection criterion when \code{tune = "cv"},
#' can be \code{"lambda.min"} or \code{"lambda.1se"}.
#' See \code{\link[glmnet]{cv.glmnet}} for details.
#' @param ebic.gamma Parameter for Extended BIC penalizing
#' size of the model space when \code{tune = "ebic"},
#' default is \code{1}. For details, see Chen and Chen (2008).
#' @param lower.limits Lower limits for coefficients.
#' Default is \code{-Inf}. For details, see \code{\link[glmnet]{glmnet}}.
#' @param upper.limits Upper limits for coefficients.
#' Default is \code{Inf}. For details, see \code{\link[glmnet]{glmnet}}.
#' @return List of model coefficients, \code{glmnet} model object,
#' and the optimal parameter set.
#'
#' @author Ziqi Zang
#'
#' @references
#' Zou, Hui, and Hao Helen Zhang. (2009).
#' On the adaptive elastic-net with a diverging number of parameters.
#' \emph{The Annals of Statistics} 37(4), 1733--1751.
#'
#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#'
#' @export AdaElasticNet
#'
#' @examples
#'
#' x=matrix(rnorm(100*20),100,20)
#' y=rnorm(100)
#' fit=AdaElasticNet(x,y,alphas = seq(0.2, 0.8, 0.2), seed = 10)
#' print(fit)
#' plot(fit)
tune.glmnet = function (x, y, family, alphas, tune, nfolds, rule, ebic.gamma,
          lower.limits, upper.limits, seed, parallel, ...)
{
  if (tune == "cv") {
    if (!parallel) {
      model.list = vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] = cv.glmnet(x = x, y = y, family = family,
                                    nfolds = nfolds, alpha = alphas[i], lower.limits = lower.limits,
                                    upper.limits = upper.limits, ...)
      }
    }
    else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.glmnet(x = x, y = y, family = family, nfolds = nfolds,
                  alpha = alphas, lower.limits = lower.limits,
                  upper.limits = upper.limits, ...)
      }
    }
    errors = unlist(lapply(model.list, function(x) min(sqrt(x$cvm))))
    errors.min.idx = which.min(errors)
    best.model = model.list[[errors.min.idx]]
    best.alpha = alphas[errors.min.idx]
    if (rule == "lambda.min")
      best.lambda = best.model$lambda.min
    if (rule == "lambda.1se")
      best.lambda = best.model$lambda.1se
    step.criterion = errors[errors.min.idx]
  }
  else {
    if (!parallel) {
      model.list = vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] = glmnet(x = x, y = y, family = family,
                                 alpha = alphas[i], lower.limits = lower.limits,
                                 upper.limits = upper.limits, ...)
      }
    }
    else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        glmnet(x = x, y = y, family = family, alpha = alphas,
               lower.limits = lower.limits, upper.limits = upper.limits,
               ...)
      }
    }
    if (tune == "aic") {
      ics.list = mapply(.aic, deviance = lapply(model.list,
                                                .deviance), df = lapply(model.list, .df), SIMPLIFY = FALSE)
    }
    if (tune == "bic") {
      ics.list = mapply(.bic, deviance = lapply(model.list,
                                                .deviance), df = lapply(model.list, .df), nobs = lapply(model.list,
                                                                                                        .nobs), SIMPLIFY = FALSE)
    }
    if (tune == "ebic") {
      ics.list = mapply(.ebic, deviance = lapply(model.list,
                                                 .deviance), df = lapply(model.list, .df), nobs = lapply(model.list,
                                                                                                         .nobs), nvar = lapply(model.list, .nvar), gamma = ebic.gamma,
                        SIMPLIFY = FALSE)
    }
    ics = sapply(ics.list, function(x) min(x))
    ics.min.idx = which.min(ics)
    best.model = model.list[[ics.min.idx]]
    best.alpha = alphas[ics.min.idx]
    best.ic.min.idx = which.min(ics.list[[ics.min.idx]])
    best.lambda = best.model$lambda[[best.ic.min.idx]]
    step.criterion = ics.list[[ics.min.idx]][[best.ic.min.idx]]
  }
  list(best.model = best.model, best.alpha = best.alpha, best.lambda = best.lambda,
       step.criterion = step.criterion)
}
