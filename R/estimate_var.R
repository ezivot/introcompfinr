#' Estimate value at risk (VaR).
#'
#' This function estimates VaR assuming either a normal distribution,
#' log-normal distribution (i.e, normal distribution for continuously compounded returns),
#' or an empirical distribution. The VaR is estimated for a given probability level
#' and initial investment value. VaR is interpreted as the dollar loss of an initial investment
#' over a specified time period with a given probability level. For example, the 5% VaR of a $100,000
#' investment over a one-month period is interpreted as "there is a 5% chance that the portfolio will lose
#' VaR *or more* over the next month."
#'
#'
#' Normal VaR for simple returns is estimated using the formula the left-tail \eqn{\alpha}-quantile
#' of a normal distribution:
#' \deqn{\widehat{\mathrm{VaR}}_{\alpha} = -(\hat{\mu} + \hat{\sigma} \times z_{\alpha})\times w,}
#' where \eqn{\hat{\mu}} and \eqn{\hat{\sigma}} are the sample mean and sample standard
#' deviation of the returns, \eqn{z_\alpha} is the \eqn{\alpha}-quantile of the
#' standard normal distribution, and \eqn{w} is the initial investment value.
#' Log-normal VaR (normal VaR for cc returns) is estimated as
#' \deqn{\widehat{\mathrm{VaR}}_{\alpha} = -(e^{\hat{\mu} + \hat{\sigma} z_{\alpha}} - 1) \times w.}
#' Empirical VaR is estimated as the negative of the empirical quantile of the returns
#' multiplied by the initial investment value:
#' \deqn{\widehat{\mathrm{VaR}}_{\alpha} = -\hat{q}_{\alpha} \times w,}
#' where \eqn{\hat{q}_{\alpha}} is the empirical \eqn{\alpha}-quantile of the returns.
#'
#' @param x, a vector or matrix of returns. If a matrix, each column is
#' treated as a separate asset.
#'
#' @param alpha, scalar tail probability \eqn{\alpha}.
#'
#' @param w, scalar initial investment value.
#'
#' @param var_method, character string specifying the method to use for estimating VaR.
#' Options are "normal" or "empirical".
#'
#' @param var_type, character string specifying the type of VaR to estimate. Options ar
#' "dollar" or "return". If var_type is "dollar", then the VaR is estimated in dollar terms
#' by multiplying the estimated return VaR by the initial investment value. If var_type is
#' "return", then the VaR is estimated as the percentage loss in return terms.
#'
#' @param return_type, character string specifying the type of returns. Options are
#' "simple" or "cc" (continuously compounded). If return_type is "cc", then the
#' log-normal VaR formula is used for normal VaR estimation.
#'
#' @returns
#' A VaRhat object with the following components:
#'
#' * `call` The function call
#' * `date_range` The date range of the data
#' * `alpha` The probability level
#' * `w` The initial investment value
#' * `var_method` The method used for estimating VaR
#' * `var_type` The type of VaR estimated
#' * `return_type` The type of returns used for estimating VaR
#' * `VaRhat` The estimated value at risk \eqn{\widehat{\mathrm{VaR}}_{\alpha}}
#' * `se_VaRhat` The estimated VaR standard error \eqn{\widehat{\mathrm{se}}(\widehat{\mathrm{VaR}}_{\alpha})}
#'
#' @importFrom stats qnorm quantile

#' @export
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' exampleRetC = log(1 + exampleRetS)
#' var_hat = estimate_var(exampleRetS, alpha=0.05, w=100000,
#'                       var_method="normal", var_type ="return",
#'                       return_type="simple")
#' names(var_hat)

estimate_var = function(x, alpha=0.05, w=100000,
                        var_method=c("normal", "empirical"),
                        var_type=c("return", "dollar"),
                        return_type=c("simple", "cc")) {
  ## To-do
  ## 1. Add standard errors for empirical VaR estimates using bootstrap method or
  ## kernel density estimation. This will require additional arguments for the number
  ## of bootstrap samples and the bandwidth for kernel density estimation.
  ##
  call <- match.call()
  var_method=var_method[1]
  var_type = var_type[1]
  return_type=return_type[1]
  # extract date range if x is an xts object
  if(xts::is.xts(x)) {
    date_range = range(zoo::index(x))
  }

  x = as.matrix(x)
  muhat = colMeans(x)
  n_obs = nrow(x)
  sigmahat = matrixStats::colSds(x)
  se_muhat = sigmahat/sqrt(n_obs)
  # estimate VaR using normal or empirical distribution and compute estimated
  # standard errors using the delta method
  if (var_method == "normal") {
    qhat = muhat + sigmahat*qnorm(alpha)
    se_qhat = (sigmahat/sqrt(n_obs))*sqrt(1 + 0.5*qnorm(alpha)^2)
  } else {
    qhat = apply(x, 2, quantile, alpha)
    se_qhat = NULL
  }
  if (return_type == "simple") {
    VaRhat = -qhat
    se_VaRhat = se_qhat
  } else {
    VaRhat = -(exp(qhat) - 1)
    se_VaRhat = se_qhat*exp(qhat)
  }
  if (var_type == "dollar") {
    VaRhat = VaRhat*w
    se_VaRhat = w*se_qhat
  }
  ans = list("call"=call,
             "date_range"=date_range,
             "alpha"=alpha,
             "w"=w,
             "var_method"=var_method,
             "var_type"=var_type,
             "return_type"=return_type,
             "VaRhat"=VaRhat,
             "se_VaRhat"=se_VaRhat
             )
  class(ans) = "VaRhat"
  return(ans)
}


#' @title Print method of class VaRhat
#'
#' @description Print method for objects of class `VaRhat`. Show the estimated value at risk,
#' the probability level, the initial investment value, and the method used for
#' estimating VaR.
#'
#' @author Eric Zivot
#'
#' @param x object of class VaRhat
#'
#' @param ... additional arguments passed to \samp{print()}
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' exampleRetC = log(1 + exampleRetS)
#' var_hat = estimate_var(exampleRetS, alpha=0.05, w=100000,
#'                       var_method="normal", var_type ="return",
#'                       return_type="simple")
#' var_hat
#'
#' @export
print.VaRhat <-
  function(x, ...)
  {
    if(!is.null(x$date_range)) {
      cat("\nDate range: ", format(x$date_range[1], ...), " to ",
          format(x$date_range[2], ...), "\n")
    }
    cat("\nVaR method: ", format(x$var_method, ...), "\n")
    cat("VaR type: ", format(x$var_type, ...), "\n")
    cat("Return type: ", format(x$return_type, ...), "\n")
    cat("Probability level: ", format(x$alpha, ...), "\n")
    if(x$var_type == "dollar") {
    cat("Initial investment: ", format(x$w, ...), "\n")
    }
    cat("Estimated VaR:", "\n")
    print(round(x$VaRhat,4), ...)
    invisible(x)
  }


# Add summary method
#' @title Summary method of class VaRhat
#'
#' @description Show estimated VaR together with estimated standard errors
#' and 95% confidence intervals.
#'
#' @author Eric Zivot
#'
#' @param object, object of class VaRhat
#'
#' @param ... additional arguments passed to \samp{summary()}
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' exampleRetC = log(1 + exampleRetS)
#' var_hat = estimate_var(exampleRetS, alpha=0.05, w=100000,
#'                       var_method="normal", var_type ="return",
#'                       return_type="simple")
#' summary(var_hat)
#' @export
summary.VaRhat <-
  function(object, ...)
  {
    # compute approx 95% confidence interval for VaR estimate
    lower_var = object$VaRhat - 2*object$se_VaRhat
    upper_var = object$VaRhat + 2*object$se_VaRhat
    width = upper_var - lower_var
    ans = cbind(object$VaRhat, object$se_VaRhat, lower_var, upper_var, width)
    colnames(ans) = c("Estimates", "Std Errors", "2.5%", "97.5%", "Width")
    cat("\nVaR method: ", format(object$var_method, ...), "\n")
    cat("VaR type: ", format(object$var_type, ...), "\n")
    cat("Return type: ", format(object$return_type, ...), "\n")
    cat("Probability level: ", format(object$alpha, ...), "\n")
    if(object$var_type == "dollar") {
      cat("Initial investment: ", format(object$w, ...), "\n")
    }
    print(round(ans,4), ...)
    invisible(object)
  }
# Add plot method
# plot probability distribution together with VaR and confidence intervals

