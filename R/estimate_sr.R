#' Estimate Sharpe Ratio
#'
#' This function estimates the Sharpe ratio given a vector or matrix of returns.
#' The Sharpe ratio is estimated as the ratio of the mean excess return to the
#' standard deviation of excess returns. The function also computes the standard
#' error of the Sharpe ratio estimate using the delta method.
#'
#' The estimated Sharpe ratio of an asset is defined as
#' \deqn{\widehat{\mathrm{SR}} = \frac{\hat{\mu} - r_f}{\hat{\sigma}}}
#' where \eqn{\hat{\mu}} is the estimated mean return, \eqn{r_f} is the risk-free rate, and
#' \eqn{\hat{\sigma}} is the estimated standard deviation of the returns. The estimated
#' standard error of the Sharpe ratio is computed using the delta method assuming that
#' returns are described by the gaussian white noise model.
#'
#' @param x, a vector or matrix of returns. If a matrix, each column is
#' treated as a separate asset.
#'
#' @param r_f, scalar risk-free rate. The default is 0.
#'
#' @returns
#' A SRhat object with the following components:
#'
#' * `call` The function call
#' * `date_range` The date range of the data
#' * `srhat` The estimated Sharpe ratio \eqn{\widehat{\mathrm{SR}}}
#' * `se_srhat` The estimated standard error of the Sharpe ratio \eqn{\hat{\mathrm{se}}(\widehat{\mathrm{SR}})}
#' * `r_f` The risk-free rate \eqn{r_f}
#' * `n.obs` The number of observations \eqn{n}
#'
#' @export
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' sr_hat = estimate_sr(exampleRetS)
#' names(sr_hat)

estimate_sr = function(x, r_f=0) {
  ## To-do
  ## 1.
  ##
  call <- match.call()
  # extract date range if x is an xts object
  if(xts::is.xts(x)) {
    date_range = range(zoo::index(x))
  }

  x = as.matrix(x)
  muhat = colMeans(x)
  n_obs = nrow(x)
  sigmahat = matrixStats::colSds(x)
  se_muhat = sigmahat/sqrt(n_obs)
  # estimate SR and compute estimated standard errors using the delta method
  SRhat = (muhat - r_f)/sigmahat
  se_SRhat = sqrt((1+0.5*SRhat^2)/n_obs)

  ans = list("call"=call,
             "date_range"=date_range,
             "SRhat"=SRhat,
             "se_SRhat"=se_SRhat,
             "r_f"=r_f,
             "n.obs"=n_obs)
  class(ans) = "SRhat"
  return(ans)
}


#' @title Print method of class SRhat
#'
#' @description Print method for objects of class `SRhat`. Show the estimated
#' Sharpe ratio, the risk-free rate, and the date range.
#'
#' @author Eric Zivot
#'
#' @param x object of class SRhat
#'
#' @param ... additional arguments passed to `print()`
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' sr_hat = estimate_sr(exampleRetS)
#' sr_hat
#'
#' @export
print.SRhat <-
  function(x, ...)
  {
    if(!is.null(x$date_range)) {
      cat("\nDate range: ", format(x$date_range[1], ...), " to ",
          format(x$date_range[2], ...), "\n")
    }
    cat("Number of observations: ", x$n.obs, "\n")
    cat("Risk-free rate: ", format(x$r_f, ...), "\n")
    cat("Estimated Sharpe Ratio:", "\n")
    print(round(x$SRhat,4), ...)
    invisible(x)
  }


# Add summary method
#' @title Summary method of class SRhat
#'
#' @description Show estimated Sharpe ratio with estimated standard errors
#' and 95% confidence intervals.
#'
#' @author Eric Zivot
#'
#' @param object, object of class SRhat
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
#' sr_hat = estimate_sr(exampleRetS)
#' summary(sr_hat)
#' @export
summary.SRhat <-
  function(object, ...)
  {
    # compute approx 95% confidence interval for VaR estimate
    lower_sr = object$SRhat - 2*object$se_SRhat
    upper_sr = object$SRhat + 2*object$se_SRhat
    width = upper_sr - lower_sr
    ans = cbind(object$SRhat, object$se_SRhat, lower_sr, upper_sr, width)
    colnames(ans) = c("Estimates", "Std Errors", "2.5%", "97.5%", "Width")
    if(!is.null(object$date_range)) {
      cat("\nDate range: ", format(object$date_range[1], ...), " to ",
          format(object$date_range[2], ...), "\n")
    }
    cat("Number of observations: ", object$n.obs, "\n")
    cat("Risk-free rate: ", format(object$r_f, ...), "\n")
    print(round(ans,4), ...)
    invisible(object)
  }
# Add plot method

