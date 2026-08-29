#' Estimate Gaussian White Noise Model
#'
#' Estimate Gaussian White Noise (GWN) model parameters from a vector or matrix of data.
#' The function estimates the mean vector, covariance matrix, and correlation matrix of the data,
#' as well as their standard errors. The function returns an object of class `gwn_model` which is a
#' list containing the estimated parameters and their standard errors, along with the number of observations.
#' There are corresponding `print`, `summary`, and `plot` methods for objects of class `gwn_model`.
#'
#' The GWN model assumes that the returns are independent and identically distributed (i.i.d.) with a
#' multivariate normal distribution:
#' \deqn{\mathbf{R}_t \sim N(\boldsymbol{\mu}, \boldsymbol{\Sigma}),}
#' where \eqn{\mathbf{R}_t} is the \eqn{n \times 1} vector of returns at time \eqn{t}, \eqn{\boldsymbol{\mu}} is the
#' \eqn{n \times 1} mean vector, and \eqn{\boldsymbol{\Sigma}} is the \eqn{n \times n} covariance matrix.
#' The estimators of \eqn{\boldsymbol{\mu}} and \eqn{\boldsymbol{\Sigma}} are the sample mean and sample
#' covariance matrix, respectively.
#'
#' @param x `n.obs x n.vars` xts time series, matrix or vector of data. If a vector is provided,
#' it will be converted to a matrix with one column.`
#'
#' @returns An object of class `gwn_model` with the following components:
#'
#' * `call` The matched function call
#' * `muhat` Estimated mean vector \eqn{\hat{\boldsymbol{\mu}}}
#' * `se_muhat` Estimated standard errors of the mean vector \eqn{\hat{\mathrm{se}}(\hat{\mu}_i),~ i=1,\ldots,n}
#' * `Sigmahat` Estimated covariance matrix \eqn{\hat{\boldsymbol{\Sigma}}}
#' * `Rhat` Estimated correlation matrix \eqn{\hat{\mathbf{R}}}
#' * `sigma2hat` Estimated variances \eqn{\hat{\sigma}^2_i,~ i=1,\ldots,n}
#' * `se_sigma2hat` Estimated standard errors of the sample variances \eqn{\hat{\mathrm{se}}(\hat{\sigma}^2_i),~ i=1,\ldots,n}
#' * `sigmahat` Estimated standard deviations \eqn{\hat{\sigma}_i,~ i=1,\ldots,n}
#' * `se_sigmahat` Estimated standard errors of the sample standard deviations \eqn{\hat{\mathrm{se}}(\hat{\sigma}_i),~ i=1,\ldots,n}
#' * `covhat` Estimated covariances \eqn{\hat{\sigma}_{ij},~ i,j=1,\ldots,n}
#' * `se_cov` Estimated standard errors of the sample covariances \eqn{\hat{\mathrm{se}}(\hat{\sigma}_{ij}),~ i,j=1,\ldots,n}
#' * `rhohat` Estimated correlations \eqn{\hat{\rho}_{ij},~ i,j=1,\ldots,n}
#' * `se_rhohat` Estimated standard errors of the sample correlations \eqn{\hat{\mathrm{se}}(\hat{\rho}_{ij}),~ i,j=1,\ldots,n}
#' * `n.obs` Number of observations \eqn{n}
#' * `ehat` Estimated residuals \eqn{\hat{\mathbf{e}}_t = \mathbf{R}_t - \hat{\boldsymbol{\mu}}}
#'
#' @references See Zivot, E. (2026) *Introduction to Computational Finance and Financial Econometrics with R*,
#' chapter 7 for a discussion of the Gaussian White Noise model and its estimation.
#'
#' @importFrom stats cov cor setNames var
#'
#' @importFrom utils combn
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
#' gwn_fit = estimate_gwn_model(exampleRetS[,"AMZN"])
#' gwn_fit2 = estimate_gwn_model(exampleRetS)
estimate_gwn_model <- function(x) {
  ## To-do:
  ## 1. Correct formula for standard errors for covariances
  ## 2. Add vcov method to return covariance matrix of estimated parameters
  ##    vcov for means, vcov for volatilities, vcov for both means and volatilities.
  ##    utilize block diagonal structure of vcov matrix for means and volatilities.
  ##
  call <- match.call()
  if(xts::is.xts(x)) {
    date_range = range(zoo::index(x))
    x_mat = zoo::coredata(x)
  }
  asset_names = colnames(x)
  if(ncol(x) > 1) {
    x_mat = as.matrix(x)
  } else {
    x_mat = matrix(x, ncol=1)
  }
  if(ncol(x) > 1) {
    n.obs = nrow(x_mat)
    # estimate mean vector and covariance matrix
    muhat <- colMeans(x_mat)
    Sigmahat <- cov(x_mat)
    Rhat = cor(x_mat)
    sigma2hat = diag(Sigmahat)
    sigmahat = sqrt(diag(Sigmahat))
    # extract pairwise covariances and correlations from Sigmahat and Rhat
    pair_names = combn(asset_names, 2, paste, collapse = ".")
    covhat = setNames(Sigmahat[lower.tri(Sigmahat)], pair_names)
    rhohat = setNames(Rhat[lower.tri(Rhat)], pair_names)

    # estimate standard errors
    se_muhat = sigmahat/sqrt(n.obs)
    se_sigma2hat = sigma2hat*sqrt(2/n.obs)
    se_sigmahat = se_sigma2hat/(2*sigmahat)
    pair_idx = combn(length(asset_names), 2)
    sigma2_products = sigma2hat[pair_idx[1, ]] * sigma2hat[pair_idx[2, ]]
    se_cov = setNames(sqrt((1 + rhohat^2) * sigma2_products / n.obs), pair_names)
    se_rhohat = setNames(sqrt((1 - rhohat^2)^2/n.obs), pair_names)

    # extract residuals
    ehat = sweep(x, 2, muhat)
  } else {
    # univariate case
    n.obs = nrow(x_mat)
    # estimate mean, variance and standard deviation
    muhat <- mean(x_mat)
    sigma2hat = as.numeric(var(x_mat))
    sigmahat = sqrt(sigma2hat)
    Sigmahat <- NULL
    Rhat = NULL
    covhat = NULL
    rhohat = NULL

    # estimate standard errors
    se_muhat = sigmahat/sqrt(n.obs)
    se_sigma2hat = sigma2hat*sqrt(2/n.obs)
    se_sigmahat = se_sigma2hat/(2*sigmahat)
    se_cov = NULL
    se_rhohat = NULL

    # extract residuals
    ehat = x - muhat
    colnames(ehat) = asset_names
  }
  # create gwn_model object

  ans = list("call"=call,
             "date_range"=date_range,
             "muhat"=muhat,
             "se_muhat"=se_muhat,
             "Sigmahat"=Sigmahat,
             "Rhat"=Rhat,
             "sigma2hat"=sigma2hat,
             "se_sigma2hat"=se_sigma2hat,
             "sigmahat"=sigmahat,
             "se_sigmahat"=se_sigmahat,
             "covhat"=covhat,
             "se_cov"=se_cov,
             "rhohat"=rhohat,
             "se_rhohat"=se_rhohat,
             "n.obs"=n.obs,
             "ehat"=ehat
             )
  class(ans) = "gwn_model"
  return(ans)

}

#' Print method for objects of class `gwn_model`
#'
#' @param x object of class `gwn_model`
#' @param ... ... additional arguments passed to \samp{print()}
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
#' gwn_fit = estimate_gwn_model(exampleRetS[,"AMZN"])
#' gwn_fit2 = estimate_gwn_model(exampleRetS)
#' gwn_fit
#' gwn_fit2
print.gwn_model <- function(x, ...) {
  cat("Call:\n")
  print(x$call, ...)
  if(!is.null(x$date_range)) {
    cat("\nDate range: ", format(x$date_range[1], ...), " to ",
        format(x$date_range[2], ...), "\n")
  }
  cat("Number of observations: ", x$n.obs, "\n")
  cat("\n")
  coefs = cbind(x$muhat, x$sigmahat)
  colnames(coefs) = c("Mean", "Volatility")
  rownames(coefs) = colnames(x$ehat)
  print(round(coefs, 4), ...)
  invisible(x)
}

#' @title Summary method of class gwn_model
#'
#' @description Show estimated mean and volatilies together with estimated standard errors
#' and 95% confidence intervals.
#'
#' @author Eric Zivot
#'
#' @param object, object of class gwn_fit
#'
#' @param show_correlations logical, whether to show pairwise correlations
#'
#' @param ... additional arguments passed to \samp{summary()}
#'
#' @references See Zivot, E. (2026) *Introduction to Computational Finance and Financial Econometrics with R*,
#' chapter 7 for a discussion of the Gaussian White Noise model and its estimation.
#'
#' @examples
#' examplePrices = xts::merge.xts(amznDailyPrices, baDailyPrices, costDailyPrices)
#' examplePrices = xts::to.monthly(examplePrices,OHLC=FALSE)
#' examplePrices=examplePrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' gwn_fit = estimate_gwn_model(exampleRetS[,"AMZN"])
#' gwn_fit2 = estimate_gwn_model(exampleRetS)
#' summary(gwn_fit)
#' summary(gwn_fit2)
#'
#' @export
summary.gwn_model <-
  function(object, show_correlations=FALSE, ...)
  {
    # To-do: add summary for covariances and correlations
    # compute approx 95% confidence interval for mean and volatility estimates
    lower_muhat = object$muhat - 2*object$se_muhat
    upper_muhat = object$muhat + 2*object$se_muhat
    width_muhat = upper_muhat - lower_muhat

    lower_sigmahat = object$sigmahat - 2*object$se_sigmahat
    upper_sigmahat = object$sigmahat + 2*object$se_sigmahat
    width_sigmahat = upper_sigmahat - lower_sigmahat

    ans = cbind(object$muhat, object$se_muhat, lower_muhat, upper_muhat, width_muhat,
                 object$sigmahat, object$se_sigmahat, lower_sigmahat, upper_sigmahat, width_sigmahat)
    colnames(ans) = c("Mean", "Std Error", "2.5%", "97.5%", "Width",
                      "Volatility", "Std Error", "2.5%", "97.5%", "Width")
    rownames(ans) = colnames(object$ehat)

    cat("Call:\n")
    print(object$call,...)
    if(!is.null(object$date_range)) {
      cat("\nDate range: ", format(object$date_range[1], ...), " to ",
          format(object$date_range[2], ...), "\n")
    }
    cat("Number of observations: ", object$n.obs, "\n")
    cat("\n")
    print(round(ans,4), ...)

    if(show_correlations) {
      ans_corr = cbind(object$rhohat, object$se_rhohat)
      lower_rhohat  = object$rhohat - 2*object$se_rhohat
      upper_rhohat  = object$rhohat + 2*object$se_rhohat
      width_rhohat  = upper_rhohat - lower_rhohat
      ans_corr = cbind(ans_corr, lower_rhohat, upper_rhohat,
                       width_rhohat)
      colnames(ans_corr) = c("Correlation", "Std Error", "2.5%",
                             "97.5%", "Width")
      cat("\nEstimated pairwise correlations:\n")
      print(round(ans_corr,4), ...)

    }
    invisible(object)
  }

#' Plot method of class gwn_model
#'
#' @param x object of class gwn_model
#'
#' @param asset character, name of asset to plot. If NULL, the first asset is used.
#'
#' @param which.plot character, one of "residuals", "fitted", "histogram", "acf",
#' "qq", "boxplot", "correlation", or "risk-return".
#'
#' @param ... additional arguments passed to the plotting functions
#'
#' @returns A plot of the specified type for the specified asset or for all assets in the model.
#'
#' @importFrom stats acf qqnorm qqline
#'
#' @importFrom graphics boxplot hist
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
#' gwn_fit = estimate_gwn_model(exampleRetS[,"AMZN"])
#' gwn_fit2 = estimate_gwn_model(exampleRetS)
plot.gwn_model <- function(x, asset = NULL, which.plot=c("residuals", "fitted", "histogram",
                                           "acf", "qq", "boxplot", "correlation",
                                           "risk-return"), ...) {
  which.plot = which.plot[1]
  asset.names = colnames(x$ehat)


  if(!(which.plot %in% c("residuals", "fitted", "histogram", "acf", "qq",
                         "boxplot", "correlation", "risk-return"))) {
    stop("which.plot must be one of 'residuals', 'fitted', 'histogram', 'acf',
         'qq', 'boxplot', 'correlation' or 'risk-return'")
  }

  if(which.plot == "correlation") {
    if(ncol(x$ehat) < 2) {
      stop("Correlation plot requires at least two assets")
    }
    corrplot::corrplot(x$Rhat, method="color", type="upper",
                       tl.col="black", tl.srt=45, ...)
  } else if(which.plot == "risk-return") {
    if(ncol(x$ehat) < 2) {
      stop("Risk-return plot requires at least two assets")
    }
    plot(x$sigmahat, x$muhat, xlab="Volatility", ylab="Mean Return",
         main="Risk-Return Plot", pch=19, col="blue", ...)
    text(x$sigmahat, x$muhat, labels=asset.names, pos=4, cex=0.8)
  } else {
    if(is.null(asset)) {
      asset = asset.names[1]
    }
    if(!(asset %in% asset.names)) {
      stop(paste("Asset", asset, "not found in model"))
    }
    if(ncol(x$ehat) == 1) {
      ehat = x$ehat
      fitted = xts::xts(rep(x$muhat, nrow(ehat)), order.by=zoo::index(ehat))
      returns = x$ehat + fitted
    } else {
      ehat = x$ehat[, asset]
      fitted = xts::xts(rep(x$muhat[asset], nrow(ehat)), order.by=zoo::index(ehat))
      returns = ehat + fitted
    }
    if(which.plot == "residuals") {
      plot(ehat, main=paste("Residuals for", asset), ylab="Residuals", ...)
    } else if(which.plot == "fitted") {
      data_to_plot = xts::merge.xts(fitted, returns)
      plot(data_to_plot, main=paste("Fitted Values for", asset), ylab="Fitted Values", ...)
    } else if(which.plot == "histogram") {
      hist(ehat, main=paste("Histogram of Residuals for", asset), xlab="Residuals", ...)
    } else if(which.plot == "acf") {
      acf(zoo::coredata(ehat), main=paste("ACF of Residuals for", asset), ...)
    } else if(which.plot == "qq") {
      qqnorm(ehat, main=paste("QQ Plot of Residuals for", asset), ...)
      qqline(ehat)
    } else if(which.plot == "boxplot") {
      boxplot(ehat, main=paste("Boxplot of Residuals for", asset), ylab="Residuals", ...)
    }
  }

}


#' Extract Model Coefficients from GWN Model
#'
#' @param object of class gwn_model
#' @param which_coef character, one of "mean", "volatility", "variance", "covariance",
#' or "correlation"
#'
#' @param ... additional arguments passed to the function
#'
#' @description Extract the specified coefficients from a GWN model.
#'
#' @author Eric Zivot
#'
#' @returns Coefficients extracted from the model `object`
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
#' gwn_fit = estimate_gwn_model(exampleRetS[,"AMZN"])
#' gwn_fit2 = estimate_gwn_model(exampleRetS)
#' coef(gwn_fit, which_coef="mean")

coef.gwn_model <- function(object, which_coef = c("mean", "volatility", "variance",
                                                   "covariance", "correlation"), ...) {
  which_coef = which_coef[1]
  if(!(which_coef %in% c("mean", "volatility", "variance", "covariance", "correlation"))) {
    stop("which_coef must be one of 'mean', 'volatility', 'variance', 'covariance', or 'correlation'")
  }
  if(which_coef == "mean") {
    return(object$muhat)
  } else if(which_coef == "volatility") {
    return(object$sigmahat)
  } else if(which_coef == "variance") {
    return(object$sigma2hat)
  } else if(which_coef == "covariance") {
    return(object$covhat)
  } else if(which_coef == "correlation") {
    return(object$rhohat)
  }
}
