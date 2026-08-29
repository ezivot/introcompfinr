#' @title Create portfolio object
#'
#' @author Eric Zivot
#'
#' @description
#' Create a portfolio object from expected return vector \eqn{\boldsymbol{\mu}},
#' covariance matrix \eqn{\boldsymbol{\Sigma}}, and weight vector \eqn{\mathbf{x}}.
#'
#' @details
#' To specify a portfolio, an expected return vector and covariance matrix for the assets under
#' consideration as well as a vector of portfolio weights are needed. The result of
#' `get_portfolio` is a `portfolio` object, which is list with components for the portfolio
#' expected return, portfolio standard deviation, and portfolio weights. There are `print`,
#' `summary` and `plot` methods.
#'
#' @param `er_vec` \eqn{N \times 1} vector of expected returns \eqn{\boldsymbol{\mu}}
#' @param `cov_mat` \eqn{N \times N} return covariance matrix \eqn{\boldsymbol{\Sigma}}
#' @param `weight_vec` \eqn{N \times 1} vector of portfolio weights \eqn{\mathbf{x}}
#'
#' @return An object of class `portfolio` with the following components:`
#'  * `call` captures function call
#'  * `er_port` portfolio expected return \eqn{\mu_p = \mathbf{x}'\boldsymbol{\mu}}
#'  * `sd_port` portfolio standard deviation \eqn{\sigma_p = \sqrt{\mathbf{x}'\boldsymbol{\Sigma}\mathbf{x}}}
#'  * `weight_vec` \eqn{N \times 1} vector of portfolio weights \eqn{\mathbf{x}}
#'  * `er_vec` \eqn{N \times 1} vector of expected returns \eqn{\boldsymbol{\mu}}
#'  * `cov_mat` \eqn{N \times N} return covariance matrix \eqn{\boldsymbol{\Sigma}}
#'
#' @examples
#' # construct the data
#' asset_names = c("MSFT", "NORD", "SBUX")
#' er = c(0.0427, 0.0015, 0.0285)
#' names(er) = asset_names
#' covmat = matrix(c(0.0100, 0.0018, 0.0011,
#'                   0.0018, 0.0109, 0.0026,
#'                   0.0011, 0.0026, 0.0199),
#'                 nrow=3, ncol=3)
#' r_free = 0.005
#' dimnames(covmat) = list(asset_names, asset_names)
#'
#' # compute equally weighted portfolio
#' ew = rep(1,3)/3
#' equal_weight_portfolio = get_portfolio(er_vec=er,cov_mat=covmat,weight_vec=ew)
#' class(equal_weight_portfolio)
#' names(equal_weight_portfolio)
#' equal_weight_portfolio
#' summary(equal_weight_portfolio)
#' plot(equal_weight_portfolio, col="blue")
#'
#' @export get_portfolio

get_portfolio <-
  function(er_vec, cov_mat, weight_vec)
  {
    call <- match.call()

    #
    # check for valid inputs
    #
    if (!is.vector(er_vec) || !is.vector(weight_vec))
      stop("er_vec and weight_vec must be vectors")
    if(length(er_vec) != length(weight_vec))
      stop("dimensions of er_vec and weight_vec do not match")
    cov_mat <- as.matrix(cov_mat)
    if(length(er_vec) != nrow(cov_mat))
      stop("dimensions of er_vec and cov_mat do not match")
    if(any(diag(chol(cov_mat)) <= 0))
      stop("Covariance matrix not positive definite")

    if(is.null(names(er_vec)))
      names(er_vec) <- paste("Asset", 1:length(er_vec), sep="_")
    asset_names <- names(er_vec)

    names(weight_vec) = asset_names
    #
    # compute portfolio expected return and volatility
    #
    er_port <- crossprod(er_vec,weight_vec)
    sd_port <- sqrt(weight_vec %*% cov_mat %*% weight_vec)
    #
    # create portfolio object
    #
    ans <- list("call" = call,
                "er_port" = as.numeric(er_port),
                "sd_port" = as.numeric(sd_port),
                "weight_vec" = weight_vec,
                "er_vec" = er_vec,
                "cov_mat" = cov_mat)
    class(ans) <- "portfolio"
    return(ans)
  }

#' @title Print method of class portfolio
#'
#' @author Eric Zivot
#'
#' @description
#' Print method for objects of class `portfolio`.
#'
#' @param `x` object of class portfolio
#' @param `...` additional arguments passed to `print()`
#'
#' @examples
#' # construct the data
#' asset_names = c("MSFT", "NORD", "SBUX")
#' er = c(0.0427, 0.0015, 0.0285)
#' names(er) = asset_names
#' covmat = matrix(c(0.0100, 0.0018, 0.0011,
#'                   0.0018, 0.0109, 0.0026,
#'                   0.0011, 0.0026, 0.0199),
#'                 nrow=3, ncol=3)
#' r.free = 0.005
#' dimnames(covmat) = list(asset_names, asset_names)
#'
#' # compute equally weighted portfolio
#' ew = rep(1,3)/3
#' equal_weight_portfolio = get_portfolio(er=er,cov_mat=covmat,weight_vec=ew)
#' print(equal_weight_portfolio)
#'
#' @export

print.portfolio <-
  function(x, ...)
  {
    cat("Call:\n")
    print(x$call, ...)
    cat("\nPortfolio expected return:    ", format(x$er_port, ...), "\n")
    cat("Portfolio standard deviation: ", format(x$sd_port, ...), "\n")
    cat("Portfolio weights:\n")
    print(round(x$weight_vec,4), ...)
    invisible(x)
  }

#' @title Summary method of class portfolio
#'
#' @author Eric Zivot
#'
#' @description
#' Summary method for objects of class `portfolio`. The output is the same as the `print`.
#' If `risk_free` is specified then the portfolio Sharpe ratio is also returned.
#'
#' @param `object` object of class portfolio
#' @param `risk_free` numeric, risk free rate
#' @param `...` additional arguments passed to `summary()`
#'
#' @examples
#' # construct the data
#' asset_names = c("MSFT", "NORD", "SBUX")
#' er = c(0.0427, 0.0015, 0.0285)
#' names(er) = asset_names
#' covmat = matrix(c(0.0100, 0.0018, 0.0011,
#'                   0.0018, 0.0109, 0.0026,
#'                   0.0011, 0.0026, 0.0199),
#'                 nrow=3, ncol=3)
#' r.free = 0.005
#' dimnames(covmat) = list(asset_names, asset_names)
#'
#' # compute equally weighted portfolio
#' ew = rep(1,3)/3
#' equal_weight_portfolio = get_portfolio(er=er,cov_mat=covmat,weight_vec=ew)
#' summary(equal_weight_portfolio)
#'
#' @export

summary.portfolio <-
  function(object, risk_free=NULL, ...)
  {
    cat("Call:\n")
    print(object$call)
    cat("\nPortfolio expected return:    ", format(object$er_port, ...), "\n")
    cat(  "Portfolio standard deviation: ", format(object$sd_port, ...), "\n")
    if(!is.null(risk_free)) {
      SharpeRatio <- (object$er_port - risk_free)/object$sd_port
      cat("Portfolio Sharpe Ratio:       ", format(SharpeRatio), "\n")
    }
    cat("Portfolio weights:\n")
    print(round(object$weight_vec,4), ...)
    invisible(object)
  }

#' @title Plot method of class portfolio
#'
#' @author Eric Zivot
#'
#' @description
#' The `plot()` method shows a bar chart of the portfolio weights.
#'
#' @param `x` object of class portfolio
#' @param `...` additional arguments passed to `barplot()`
#'
#' @examples
#' # construct the data
#' asset_names = c("MSFT", "NORD", "SBUX")
#' er = c(0.0427, 0.0015, 0.0285)
#' names(er) = asset_names
#' covmat = matrix(c(0.0100, 0.0018, 0.0011,
#'                   0.0018, 0.0109, 0.0026,
#'                   0.0011, 0.0026, 0.0199),
#'                 nrow=3, ncol=3)
#' r_free = 0.005
#' dimnames(covmat) = list(asset_names, asset_names)
#'
#' # compute equally weighted portfolio
#' ew = rep(1,3)/3
#' equal_weight_portfolio = get_portfolio(er=er,cov_mat=covmat,weight_vec=ew)
#' plot(equal_weight_portfolio, col="blue")
#'
#' @export
#' @importFrom graphics barplot

plot.portfolio <-
  function(x, ...)
  {
    asset_names <- names(x$weight_vec)
    barplot(x$weight_vec, names=asset_names,
            xlab="Assets", ylab="Weight", main="Portfolio Weights", ...)
    invisible()
  }
