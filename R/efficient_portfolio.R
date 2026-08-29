#' @title Compute minimum variance portfolio subject to target return
#'
#' @author Eric Zivot
#'
#' @description
#' Compute minimum variance portfolio subject to target return either allowing all assets to be sold
#' short or not allowing any asset to be sold short. The returned object is of class
#' `portfolio`.
#'
#' @details
#' A mean-variance efficient portfolio \eqn{\mathbf{x}} allowing short sales (negative weights) that achieves
#' a target expected return \eqn{\mu_0} solves the optimization problem:
#' \deqn{\min_{\mathbf{x}} \mathbf{x}'\boldsymbol{\Sigma} \mathbf{x} ~ s.t.~
#' \mathbf{x}'\mathbf{1}=1 ~ \mathrm{and} ~  \mathbf{x}'\boldsymbol{\mu}=\mu_0,}
#' for which there is an analytic solution using matrix
#' algebra. If short sales are not allowed then the portfolio is computed numerically using the
#' function `solve.QP()` from the `quadprog` package.
#'
#' @param `er_vec` \eqn{N \times 1} vector of expected returns
#' @param `cov_mat` \eqn{N \times N} return covariance matrix
#' @param `target_return` scalar, target expected return
#' @param `shorts` logical, if `TRUE` then short sales (negative portfolio weights)
#' are allowed. If `FALSE` then no asset is allowed to be sold short.
#'
#' @return
#' A `portfolio` object with the following components:
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
#' # compute efficient portfolio subject to target return
#' target_return = er["MSFT"]
#' e_port_msft = efficient_portfolio(er, covmat, target_return)
#' e_port_msft
#' summary(e_port_msft, risk_free=r_free)
#' plot(e_port_msft, col="blue")
#'
#' # compute efficient portfolio subject to target return with no short sales
#' target_return = er["MSFT"]
#' e_port_msft_ns = efficient_portfolio(er, covmat, target_return, shorts=FALSE)
#' e_port_msft_ns
#' summary(e_port_msft_ns, risk_free=r_free)
#' plot(e_port_msft_ns, col="blue")
#'
#' @export efficient_portfolio

efficient_portfolio <-
  function(er_vec, cov_mat, target_return, shorts=TRUE)
  {
    call <- match.call()

    #
    # check for valid inputs
    #

    N <- length(er_vec)
    cov_mat <- as.matrix(cov_mat)
    if(N != nrow(cov_mat))
      stop("er_vec and cov_mat are not comformable.")
    if(any(diag(chol(cov_mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov_mat is positive semi-definite
    if(is.null(names(er_vec)))
      names(er_vec) <- paste("Asset", 1:length(er_vec), sep="_")
    asset_names <- names(er_vec)

    #
    # compute efficient portfolio
    #
    if(shorts==TRUE){
      ones <- rep(1, N)
      top <- cbind(2*cov_mat, er_vec, ones)
      bot <- cbind(rbind(er_vec, ones), matrix(0,2,2))
      A <- rbind(top, bot)
      b_target <- as.matrix(c(rep(0, N), target_return, 1))
      x <- solve(A, b_target)
      w <- x[1:N]
    } else if(shorts==FALSE){
      Dmat <- 2*cov_mat
      dvec <- rep.int(0, N)
      Amat <- cbind(rep(1,N), er_vec, diag(1,N))
      bvec <- c(1, target_return, rep(0,N))
      result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
      w <- round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }

    #
    # compute portfolio expected returns and variance
    #
    names(w) <- asset_names
    ans <- get_portfolio(er_vec=er_vec, cov_mat = cov_mat,
                         weight_vec=w)
    return(ans)
  }
