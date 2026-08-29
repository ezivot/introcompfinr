#' @title Compute global minimum variance portfolio
#'
#' @author Eric Zivot
#'
#' @description
#' Compute global minimum variance portfolio given expected return vector \eqn{\boldsymbol{\mu}}
#' and covariance matrix \eqn{\boldsymbol{\Sigma}}. The resulting
#' portfolio \eqn{\mathbf{m}} can allow all assets to be shorted or not allow any assets to be shorted.
#' The returned  object is of class `portfolio`.
#'
#' @details
#' The global minimum variance portfolio \eqn{\mathbf{m}} allowing for short sales solves the optimization
#' problem:
#' \deqn{\min_{\mathbf{m}} \mathbf{m}'\boldsymbol{\Sigma} \mathbf{m} ~ s.t.~ \mathbf{m}'\mathbf{1}=1,}
#' for which there is an analytic solution using matrix algebra:
#' \deqn{\mathbf{m} = \frac{\boldsymbol{\Sigma}^{-1}\mathbf{1}}{\mathbf{1}'\boldsymbol{\Sigma}^{-1}\mathbf{1}}.}
#' If short sales are not allowed then the portfolio is computed numerically using
#' the function `solve.QP()` from the **quadprog** package.
#'
#' @param `er_vec` \eqn{N \times 1} vector of expected returns \eqn{\boldsymbol{\mu}}
#' @param `cov_mat` \eqn{N \times N} return covariance matrix \eqn{\boldsymbol{\Sigma}}
#' @param `shorts` logical, if `TRUE` then short sales (negative portfolio weights)
#' are allowed. If `FALSE` then no asset is allowed to be sold short.
#'
#' @return
#' A `portfolio` object with the following components:
#'  * `call` captures function call
#'  * `er_port` portfolio expected return \eqn{\mu_m = \mathbf{m}'\boldsymbol{\mu}}`
#'  * `sd_port` portfolio standard deviation \eqn{\sigma_m = \sqrt{\mathbf{m}'\boldsymbol{\Sigma}\mathbf{m}}}
#'  * `weight_vec` \eqn{N \times 1} vector of portfolio weights \eqn{\mathbf{m}}
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
#' # compute global minimum variance portfolio
#' gmin_port = globalmin_portfolio(er, covmat)
#' attributes(gmin_port)
#' print(gmin_port)
#' summary(gmin_port, risk_free=r_free)
#' plot(gmin_port, col="blue")
#'
#' # compute global minimum variance portfolio with no short sales
#' gmin_port_ns = globalmin_portfolio(er, covmat, shorts=FALSE)
#' attributes(gmin_port_ns)
#' print(gmin_port_ns)
#' summary(gmin_port_ns, risk_free=r_free)
#' plot(gmin_port_ns, col="blue")
#'
#' @export globalmin_portfolio

globalmin_portfolio <-
  function(er_vec, cov_mat, shorts=TRUE)
  {
    call <- match.call()

    #
    # check for valid inputs
    #
    cov_mat <- as.matrix(cov_mat)
    N <- length(er_vec)
    if(N != nrow(cov_mat))
      stop("invalid inputs")
    if(any(diag(chol(cov_mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov_mat is positive semi-definite
    if(is.null(names(er_vec)))
      names(er_vec) <- paste("Asset", 1:length(er_vec), sep="_")
    asset_names <- names(er_vec)
    #
    # compute global minimum portfolio
    #
    if(shorts==TRUE){
      cov_mat_inv <- solve(cov_mat)
      one_vec <- rep(1,N)
      w_gmin <- rowSums(cov_mat_inv) / sum(cov_mat_inv)
      w_gmin <- as.vector(w_gmin)
    } else if(shorts==FALSE){
      Dmat <- 2*cov_mat
      dvec <- rep.int(0, N)
      Amat <- cbind(rep(1,N), diag(1,N))
      bvec <- c(1, rep(0,N))
      result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w_gmin <- round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }

    names(w_gmin) <- asset_names
    er_gmin <- crossprod(w_gmin,er_vec)
    sd_gmin <- sqrt(t(w_gmin) %*% cov_mat %*% w_gmin)
    gmin_port <- get_portfolio(er_vec=er_vec, cov_mat=cov_mat,
                               weight_vec=w_gmin)
    return(gmin_port)
  }
