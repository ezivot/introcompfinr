#' @title Compute tangency portfolio
#'
#' @author Eric Zivot
#'
#' @description
#' Compute tangency (maximum Sharpe ratio) portfolio. The portfolio can allow all assets to be
#' shorted or not allow any assets to be shorted.
#'
#' @details
#' The tangency portfolio \eqn{\mathbf{t}} is the portfolio of risky assets with the highest Sharpe's slope
#' and solves the optimization problem:
#' \deqn{ \max_{\mathbf{t}} \frac{\mathbf{t}'\boldsymbol{\mu}-r_f}{(\mathbf{t}'\boldsymbol{\Sigma} \mathbf{t})^{1/2}}}
#' s.t. \eqn{\mathbf{t}'\mathbf{1}=1,} where \eqn{r_f} denotes the risk-free rate.
#' If short sales are allowed then there is an analytic solution using matrix algebra of the form
#' \deqn{\mathbf{t} = \frac{\boldsymbol{\Sigma}^{-1}(\boldsymbol{\mu}-r_f)}{\mathbf{1}'\boldsymbol{\Sigma}^{-1}(\boldsymbol{\mu}-r_f))}.}
#' If short sales are not allowed then the maximum Sharpe ratio portfolio must be computed numerically
#' using quadratic programming.
#'
#' @param er_vec \eqn{N \times 1} vector of expected returns \eqn{\boldsymbol{\mu}}
#' @param cov_mat \eqn{N \times N} return covariance matrix \eqn{\boldsymbol{\Sigma}}
#' @param risk_free numeric, risk free rate \eqn{r_f} used to compute excess returns
#' @param shorts logical, if `TRUE` then short sales (negative portfolio weights)
#' are allowed. If `FALSE` then no asset is allowed to be sold short.
#'
#' @return
#' A `portfolio` object with the following components:
#'  * `call` captures function call
#'  * `er_port` portfolio expected return \eqn{\mu_t \ \mathbf{t}'\boldsymbol{\mu}}
#'  * `sd_port` portfolio standard deviation \eqn{\sigma_t = \sqrt{\mathbf{t}'\boldsymbol{\Sigma} \mathbf{t}}}
#'  * `weight_vec` \eqn{N \times 1} vector of tangency portfolio weights \eqn{\mathbf{t}}
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
#' # compute tangency portfolio
#' tan_port <- tangency_portfolio(er, covmat, r_free)
#' tan_port
#' summary(tan_port, risk_free=r_free)
#' plot(tan_port, col="blue")
#'
#' # compute tangency portfolio with no short sales
#' tan_port_ns <- tangency_portfolio(er, covmat, r_free, shorts=FALSE)
#' tan_port_ns
#' summary(tan_port_ns, risk_free=r_free)
#' plot(tan_port_ns, col="blue")
#'
#' @export tangency_portfolio

tangency_portfolio <-
  function(er_vec,cov_mat,risk_free, shorts=TRUE)
  {
    call <- match.call()

    #
    # check for valid inputs
    #
    if(risk_free < 0)
      stop("Risk-free rate must be positive")
    er <- as.vector(er)
    cov_mat <- as.matrix(cov_mat)
    N <- length(er)
    if(N != nrow(cov_mat))
      stop("invalid inputs")
    if(any(diag(chol(cov_mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov_mat is positive semi-definite
    if(is.null(names(er_vec)))
      names(er_vec) <- paste("Asset", 1:length(er_vec), sep="_")
    asset_names <- names(er_vec)

    #
    # compute global minimum variance portfolio
    #
    gmin_port <- globalmin_portfolio(er_vec, cov_mat, shorts=shorts)
    if(gmin_port$er_port < risk_free)
      stop("Risk-free rate greater than avg return on global minimum variance portfolio")

    #
    # compute tangency portfolio
    #
    if(shorts==TRUE){
      cov_mat_inv <- solve(cov_mat)
      w_t <- cov_mat_inv %*% (er - risk_free) # tangency portfolio
      w_t <- as.vector(w_t/sum(w_t))          # normalize weights
    } else if(shorts==FALSE){
      Dmat <- 2*cov_mat
      dvec <- rep.int(0, N)
      er_excess <- er_vec - risk_free
      Amat <- cbind(er_excess, diag(1,N))
      bvec <- c(1, rep(0,N))
      result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w_t <- round(result$solution/sum(result$solution), 6)
    } else {
      stop("Shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }

    names(w_t) <- asset_names
    tan_port <- get_portfolio(er_vec=er_vec, cov_mat=cov_mat,
                               weight_vec=w_t)
    return(tan_port)
  }
