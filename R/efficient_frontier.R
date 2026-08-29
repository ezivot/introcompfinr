#' @title Compute efficient frontier of risky assets
#'
#' @author Eric Zivot
#'
#' @description
#' The function constructs the set of mean-variance efficient portfolios that either allow all
#' assets to be sold short or not allow any asset to be sold short. The returned object is of class
#' \samp{Markowitz} for which there are \samp{print}, \samp{summary} and \samp{plot} methods.
#'
#' @details
#' If short sales are allowed (negative weights) then the set of efficient portfolios of risky
#' assets can be computed as a convex combination of any two efficient portfolios. It is convenient
#' to use the global minimum variance portfolio as one portfolio and an efficient portfolio with
#' target expected return equal to the maximum expected return of the assets under consideration as
#' the other portfolio. Call these portfolios \eqn{m} and \eqn{x}, respectively. Then for any number
#' \samp{alpha}, another efficient portfolio can be computed as \eqn{z=\alpha m+(1-\alpha)x}. If
#' short sales are not allowed, then the set of efficient portfolios is computed by repeated calls
#' to the function \samp{efficient_portfolio()}, with \samp{shorts=FALSE}, for a grid of target
#' expected returns starting at the expected return of the global minimum variance portfolio (not
#' allowing short sales) and ending at the expected return equal to the maximum expected return of
#' the assets under consideration.
#'
#' @param er_vec \samp{N x 1} vector of expected returns
#' @param cov_mat \samp{N x N} return covariance matrix
#' @param nport scalar, number of efficient portfolios to compute
#' @param alpha_min minimum value of \samp{alpha}, default is \samp{-.5}
#' @param alpha_max maximum value of \samp{alpha}, default is \samp{1.5}
#' @param shorts logical, if \samp{TRUE} then short sales (negative portfolio weights)
#' are allowed. If \samp{FALSE} then no asset is allowed to be sold short
#'
#' @return
#' An object of class \samp{Markowitz} with the following components:
#'
#' * `call` captures function call
#' * `er_vec` `nport x 1` vector of expected returns of efficient porfolios
#' * `sd_vec` `nport x 1` vector of standard deviations of efficient portfolios
#' * `weight_mat` `nport x N` matrix of weights of efficient portfolios
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
#' # tangency portfolio
#' tan_port <- tangency_portfolio(er, covmat, r_free)
#' # compute global minimum variance portfolio
#' gmin_port = globalmin_portfolio(er, covmat)
#'
#' # compute portfolio frontier
#' ef <- efficient_frontier(er, covmat, alpha_min=-2,
#'                          alpha_max=1.5, nport=20)
#' attributes(ef)
#'
#' plot(ef)
#' plot(ef, plot.assets=TRUE, col="blue", pch=16)
#' points(gmin_port$sd_port, gmin_port$er_port, col="green", pch=16, cex=2)
#' points(tan_port$sd_port, tan_port$er_port, col="red", pch=16, cex=2)
#' text(gmin_port$sd_port, gmin_port$er_port, labels="GLOBAL MIN", pos=2)
#' text(tan_port$sd_port, tan_port$er_port, labels="TANGENCY", pos=2)
#' sr_tan = (tan_port$er_port - r_free)/tan_port$sd_port
#' abline(a=r_free, b=sr_tan, col="green", lwd=2)
#'
#' @export efficient_frontier

efficient_frontier <-
  function(er_vec, cov_mat, nport=20, alpha_min=-0.5,
           alpha_max=1.5, shorts=TRUE)
  {
    call <- match.call()

    #
    # check for valid inputs
    #

    N <- length(er_vec)
    cov_mat <- as.matrix(cov_mat)
    if(N != nrow(cov_mat))
      stop("invalid inputs")
    if(any(diag(chol(cov_mat)) <= 0))
      stop("Covariance matrix not positive definite")
    if(is.null(names(er_vec)))
      names(er_vec) <- paste("Asset", 1:length(er_vec), sep="_")
    asset_names <- names(er_vec)
    #
    # create portfolio names
    #
    port_names <- rep("port",nport)
    ns <- seq(1,nport)
    port_names <- paste(port_names,ns)

    #
    # compute global minimum variance portfolio
    #
    cov_mat_inv <- solve(cov_mat)
    one.vec <- rep(1, N)
    port_gmin <- globalmin_portfolio(er_vec, cov_mat, shorts)
    w_gmin <- port_gmin$weight_vec

    if(shorts==TRUE){
      # compute efficient frontier as convex combinations of two efficient portfolios
      # 1st efficient port: global min var portfolio
      # 2nd efficient port: min var port with ER = max of ER for all assets
      er_max <- max(er_vec)
      port_max <- efficient_portfolio(er_vec,cov_mat,er_max)
      w_max <- port_max$weight_vec
      a <- seq(from=alpha_min,to=alpha_max,length=nport) # convex combinations
      we_mat <- (1-a) %o% w_gmin + a %o% w_max	         # rows are efficient portfolios
      er_e <- we_mat %*% er_vec					                 # expected returns of efficient portfolios
      er_e <- as.vector(er_e)
    } else if(shorts==FALSE){
      we_mat <- matrix(0, nrow=nport, ncol=N)
      we_mat[1,] <- w_gmin
      we_mat[nport, which.max(er_vec)] <- 1
      er_e <- as.vector(seq(from=port_gmin$er_port, to=max(er_vec), length=nport))
      for(i in 2:(nport-1))
        we_mat[i,] <- efficient_portfolio(er_vec, cov_mat, er_e[i], shorts)$weight_vec
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }

    names(er_e) <- port_names
    cov_e <- we_mat %*% cov_mat %*% t(we_mat) # cov mat of efficient portfolios
    sd_e <- sqrt(diag(cov_e))					        # std devs of efficient portfolios
    sd_e <- as.vector(sd_e)
    names(sd_e) <- port_names
    dimnames(we_mat) <- list(port_names,asset_names)

    #
    # summarize results
    #
    ans <- list("call" = call,
                "er_vec" = er_e,
                "sd_vec" = sd_e,
                "weight_mat" = we_mat)
    class(ans) <- "Markowitz"
    ans
  }

#' @title Print efficient frontier
#'
#' @author Eric Zivot
#'
#' @description
#' Print method for \samp{Markowitz} objects.
#'
#' @param x object of class Markowitz
#' @param ... additional arguments passed to \samp{print()}
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
#' # tangency portfolio
#' tan_port <- tangency_portfolio(er, covmat, r_free)
#' # compute global minimum variance portfolio
#' gmin_port = globalmin_portfolio(er, covmat)
#'
#' # compute portfolio frontier
#' ef <- efficient_frontier(er, covmat, alpha_min=-2,
#'                          alpha_max=1.5, nport=20)
#' attributes(ef)
#' print(ef)
#'
#' @export

print.Markowitz <-
  function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    xx <- rbind(x$er_vec,x$sd_vec)
    dimnames(xx)[[1]] <- c("ER","SD")
    cat("\nFrontier portfolios' expected returns and standard deviations\n")
    print(round(xx,4), ...)
    invisible(x)
  }

#' @title Summary method of class Markowitz
#'
#' @author Eric Zivot
#'
#' @description
#' Summary method for objects of class \samp{Markowitz}. For all portfolios on the efficient
#' frontier, the expected return, standard deviation and asset weights are shown. If
#' \samp{risk_free} is given then efficient portfolios that are combinations of the risk free asset
#' and the tangency portfolio are computed. The class \samp{summary_Markozitz} will be created.
#'
#' @param object object of class Markowitz
#' @param risk_free numeric, risk free rate
#' @param ... additional arguments passed to \samp{summary()}
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
#' # tangency portfolio
#' tan_port <- tangency_portfolio(er, covmat, r_free)
#' # compute global minimum variance portfolio
#' gmin_port = globalmin_portfolio(er, covmat)
#'
#' # compute portfolio frontier
#' ef <- efficient_frontier(er, covmat, alpha_min=-2,
#'                          alpha_max=1.5, nport=20)
#' attributes(ef)
#' summary(ef)
#'
#' @export

summary.Markowitz <-
  function(object, risk_free=NULL, ...)
  {
    call <- object$call
    asset_names <- colnames(object$weight_mat)
    port_names <- rownames(object$weight_mat)
    if(!is.null(risk_free)) {
      # compute efficient portfolios with a risk-free asset
      # this doesn't work
      nport <- length(object$er_vec)
      sd_max <- object$sd_vec[1]
      sd_e <- seq(from=0,to=sd_max,length=nport)
      names(sd_e) <- port_names

      #
      # get original er_vec and cov_mat data from call
      # these values are now part of the object but we want to get the original values from the call
      er_vec <- eval(object$call$er_vec)
      cov_mat <- eval(object$call$cov_mat)

      #
      # compute tangency portfolio
      tan_port <- tangency_portfolio(er_vec,cov_mat,risk_free)
      x_t <- sd_e/tan_port$sd_port		                        # weights in tangency port
      x_f <- 1 - x_t			                                # weights in t-bills
      er_e <- risk_free + x_t*(tan_port$er_vec - risk_free)
      names(er_e) <- port_names
      we_mat <- x_t %o% tan_port$weight_vec	              # rows are efficient portfolios
      dimnames(we_mat) <- list(port_names, asset_names)
      we_mat <- cbind(x_f,we_mat)
    }
    else {
      er_e <- object$er_vec
      sd_e <- object$sd_vec
      we_mat <- object$weight_mat
    }
    # Why create an object of class summary.Markowitz?
    # Because we want to have a different print method for the summary of the efficient frontier. The summary method will show the expected return, standard deviation and weights of each efficient portfolio. The print method for the summary.Markowitz class will show the expected return and standard deviation of each efficient portfolio.
    ans <- list("call" = call,
                "er"=er_e,
                "sd"=sd_e,
                "weights"=we_mat)
    class(ans) <- "summary.Markowitz"
    ans
  }

#' @title Plot method of class Markowitz
#'
#' @author Eric Zivot
#'
#' @description
#' Plot efficient frontier. The efficient frontier is a plot of portfolio expected return vs.
#' portfolio standard deviation for a collection of mean-variance efficient portfolios - portfolios
#' that minimize variance subject to a target expected return.
#'
#' @param x object of class Markowitz
#' @param plot.assets if \samp{TRUE} then plot asset \samp{sd} and \samp{er} with asset name labels
#' @param ... additional arguments passed to \samp{plot()}
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
#' # tangency portfolio
#' tan_port <- tangency_portfolio(er, covmat, r_free)
#' # compute global minimum variance portfolio
#' gmin_port = globalmin_portfolio(er, covmat)
#'
#' # compute portfolio frontier
#' ef <- efficient_frontier(er, covmat, alpha_min=-2,
#'                          alpha_max=1.5, nport=20)
#' attributes(ef)
#'
#' plot(ef)
#' plot(ef, plot.assets=TRUE, col="blue", pch=16)
#' points(gmin_port$sd_port, gmin_port$er_port, col="green", pch=16, cex=2)
#' points(tan_port$sd_port, tan_port$er_port, col="red", pch=16, cex=2)
#' text(gmin_port$sd_port, gmin_port$er_port, labels="GLOBAL MIN", pos=2)
#' text(tan_port$sd_port, tan_port$er_port, labels="TANGENCY", pos=2)
#' sr_tan = (tan_port$er_port - r_free)/tan_port$sd_port
#' abline(a=r_free, b=sr_tan, col="green", lwd=2)
#'
#' @export
#' @importFrom graphics text

plot.Markowitz <-
  function(x, plot.assets=FALSE, ...)
  {
    if (!plot.assets) {
      y.lim=c(0,max(x$er_vec))
      x.lim=c(0,max(x$sd_vec))
      plot(x$sd_vec,x$er_vec,type="b",xlim=x.lim, ylim=y.lim,
           xlab="Portfolio SD", ylab="Portfolio ER",
           main="Efficient Frontier", ...)
    }
    else {
      call = x$call
      mu.vals = eval(call$er_vec)
      sd.vals = sqrt( diag( eval(call$cov_mat) ) )
      y.lim = range(c(0,mu.vals,x$er_vec))
      x.lim = range(c(0,sd.vals,x$sd_vec))
      plot(x$sd_vec,x$er_vec,type="b", xlim=x.lim, ylim=y.lim,
           xlab="Portfolio SD", ylab="Portfolio ER",
           main="Efficient Frontier", ...)
      text(sd.vals, mu.vals, labels=names(mu.vals))
    }
    invisible()
  }
