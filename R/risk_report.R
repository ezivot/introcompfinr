#' Compute a Volatility Risk Report
#'
#' Generic function to compute an asset-level risk report based on portfolio
#' volatility. A method is defined for objects of class `portfolio`.
#'
#' @param `object` An object for which a risk report is defined.
#'
#' @param `...` Additional arguments passed to methods.
#'
#' @returns A numeric matrix with one row per asset plus a `"Portfolio"` summary
#'   row. Returns invisibly.
#'
#' @export
risk_report <- function(object, ...) {
  UseMethod("risk_report")
}

#' Volatility Risk Report for a Portfolio
#'
#' Compute asset-level risk contributions based on portfolio volatility for an
#' object of class `portfolio`. The report includes dollar allocations, weights,
#' individual asset volatilities, marginal contributions to risk (MCR),
#' contributions to risk (CR), percentage contributions to risk (PCR),
#' asset-portfolio correlations, and portfolio betas.
#'
#' @details
#' The risk quantities are defined as follows, where \eqn{\mathbf{x}} is the
#' weight vector, \eqn{\boldsymbol{\Sigma}} is the covariance matrix, and
#' \eqn{\sigma_p = \sqrt{\mathbf{x}'\boldsymbol{\Sigma}\mathbf{x}}} is portfolio volatility:
#'
#' \describe{
#'   \item{MCR}{\eqn{(\boldsymbol{\Sigma}\mathbf{x}) / \sigma_p} — marginal contribution to risk}
#'   \item{CR}{\eqn{x_i \cdot \text{MCR}_i} — contribution to risk; sums to \eqn{\sigma_p}}
#'   \item{PCR}{\eqn{\text{CR}_i / \sigma_p} — percentage contribution to risk; sums to 1}
#'   \item{Rho}{\eqn{\text{MCR}_i / \sigma_i} — correlation of asset \eqn{i} with the portfolio}
#'   \item{Beta}{\eqn{\text{PCR}_i / x_i} — beta of asset \eqn{i} with respect to the portfolio}
#' }
#'
#' @param object Object of class `portfolio`, as returned by [get_portfolio()].
#' @param W0 Numeric scalar, initial portfolio dollar value. Default is `1`.
#' @param ... Additional arguments (currently unused).
#'
#' @returns A numeric matrix with one row per asset plus a `"Portfolio"` summary
#'   row, and columns `Dollar`, `Weight`, `Vol`, `MCR`, `CR`, `PCR`, `Rho`,
#'   `Beta`. Returns invisibly.
#'
#' @examples
#' asset_names = c("MSFT", "NORD", "SBUX")
#' er = c(0.0427, 0.0015, 0.0285)
#' names(er) = asset_names
#' covmat = matrix(c(0.0100, 0.0018, 0.0011,
#'                   0.0018, 0.0109, 0.0026,
#'                   0.0011, 0.0026, 0.0199),
#'                 nrow = 3, ncol = 3,
#'                 dimnames = list(asset_names, asset_names))
#' ew = rep(1, 3) / 3
#' names(ew) = asset_names
#' port = get_portfolio(er_vec = er, cov_mat = covmat, weight_vec = ew)
#' risk_report(port, W0 = 100000)
#'
#' @export
risk_report.portfolio <- function(object, W0 = 1, ...) {
  x       <- object$weight_vec
  cov_mat <- object$cov_mat
  sd_port <- object$sd_port

  sd_assets <- sqrt(diag(cov_mat))

  d    <- x * W0
  MCR  <- as.numeric((cov_mat %*% x) / sd_port)
  names(MCR) <- names(x)
  CR   <- x * MCR
  PCR  <- CR / sd_port
  rho  <- MCR / sd_assets
  beta <- PCR / x

  report <- cbind(d, x, sd_assets, MCR, CR, PCR, rho, beta)
  colnames(report) <- c("Dollar", "Weight", "Vol", "MCR", "CR", "PCR", "Rho", "Beta")

  PORT <- c(W0, 1, NA, NA, sum(CR), sum(PCR), 1, 1)
  report <- rbind(report, PORT)
  rownames(report)[nrow(report)] <- "Portfolio"

  print(round(report, 4))
  invisible(report)
}
