#' Create four-panel diagnostic plot of returns
#'
#' Four-panel plot showing histogram with normal curve overlaid, boxplot, sample ACF, and normal
#' qq-plot. Used for accessing the appropriateness of the Gaussian White Noise model for
#' returns. Inspired by the four-panel plot in *Statistical Analysis of Financial Data with R* by
#' Rene Carmona.
#'
#' @param ret Single data object (xts, zoo, matrix, data.frame) of returns. It is assumed that the
#' column has a name
#'
#' @importFrom graphics boxplot hist legend lines par
#'
#' @importFrom stats acf density dnorm qqnorm qqline sd
#'
#' @export
#'
#' @examples
#' examplePrices=amznDailyPrices["1998::2002"]
#' exampleRetS = PerformanceAnalytics::Return.calculate(examplePrices,
#'                                     method="discrete")
#' exampleRetS = na.omit(exampleRetS)
#' four_panel_plot(exampleRetS)
four_panel_plot <-
  function(ret) {
    ret = PerformanceAnalytics::checkData(ret, "matrix")
    if(ncol(ret) > 1) stop("four_panel_plot() only works for a single column of returns")
    retName = colnames(ret)
    if (is.null(retName)) retName = "Asset"
    ret.den = density(ret)
    par(mfrow=c(2,2))
    hist(ret, main=paste(retName, "Returns", sep=" "),
         xlab=retName, probability=T, col="cornflowerblue")
    # overlay normal distribution on smoothed density
    lines(ret.den$x, dnorm(ret.den$x, mean=mean(ret), sd=sd(ret)),
          col="black", lwd=2)
    legend(x="topleft", legend=c("Normal Curve"),
           lty=1, col="black", lwd=2, bty="n")
    boxplot(ret, main="Boxplot of Returns", outchar=T, col="cornflowerblue")
    # autocorrelations
    acf(ret, lwd=2, main="ACF of Returns", col="cornflowerblue")
    # qq plot
    qqnorm(ret, col="cornflowerblue", pch=16)
    qqline(ret)
    par(mfrow=c(1,1))
  }
