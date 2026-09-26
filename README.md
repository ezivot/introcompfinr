
<!-- README.md is generated from README.Rmd. Please edit that file -->

# introcompfinr

<!-- badges: start -->

<!-- badges: end -->

The package **introcompfinr** accompanies the book *Introduction to
Computational Finance and Financial Econometrics with R* by Eric Zivot.
The package contains functions for portfolio analysis, risk management,
and financial econometrics.

## Installation

You can install the development version of **introcompfinr** from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ezivot/introcompfinr")
```

## Usage

The package **introcompfinr** contains

- All of the data sets used in the book
- Functions for estimating econometric models for asset returns
- Functions for estimating risk and performance measures
- Functions for portfolio analysis and risk management

The package is designed to be used with the book, but can also be used
independently.

## Functions for estimating econometric models for asset returns

The functions for estimating econometric models for asset returns
summarized in the table below:

| **Function** | **Description** |
|:---|:---|
| `estimate_gwn_model` | Estimate the Gaussian white noise model for asset returns |
| `estimate_si_model` | Estimate Sharpe’s Single Index (Market) model for asset returns |
| `estimate_factor_model` | Estimate linear factor model for asset returns |

Note: `estimate_si_model` and `estimate_factor_model` are not yet
implemented in the package.

## Function for estimating risk and performance measures

The functions for estimating risk and performance measures summarized in
the table below:

| **Function**   | **Description**                                |
|:---------------|:-----------------------------------------------|
| `estimate_var` | Estimate value at risk (VaR) for asset returns |
| `estimate_sr`  | Estimate Sharpe ratio (SR) for asset returns   |
| `risk_report`  | Estimate portfolio risk report                 |

## Portfolio analysis and risk management functions

The package contains a few R functions for computing Markowitz
mean-variance efficient portfolios allowing for short sales using matrix
algebra computations, and not allowing short sales using quadratic
programming optimization methods. These functions allow for the easy
computation of the global minimum variance portfolio, an efficient
portfolio with a given target expected return, the tangency portfolio,
and the efficient frontier. These functions are summarized in the table
below:

| **Function** | **Description** |
|:---|:---|
| `get_portfolio` | create portfolio object |
| `globalmin_portfolio` | compute global minimum variance portfolio |
| `efficient_portfolio` | compute minimum variance portfolio subject to target return |
| `tangency_portfolio` | compute tangency portfolio |
| `efficient_frontier` | compute efficient frontier of risky assets |
