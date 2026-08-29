
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

The package **introcompfinr** contains a few R functions for computing
Markowitz mean-variance efficient portfolios allowing for short sales
using matrix algebra computations, and not allowing short sales using
quadratic programming optimization methods. These functions allow for
the easy computation of the global minimum variance portfolio, an
efficient portfolio with a given target expected return, the tangency
portfolio, and the efficient frontier. These functions are summarized in
the table below:

| **Function** | **Description** |
|:---|:---|
| `get_portfolio` | create portfolio object |
| `globalmin_portfolio` | compute global minimum variance portfolio |
| `efficient_portfolio` | compute minimum variance portfolio subject to target return |
| `tangency_portfolio` | compute tangency portfolio |
| `efficient_frontier` | compute efficient frontier of risky assets |
