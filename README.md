
<!-- README.md is generated from README.Rmd. Please edit that file -->

# levelnet

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/schochastics/levelnet.svg?branch=master)](https://travis-ci.org/schochastics/levelnet)

levelnet is, so far, an early-stage R package that can be used to
analyse two-mode networks and specifically their one-mode projection.
The main purpose is to make several methods available that extract the
*binary backbone* of a one-mode projection. The best documented and
tested method is the *stochastic degree sequence model* (SDSM) by Z.
Neal
([link](https://www.sciencedirect.com/science/article/pii/S0378873314000343)).
It also implements the fairly unknown (but in the case of SDSM useful)
scobit model.

The package still contains a lot of undocumented functions, which can be
ignored for now.

## Installation

You can install the developers version of levelnet with:

``` r
# install.packages("devtools")
devtools::install_github("schochastics/levelnet")
```

## Example of `sdsm`

The package includes bill-cosponsorship data of the 115th Senate.
Similar data was used in the Paper by Neal an in [follow up
work](http://www.sciencedirect.com/science/article/pii/S0378873317303039).
The function `bipartite_from_data_frame` can be used to turn the data
frame into a two-mode network (bipartite graph).

``` r
library(igraph)
library(levelnet)
data("cosponsor_senate_15")
head(cosponsor_senate_15)
#>                  name party     bill
#> 1    Enzi, Michael B.     R    115s1
#> 2 Cardin, Benjamin L.     D   115s10
#> 3    Wicker, Roger F.     R   115s10
#> 4    Alexander, Lamar     R  115s100
#> 5         Franken, Al     D 115s1000
#> 6       Murray, Patty     D 115s1000

g <- bipartite_from_data_frame(cosponsor_senate_15,"name","bill")
```

The function `sdsm_diagnostic` checks the performance of several link
functions. The default optimization method for the MLE of the scobit
model is “BFGS”. ~~However, I noticed that it often runs into
convergence issues (like for this dataset). A more robust method is
“Nelder-Mead” which is a little bit more time consuming. You may still
run into convergence problems. If so, try to adjust the parameter
vector.~~

**Version 0.2.0 implements the gradients of the scobit function
explicitly which facilitates the optimization**.

``` r
params <- list(b0 = 1e-5,b1 = 1e-5,b2 = 1e-5,b3 = 1e-5,a = 0.8)
sdsm_diagnostic(g,verbose = FALSE,params = params)
```

<img src="man/figures/README-diagnostics-1.png" width="100%" />

    #>      name rmse_row rmse_col   time
    #> 1   logit 32.49987 3.539514  1.977
    #> 2  probit 28.82683 3.265757  1.721
    #> 3 cloglog 38.65148 4.156065  5.482
    #> 4  scobit 14.45308 2.731727 35.594

Note that there is no significant speed up for the scobit model to the
previous version (runtime was ~50 sec). The gradient function is rather
complicated and needs to be evaluated many times during the optimization
process. This and the scobit function evaluation are the main
bottleneck.

As was noted in the paper, the scobit model produces the best fit of the
data.

``` r
params <- list(b0 = 1e-5,b1 = 1e-5,b2 = 1e-5,b3 = 1e-5,a = 0.8)
l <- sdsm(g,proj = "true",model = "scobit",params = params)
```

<img src="man/figures/README-graph.png" width="80%" style="display: block; margin: auto;"/>
