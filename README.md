
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psacomp

<!-- badges: start -->
<!-- badges: end -->

*psacomp* is a package for computing *Principal Subsimplex Analysis*
introduced in *Principal Subsimplex Analysis* (2025+), Hyeon Lee, Kassel
Liam Hingee, Janice L. Scealy, Andrew T. A. Wood, Eric Grunsky, and J.
S. Marron.

## Installation

The package `psacomp` can be installed by running the following code.

``` r
library(devtools)
devtools::install_github("haneone33/psacomp")
```

## PSA Basics

First we create an example compositional data set with 7 variables and
20 data points.

``` r
set.seed(1)
X = matrix(runif(20*7),20,7)
X = sweep(X, 1, rowSums(X), '/')
```

Main function of *psacomp* is `psa()` which takes `type` and `X` as
arguments. Available options for `type` are `type = 's'` for PSA-S and
`type = 'o'` for PSA-O.

``` r
library(psacomp)
X.psas = psa(type = 's', X = X)
#> Finished creating dimension 6.
#> Finished creating dimension 5.
#> Finished creating dimension 4.
#> Finished creating dimension 3.
#> Finished creating dimension 2.
#> Finished creating dimension 1.
X.psao = psa(type = 'o', X = X)
#> Finished creating dimension 6.
#> Finished creating dimension 5.
#> Finished creating dimension 4.
#> Finished creating dimension 3.
#> Finished creating dimension 2.
#> Finished creating dimension 1.
```

## Output Details

`psa()` outputs a list of information calculated through PSA. The
information is summarized in the following table. We assume `X` is a
data matrix with $n$ rows and $D$ columns.

| Name              |    Type    |    Dimension    | Description                                                                                                               |
|-------------------|:----------:|:---------------:|---------------------------------------------------------------------------------------------------------------------------|
| Vhat              |    list    |       $D$       | `vertices$'r=i'` is an $i\times D$ matrix of vertices of representation with $i$ vertices, each row representing a vertex |
| Xhat              |    list    |       $D$       | `pts.approx$'r=i'` is an $n   imes i$ matrix of rank $i-1$ representation, in the new vertex system of `vertices$'r=i'`   |
| Xhat_reduced      |    list    |       $D$       | `pts$'r=i'` is an $n\times D$ matrix of rank $i-1$ representation, in the original vertex system of $D$ unit vertices     |
| X                 |   matrix   |   $n\times D$   | given data matrix                                                                                                         |
| residuals         |    list    |      $D-1$      |                                                                                                                           |
| scores            |   matrix   | $n\times (D-1)$ | scores                                                                                                                    |
| RSS               |   vector   |      $D-1$      | residual sum of squares                                                                                                   |
| backwards_mean    |   vector   |       $D$       | the backwards mean                                                                                                        |
| loadings          |   matrix   | $D\times (D-1)$ | a matrix of loading vectors, each column representing a loading vector                                                    |
| construction_info | data frame |   $D\times 3$   | construction information. contains index of merged vertices and merge weights                                             |
