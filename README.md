OuRS: Outliers and Robust Structures
================

OuRS
====

`OuRS` is short for *Ou*tliers and *R*obust *S*tructures. The goal of this package is to provide a variety of multivariate outlier detection and robust subspace techniques. While there are several other packages to do so in `R`, this package was developed to meet specific needs not currently met including, for examples, outlier detection in multivariate categorical or mixed data, and simplified approaches to high dimensional (i.e., *n* &lt; &lt;*P*) data.

As of now, the OuRS package is under active development. However, several approaches we make available in `OuRS` are stable and available for use, namely:

-   Minimum covariance determinant (MCD)

-   A generalized MCD ("GenMCD") for almost any data type

-   Principal components analysis + split-half resampling (PCA+SHR)

Installation
------------

Because `OuRS` is under development it is not yet available as a complete package, though it can still be installed. `OuRS` depends on the `GSVD` package ([GSVD](https://github.com/derekbeaton/ExPosition-Family/tree/master/ExPosition2/GSVD)). For now, the simplest approach to installation is through the `devtools` package:

``` r
devtools::install_github("derekbeaton/ExPosition-Family",subdir="/ExPosition2/GSVD/Package")
devtools::install_github("derekbeaton/OuRS",subdir="/OuRS")
```

Examples
========

Here we provide examples with brief explanations of the three previously mentioned approaches: MCD, GenMCD, and PCA+SHR.

MCD
---

The MCD algorithm () is an outlier detection and robust covariance technique. The MCD's goal is to find a covariance matrix with a robust location (mean) and minimum scatter (covariance). The objective of the MCD is to find a minimum determinant (i.e., product of the eigenvalues), which reflects a minimum scatter. The MCD primarily relies on *subsampling* to find some subsample *H* that provides a minimum determinant. Generally the MCD works as follows. Given some matrix **X** with *I* rows and *J* columns, we randomly subsample **X** to have *H* rows, referred to as **X**<sub>*H*</sub>. The size of *H* is controlled by the user through the *α* parameter which must exist between .5 and 1:

1.  For **X**<sub>*H*</sub>, compute the column-wise mean as **μ**<sub>*H*</sub> and the covariance as **S**<sub>*H*</sub>.

2.  Compute the squared Mahalanobis distances for each observation in **X** as *m*<sub>*i*</sub> = (*x*<sub>*i*</sub> − **μ**<sub>*H*</sub>)**S**<sub>*H*</sub><sup>−1</sup>(*x*<sub>*i*</sub> − **μ**<sub>*H*</sub>)<sup>*T*</sup>.

3.  Set **X**<sub>*H*</sub> as the subset of *H* observations with the smallest Mahalanobis distances computed from Step 2.

4.  Compute and store the determinant of **S**<sub>*H*</sub> which is the covariance matrix of **X**<sub>*H*</sub>.

The steps above are repeated for either a pre-specified number of iterations or until convergence (i.e., the determinant in Step 4 has little to no change). Once completed, the *H* set with a minimum determinant reflects the sub-sample with the robust location (robust mean), minimum scatter (robust covariance) from the data in **X**. Furthermore, a set of robust Mahalanobis distances are computed from the robust covariance matrix. The robust Mahalanobis distances can then be used to identify outliers.

Below we provide an example of the MCD for a routinely used data set for the MCD: the `philips` dataset from the `cellWise` package. While there are numerous outputs from the algorithm, we only show a graph of the standard vs. robust Mahalanobis distances for the `philips` data.

    ## Loading required package: GSVD

    ## Loading required package: ours

    ## Loading required package: cellWise

    ## Warning: package 'cellWise' was built under R version 3.4.4

``` r
  # default parameters for now.
mcd.results <- cont.mcd(philips)
plot(mcd.results$dists$md,mcd.results$dists$rob.md,pch=20,col="mediumorchid4")
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

GenMCD
------

PCA+SHR
-------
