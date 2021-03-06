OuRS: Outliers and Robust Structures
================

# OuRS

`OuRS` is short for *Ou*tliers and *R*obust *S*tructures. The goal of
this package is to provide a variety of multivariate outlier detection
and robust subspace techniques. While there are several other packages
to do so in `R`, this package was developed to meet specific needs not
currently met by those packages. *OuRS* includes, for examples, outlier
detection in multivariate categorical or mixed data, and simplified
approaches to high dimensional (i.e., \(n << P\)) data, as well as
approaches to identify sources (variables) of outlierness (of
observations).

As of now, the OuRS package is under active development. However,
several approaches we make available in `OuRS` are stable and available
for use, namely:

  - Minimum covariance determinant (MCD)

  - A generalized MCD (“GenMCD”) for almost any data type

  - The “CorrMax” transformation (a.k.a. the Garthwaite-Koch partition)

  - A generalized approach to the “CorrMax” transformation

  - Principal components and correspondence analyses (PCA, CA)

  - Temporarily removed and to be re-introduced:
    
      - Two-fold repeated PCA and CA to identify outliers in high
        dimensional data
    
      - Robust PCA/SVD (a la Candes)

## Installation

Because `OuRS` is under development it is not yet available as a
complete package, though it can still be installed. `OuRS` depends on
the `GSVD` package ([GSVD](https://github.com/derekbeaton/GSVD)). For
now, the simplest approach to installation is through the `devtools` or
`remote` packages. **NOTE** both packages are more aggressive about
warnings and treat them as errors, and thus prevents installation (see
[here](https://github.com/r-lib/remotes/issues/403), for example). So to
install `OuRS`, you will need to set an environment variable for the
`remotes` package (there may be equivalents in `devtools` if that is
preferred):

``` r
remotes::install_github("derekbeaton/GSVD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
remotes::install_github("derekbeaton/OuRS/OuRS")
```

# Set up

Prior to running the examples, we require several packages to be loaded.
The analyses are all performed in `OuRS`, which depends on `GSVD`. We
also require `cellWise` to access data for this README.

# Examples

Here we provide examples with brief explanations of the three previously
mentioned approaches: MCD and a GenMCD.

## MCD

The [MCD algorithm (Hubert &
Debruyne, 2009)](https://onlinelibrary.wiley.com/doi/full/10.1002/wics.61)
is an outlier detection and robust covariance technique. The MCD’s goal
is to find a covariance matrix with a robust location (mean) and minimum
scatter (covariance). The objective of the MCD is to find a minimum
determinant (i.e., product of the eigenvalues), which reflects a minimum
scatter. The MCD primarily relies on *subsampling* to find some
subsample \(H\) that provides a minimum determinant. Generally the MCD
works as follows. Given some matrix \(\mathbf{X}\) with \(I\) rows and
\(J\) columns, we randomly subsample \(\mathbf{X}\) to have \(H\) rows,
referred to as \(\mathbf{X}_{H}\). The size of \(H\) is controlled by
the user through the \(\alpha\) parameter which must exist between .5
and 1:

1.  For \(\mathbf{X}_{H}\), compute the column-wise mean as
    \(\boldsymbol{\mu}_{H}\) and the covariance as \(\mathbf{S}_{H}\).

2.  Compute the squared Mahalanobis distances for each observation in
    \(\mathbf{X}\) as
    \(m_{i} = (x_{i} - \boldsymbol{\mu}_{H}) \mathbf{S}_{H}^{-1} (x_{i} - \boldsymbol{\mu}_{H})^{T}\).

3.  Set \(\mathbf{X}_{H}\) as the subset of \(H\) observations with the
    smallest Mahalanobis distances computed from Step 2.

4.  Compute and store the determinant of \(\mathbf{S}_{H}\) which is the
    covariance matrix of \(\mathbf{X}_{H}\).

The steps above are repeated for either a pre-specified number of
iterations or until convergence (i.e., the determinant in Step 4 has
little to no change). Once completed, the \(H\) set with a minimum
determinant reflects the sub-sample with the robust location (robust
mean), minimum scatter (robust covariance) from the data in
\(\mathbf{X}\). Furthermore, a set of robust Mahalanobis distances are
computed from the robust covariance matrix. The robust Mahalanobis
distances can then be used to identify outliers.

Below we provide an example of the MCD for a routinely used data set for
the MCD: the `philips` dataset from the `cellWise` package. While there
are numerous outputs from the algorithm, we only show a graph of the
standard vs. robust Mahalanobis distances for the `philips` data. The
`philips` data are used as a typical example because of the “masking”
effect, wherein there are outliers that exist but are not detectable by
Mahalanobis distance. This is because they outlyingness of these
observations exist, typically, in a lower variance subspace. Thus the
MCD and robust Mahalanobis distance “unmasks” these observations.

``` r
data("philips")
# mostly default parameters for now, but with the lowest size
# of H, i.e., alpha=.5, and many iterations
mcd.results <- continuous_mcd(philips, alpha = 0.5)
dd_plot(mcd.results, dist_transform = "sqrt")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

As can be seen in the above figure, there are observations with a large
robust Mahalanobis distance, but relatively small Mahalanobis distance;
these are the observations that were “masked”.

## GenMCD

One of the primary motivations behind the creation of the `OuRS` package
was the development of an extension to the MCD. We initially extended
the MCD so that it could be applied to categorical (and similar data).
However, our extension of the MCD can incorporate data of almost any
data type (categorical, ordinal, continuous, counts), either as a
homogeneous set or heterogeneous variables (i.e., mixed data). We have
made a [preprint](https://www.biorxiv.org/content/10.1101/333005v2) of
this available on bioRxiv.

The example we provide here is strictly categorical data. The data here
are a toy data set of genetic variables (single nucleotide
polymorphisms, a.k.a. SNPs).

``` r
load("OuRS/data/SNPS.rda")
summary(SNPS[, 1:2])
```

    ##     SNP.1              SNP.2          
    ##  Length:60          Length:60         
    ##  Class :character   Class :character  
    ##  Mode  :character   Mode  :character

``` r
# mostly default parameters for now, but with the lowest size
# of H, i.e., alpha=.5, and many iterations
catmcd.results <- categorical_mcd(SNPS, alpha = 0.5)
dd_plot(catmcd.results, dist_transform = "sqrt")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Additional materials

Please see the (Publications)\[./Publications\] and
(Presentations)\[./Presentations\] directories for more examples of some
of our outlier and robust structures work.
