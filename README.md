[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Coverage Status](https://img.shields.io/codecov/c/github/Oshlack/sargent/master.svg)](https://codecov.io/github/Oshlack/sargent?branch=master)
[![Bioc Years](https://bioconductor.org/shields/years-in-bioc/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/sargent.html)
[![Bioc Stats](https://bioconductor.org/shields/downloads/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/sargent.html)
[![Bioc Build](https://bioconductor.org/shields/build/devel/bioc/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/
sargent.html)

![Sargent logo](vignettes/sargent-logo.png)


# Sargent

Identifying the cell type of origin for single cells is a key step in scRNA-seq
data analysis. Sargent is a score-based method that uses gene set of cell
type-specific markers to assign cell identities.


## Availability 

The newest stable version of **sargent** is available in 
*[Bioconductor](https://bioconductor.org/packages/sargent)*. 


## Installation

Sargent is available from [Bioconductor][bioc] for R >=4.2.0.

It can be installed from Bioconductor with:

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("sargent")
```

If you wish to build a local version of the vignette use:

```{r}
BiocManager::install("sargent", build_vignettes=TRUE)
```

This will also build the vignette and install all suggested dependencies (which
aren't required for core functionality).


## Getting started

Once installed the best place to get started is the vignette. For most users
the most convenient way to access this is online [here][vignette].


## Note

This Github repository is meant mostly for development. Use at your own risk.


## Contact

For help and questions please contact the [sargent's maintenance group](mailto:ni.nouri@gmail.com).


## Citing Sargent

If you use Sargent please cite our paper ["Nouri N. et al. A marker gene-based
method for identifying the cell-type of origin from single-cell RNA sequencing
data. Genome Biology. 2023; doi:10.1186/s13059-017-1305-0"][paper].

```
  @Article{,
    author = {Nima Nouri},
    title = {A marker gene-based method for identifying the cell-type 
             of origin from single-cell RNA sequencing data},
    journal = {Genome Biology},
    year = {2023},
    url = {http://dx.doi.org/10.1186/s13059-017-1305-0},
    doi = {10.1186/s13059-017-1305-0},
  }
```

[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/sargent/inst/doc/Sargent-Vignette.html