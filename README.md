[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Bioc Years](https://bioconductor.org/shields/years-in-bioc/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/sargent.html)
[![Bioc Stats](https://bioconductor.org/shields/downloads/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/sargent.html)
[![Bioc Build](https://bioconductor.org/shields/build/devel/bioc/sargent.svg)](https://bioconductor.org/packages/devel/bioc/html/
sargent.html)

<p align="center" width="100%">
<img width="50%" src="vignettes/sargent-logo.png"> 
</p>

# Sargent

Identifying the cell type of origin for single cells is a key step in scRNA-seq
data analysis. Sargent is a score-based method that uses gene set of cell
type-specific markers to assign cell identities.


## Availability 

The newest stable version of **sargent** is available in [Bioconductor][bioc]. 


## Installation

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


## Getting started

Once installed the best place to get started is the [vignette][vignette].


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
[bioc]: https://bioconductor.org/packages/devel/bioc/html/sargent.html