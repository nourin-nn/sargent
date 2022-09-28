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


## Installation

**sargent** can be installed directly from this github with:

```{r}
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Sanofi-GitHub/PMCB-Sargent", 
                         auth_token="ask from sargent's maintenance team", 
                         build_vignettes=FALSE)
```

If you wish to build a local version of the vignette use:

```{r}
if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocStyle")

devtools::install_github("Sanofi-GitHub/PMCB-Sargent", 
                         auth_token="ask from sargent's maintenance team", 
                         build_vignettes=TRUE)
```


## Getting started

Once installed the best place to get started is the [vignette][vignette].


## Contact

For help and questions please contact the [sargent's maintenance team](mailto:nima.nouri@sanofi.com).


## Citing Sargent

If you use Sargent please cite our paper: Nouri N. et al. "A marker gene-based
method for identifying the cell-type of origin from single-cell RNA sequencing
data", Journal Name, 2022, [doi:XX.XXXX/XXXXX][paper].

```
  @Article{,
    author = {Nima Nouri},
    title = {A marker gene-based method for identifying the cell-type 
             of origin from single-cell RNA sequencing data},
    journal = {Journal Name},
    year = {2022},
    url = {http://dx.doi.org/XX.XXXX/XXXXX},
    doi = {XX.XXXX/XXXXX},
  }
```

[vignette]: https://github.com/Sanofi-GitHub/PMCB-Sargent/blob/master/vignettes/Sargent-Vignette.Rmd
[bioc]: https://bioconductor.org/packages/devel/bioc/html/sargent.html
[paper]: http://dx.doi.org/XX.XXXX/XXXXX
