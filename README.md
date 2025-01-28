[![bioc](http://www.bioconductor.org/shields/years-in-bioc/MPAC.svg)](http://bioconductor.org/packages/MPAC) [![DOI](https://zenodo.org/badge/664257440.svg)](https://zenodo.org/doi/10.5281/zenodo.10805479)

MPAC: inferring cancer pathway activities from multi-omic data
===========================================

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Installation](#installation)
* [Shiny App](#shinyapp)
* [Software Dependency](#softwaredependency)
* [Reference](#reference)
* [Contact](#contact)
* [License](#license)

* * *

## <a name='introduction'></a> Introduction

Multi-omic Pathway Analysis of Cell (__MPAC__) is an __R__ package that 
interprets multi-omic data through the prior knowledge of biological pathways.
The workflow of MPAC contains several steps, which are shown in 
the figure below.  MPAC has a
[vignette](https://bioconductor.org/packages/devel/bioc/vignettes/MPAC/inst/doc/MPAC.html) that describes each function in details.

<p align='center'>
    <img src="vignettes/workflow.jpg" width="800">
</p>

## <a name='installation'></a> Installation

### From GitHub

Start __R__ and enter: 

```r
devtools::install_github('pliu55/MPAC')
```

### From Bioconductor

Start __R__ and enter:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("MPAC")
```


## <a name="shinyapp"></a> Shiny App

To explore key results reported in the MPAC manuscript as well as associated
results please check out
[this R Shiny app](https://github.com/pliu55/MPAC_Shiny).


## <a name="softwaredependency"></a> Software Dependency

The `runPrd()` function requires an external software named 
[PARADIGM](https://doi.org/10.1093/bioinformatics/btq182),
which is only available for 
[Linux](
https://github.com/sng87/paradigm-scripts/tree/master/public/exe/LINUX) and 
[MacOS](
https://github.com/sng87/paradigm-scripts/tree/master/public/exe/MACOSX).
For details, please see the __Required external software__ section in vignette's __Run PARADIGM: runPrd()__.


## <a name="reference"></a> Reference

### Manuscript

__MPAC: a computational framework for inferring cancer pathway activities from multi-omic data__. Peng Liu, David Page, Paul Ahlquist, Irene M. Ong, Anthony Gitter. _bioRxiv_, 2024.06.15.599113. __doi__: [https://doi.org/10.1101/2024.06.15.599113](https://doi.org/10.1101/2024.06.15.599113)

## <a name="contact"></a> Contact

Got a question? Please report it at the [issues tab](https://github.com/pliu55/MPAC/issues) in this repository.

## <a name="license"></a> License

MPAC is licensed under the [GNU General Public License v3](LICENSE.md).
