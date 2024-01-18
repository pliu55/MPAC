<!--
![bioc](http://www.bioconductor.org/shields/years-in-bioc/MPAC.svg)](http://bioconductor.org/packages/devel/bioc/html/MPAC.html)
-->

MPAC: inferring cancer pathway activities from multi-omic data
===========================================

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Installation](#installation)
* [Shiny App](#shinyApp)
* [Reference](#reference)
* [Contact](#contact)
* [License](#license)

* * *

## <a name='introduction'></a> Introduction

Multi-omic Pathway Analysis of Cancer (__MPAC__) is an __R__ package that 
interprets multi-omic data through the prior knowledge of biological pathways.
The workflow of MPAC contains several steps, which are shown in 
the figure below.  MPAC has a
[vignette TO UPDATE!!](https://bioconductor.org/packages/devel/bioc/vignettes/pram/inst/doc/pram.pdf) that describes each function in details.

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
BiocManager::install("MPAC")
```


## <a name="reference"></a> Reference

__MPAC: a computational framework for inferring cancer pathway activities from multi-omic data__. Peng Liu, David Page, Paul Ahlquist, Irene M. Ong, Anthony Gitter. _bioRxiv_, 2024. __doi__: TO-UPDATE!!


## <a name="shinyApp"></a> Shiny App

To explore key results reported in the MPAC manuscript as well as associated
results please check out
[this R Shiny app](https://github.com/pliu55/MPAC_Shiny).


## <a name="contact"></a> Contact

Got a question? Please report it at the [issues tab](https://github.com/pliu55/MPAC/issues) in this repository.

## <a name="license"></a> License

MPAC is licensed under the [GNU General Public License v3](LICENSE.md).
