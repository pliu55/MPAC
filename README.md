<!--
![bioc](http://www.bioconductor.org/shields/years-in-bioc/MPAC.svg)](http://bioconductor.org/packages/devel/bioc/html/MPAC.html)
-->

MPAC: Inferring Cancer Pathway Activities from Multi-omic data
===========================================

Table of Contents
-----------------

* [Introduction](#Introduction)
* [Installation](#Installation)
* [Reference](#Reference)
* [Shiny App](#ShinyApp)
* [Contact](#Contact)
* [License](#License)

* * *

## <a name='Introduction'></a> Introduction

Multi-omic Pathway Analysis of Cancer (__MPAC__) is an __R__ package that 
interprets multi-omic data through the prior knowledge of biological pathways.
The workflow of MPAC contains several steps, which is shown in 
the figure below with function names and associated key parameters.  PRAM has a
[vignette TO UPDATE!!](https://bioconductor.org/packages/devel/bioc/vignettes/pram/inst/doc/pram.pdf) that describes each function in details.

<p align='center'>
    <img src="vignettes/workflow.jpg" width="700">
</p>

## <a name='Installation'></a> Installation

### From GitHub

Start __R__ and enter: 

```r
devtools::install_github('pliu55/MAPC')
```

### From Bioconductor

Start __R__ and enter:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MPAC")
```


## <a name="Reference"></a> Reference

__MAPC: a computational framework for inferring cancer pathway activities from multi-omic data__. Peng Liu, David Page, Paul Ahlquist, Irene M. Ong, Anthony Gitter. _bioRxiv_, 2023. __doi__: TO-UPDATE!!


## <a name="ShinyApp"></a> Shiny App

To explore key results reported in the MAPC manuscript as well as associated
results please check out 
[an R Shiny app](https://github.com/pliu55/MPAC_Shiny).


## <a name="Contact"></a> Contact

Got a question? Please report it at the [issues tab](https://github.com/pliu55/MPAC/issues) in this repository.

## <a name="License"></a> License

PRAM is licensed under the [GNU General Public License v3](LICENSE).
