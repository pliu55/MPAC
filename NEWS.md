# MPAC 1.1.6

- updated pseudocode in vignettes
- added new functions
  - `pltConMtf()`: to plot consensus pathway submodules
  - `pltMtfPrtIPL()`: to plot a heatmap of IPLs on proteins from consensus
    pathway submodules
  - `pltSttKM()`: to plot Kaplan-Meier curve for patient samples stratified
    by given protein(s)' pathway states

# MPAC 1.1.5

- use 'Multi-omic Pathway Analysis of Cells' for MPAC

# MPAC 1.1.4

- updated text in DESCRIPTION

# MPAC 1.1.3

- Vignettes
  - added examples on `pltOvrHm()`
  - added pseudocode for MPAC
  - changed MPAC's full name

# MPAC 1.1.2

- pltOvrHm(): use lowercase 'p' instead of 'P' in heatmap colorbar legend
- MPAC's full name is changed to 'Multi-omic Pathway Analysis of Cell'

# MPAC 1.1.1

- added new functions
  - `pltOvrHm()`: to plot a heatmap of over-represented gene sets for clustered
    samples
  - `ppRunPrd()`: to prepare required files to run PARADIGM
- updated existing functions
  - `conMtf()`: use `decompose()` instead of `decompose.graph()` for igraph 2.0
  - `ppPermInp()`: increased default permutations from 3 to 100
  - `ppRealInp()` and `ppRnaInp()`: added an option `rna_n_sd` to set standard
    deviation range to define RNA state
  - `pltNeiStt()`: fixed bugs and adjusted output figure height
  - `colPermIPL()` & `fltByPerm()`: added a `threads` option to parallelize
  - `runPrd()` & `runPermPrd()`: when no `sampleids` is provided, all the  
    samples in input `real_se` and `perml` will be used.
  - `subNtw()`: entities with NA values in all samples in the input `fltmat` 
    will be ignored.
- updated tests
  - `test-conMtf()` & `test-subNtw()`: use `as_edgelist()` instead of
    `get.edgelist()` for igraph 2.0
  - `test-conMtf()` & `test-ovrGMT()`: add `upgrade_graph()` on extdata


# MPAC 0.99.10

Added documentation in inst/scripts/ for the source and derivation of 
inst/extdata/

# MPAC 0.99.1

- Used SummarizedExperiment objects for CNA and RNA matrix in functions:   
  - ppCnInp()
  - ppRnaInp()
  - ppPermInp()
  - runPrd()
  - runPermPrd()
  - pltNeiStt()

- Added a function to prepare both CNA and RNA input: ppRealInp()
- Removed GitHub PARADIGM URLs from R code
- Removed dependence on the magrittr package

# MPAC 0.99.0

* Added a `NEWS.md` file to track changes to the package.
