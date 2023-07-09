#' @title  Cluster samples by pathway over-representation
#'
#' @param  ovrmat  A matrix of gene set over-representation adjusted p-values
#'                 with rows as gene sets and columns as samples. It is the 
#'                 output from `ovrGMT()`.
#'
#' @param n_neighbors  Number of neighbors for clustering. A larger number is
#'                     recommended the size of samples is large. Default: 10.
#'
#' @inheritParams ppRnaInp
#'
#' @return  A matrix of sample clustering result with rows as samples and a 
#'          column of cluster index
#'
#' @examples
#'
#' fovr = system.file('extdata/clSamp/ovrmat.rds', package='MPAC')
#' ovrmat = readRDS(fovr)
#'
#' clSamp(ovrmat)
#'
#' @export
#'
#' @import SingleCellExperiment
#' @importFrom scran modelGeneVar getTopHVGs denoisePCA clusterCells
#' @importFrom BiocSingular RandomParam
#' @importFrom bluster NNGraphParam
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom S4Vectors DataFrame
#'
clSamp <- function(ovrmat, n_neighbors=10, threads=1) {
    . = NULL
    n_top_hvgs = 100

    ovrmat[is.na(ovrmat)] = 1.0
    sce = SingleCellExperiment(
        list(logcounts = abs(log10(ovrmat))),
        colData = DataFrame(samp=colnames(ovrmat)), 
        rowData = DataFrame(goname=rownames(ovrmat))
    )

    dec = modelGeneVar(sce)
    top_hvgs = getTopHVGs(dec, n=n_top_hvgs)
    rowData(sce)$is_top_hvgs = (rownames(sce) %in% top_hvgs)

    sce = denoisePCA(sce, subset.row=top_hvgs, technical=dec,
        BSPARAM=RandomParam())

    NNGraphParam(k=n_neighbors, type='jaccard', cluster.fun='louvain') %>%
    clusterCells(sce, use.dimred='PCA', BLUSPARAM=.) %>% as.integer() %>%
    matrix(ncol=1, dimnames=list(colnames(sce), c('icl'))) %>%
    return()
}
