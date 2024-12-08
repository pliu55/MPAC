#' @title  Cluster samples by pathway over-representation
#'
#' @param  ovrmat  A matrix of gene set over-representation adjusted p-values
#'                 with rows as gene sets and columns as samples. It is the
#'                 output from `ovrGMT()`.
#'
#' @param n_neighbors  Number of neighbors for clustering. A larger number is
#'                     recommended if the size of samples is large. Default: 10.
#'
#' @param n_random_runs  Number of random runs. Due to randomness introduced
#'                       to the Louvain algorithm in R igraph 1.3.0
#'                       (https://github.com/igraph/rigraph/issues/539), a large
#'                       number of runs are recommended to evaluate randomness
#'                       in the clustering results. Default: 200, which shall be
#'                       safe for sample size < 50. Please increase it
#'                       accordingly for a larger sample size.
#'
#' @inheritParams ppRnaInp
#'
#' @return  A data table with each row representing one clustering result, and
#'          the first column denotes the number of occurrences of a clustering
#'          result and the rest of columns indicating each sample's cluster
#'          index. Rows are ordered by the number of occurrences from high to
#'          low.
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
#' @importFrom BiocParallel bplapply
#'
clSamp <- function(ovrmat, n_neighbors=10, n_random_runs=200, threads=1) {
    n_top_hvgs <- 100

    ovrmat <- ovrmat[, sort(colnames(ovrmat))]
    ovrmat[is.na(ovrmat)] <- 1.0
    sce <- SingleCellExperiment(
        list(logcounts = abs(log10(ovrmat))),
        colData = DataFrame(samp=colnames(ovrmat)),
        rowData = DataFrame(goname=rownames(ovrmat))
    )

    bp <- getBPPARAM(threads)
    dec <- modelGeneVar(sce, BPPARAM=bp)
    top_hvgs <- getTopHVGs(dec, n=n_top_hvgs)
    rowData(sce)$is_top_hvgs <- (rownames(sce) %in% top_hvgs)

    sce <- denoisePCA(sce, subset.row=top_hvgs, technical=dec,
        BSPARAM=RandomParam(), BPPARAM=bp)

    ngp <- NNGraphParam(k=n_neighbors, type='jaccard', cluster.fun='louvain')

    pats <- colnames(sce)

    cldt <- bplapply(seq_len(n_random_runs), function(irep) {
        clusterCells(sce, use.dimred='PCA', BLUSPARAM=ngp) |> as.integer() |>
        data.table(irep=irep, pat=pats, icl=_)
    }, BPPARAM=bp) |> rbindlist() |>
    dcast(irep ~ pat, value.var='icl')

    nreps <- NULL
    cldt[, list(nreps = .N), by=pats] |>
    _[, c('nreps', pats), with=FALSE] |> _[order(-nreps)]
}
