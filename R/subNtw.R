#' @title  Subset pathways by IPL results
#'
#' @param fltmat  A matrix contains filterd IPL with rows
#'                as 'entity' and column as samples. This is the output from
#'                `fltByPerm()`.
#'
#' @param fgmt  A gene set GMT file. This will be the same file used for the
#'              gene set over-representation calculation in the next step. It
#'              is used here to ensure output sub-pathway contains a minimum 
#'              number of genes from to-be-used gene sets.
#'
#' @param min_n_gmt_gns  Minimum number of genes from the GMT file in the output
#'                       sub-pathway. Default: 2.
#'
#' @inheritParams ppRnaInp
#' @inheritParams runPrd
#'
#' @return  A list of igraph objects representing the largest sub-pathway
#'          for each sample.
#'
#' @export
#'
#' @examples
#'
#' fflt = system.file('extdata/fltByPerm/flt_real.rds', package='MPAC')
#' fltmat = readRDS(fflt)
#' fpth = system.file('extdata/Pth/tiny_pth.txt',       package='MPAC')
#' fgmt = system.file('extdata/ovrGMT/fake.gmt',        package='MPAC')
#'
#' subNtw(fltmat, fpth, fgmt, min_n_gmt_gns=1)
#'
#' 
#' @importFrom fgsea gmtPathways
#'
subNtw <- function(fltmat, fpth, fgmt, min_n_gmt_gns=2, threads=1) {
    pthl <- getNodeEdge(fpth)
    nodedt <- pthl$nodedt
    edgedt <- pthl$edgedt

    gmt_gns <- gmtPathways(fgmt) |> do.call(c, args=_) |> unique()
    sampleids <- colnames(fltmat)

    outl <- getBPPARAM(threads) |> 
        bplapply(sampleids, getSubNtwByPat, nodedt, edgedt, fltmat, gmt_gns,
        min_n_gmt_gns, BPPARAM=_)

    names(outl) <- sampleids
    return(outl)
}

#' @import igraph
#'
getSubNtwByPat <- function(pat, in_nodedt, in_edgedt, fltmat, gmt_gns, 
    min_n_gmt_gns) {

    nents <- from <- to <- entity <- isub <- NULL

    ipls <- fltmat[, pat]
    ents <- ipls[ ! is.na(ipls) ] |> names()

    ipldt <- fltmat[, pat, drop=FALSE] |>
        as.data.table(keep.rownames='entity') |> setnames(pat, 'ipl')
    edgedt <- in_edgedt[(from %in% ents) & (to %in% ents)]
    nodedt <- in_nodedt[(entity %in% edgedt$from) | (entity %in% edgedt$to)] |>
        merge(ipldt, by='entity', all.x=TRUE)

    subl <- graph_from_data_frame(d=edgedt, directed=TRUE, vertices=nodedt) |>
        decompose() |>
        Filter(function(g) {
            n_gmt_gns <- intersect(V(g)$name, gmt_gns) |> length()
            return(n_gmt_gns >= min_n_gmt_gns)
        }, x=_)

    ndt <- length(subl) |> seq_len() |> lapply(function(isub) 
        list(isub=isub, nents=length(V(subl[[isub]])$name))) |> rbindlist()

    subntw <- subl[[ ndt[ nents == max(nents)]$isub ]]
    return(subntw)
}
