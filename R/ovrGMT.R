#' @title  Calculate over-representation of gene sets in each sample by genes
#'         from sample's largest sub-pathway
#'
#' @param  subntwlist  A list of igraph objects represented the largest 
#'                     sub-pathway for each sample. It is the output of 
#'                     `subNtw()`.
#'
#' @param  omic_genes  A vector of gene symbols to narrow down 
#'                     over-representation calculation to only those with input
#'                     genomic data. If not provided, all genes in the GMT file
#'                     will be considered. Default: NULL.
#'
#' @inheritParams subNtw
#'
#' @return  A matrix containing over-representation adjusted P with rows as 
#'          gene set names and columns as sample IDs.
#'
#' @export
#'
#' @examples
#'
#' fsubntwl  = system.file('extdata/subNtw/subntwl.rds',    package='MPAC')
#' fgmt      = system.file('extdata/ovrGMT/fake.gmt',       package='MPAC')
#' fomic_gns = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
#' subntwl  = readRDS(fsubntwl)
#' omic_gns = rownames(readRDS(fomic_gns))
#'
#' ovrGMT(subntwl, fgmt, omic_gns)
#'
#'
#' @import igraph
#' @importFrom fgsea gmtPathways
#'
ovrGMT <- function(subntwlist, fgmt, omic_genes=NULL, threads=1) {
    gmtl <- gmtPathways(fgmt)
    gmt_gns <- do.call(c, gmtl) |> unique()

    urn_balls <- NULL
    if ( is.null(omic_genes) ) {
        urn_balls <- gmt_gns
    } else {
        urn_balls <- intersect(gmt_gns, omic_genes)
    }

    getBPPARAM(threads) |> 
    bplapply(names(subntwlist), getOvrSubNtwByPat, gmtl, urn_balls, subntwlist,
        BPPARAM=_) |>
    rbindlist() |> dcast(goname ~ pat, value.var='fisher_padj') |>
    as.matrix(rownames='goname')
}

#' @importFrom stats p.adjust
getOvrSubNtwByPat <- function(pat, gmtl, urn_balls, subntwlist) {
    fisher_padj <- ipl <- name <- fisher_pval <- NULL

    subgrph <- subntwlist[[pat]]
    sel_ents <- as_data_frame(subgrph, what='vertices') |> data.table() |>
        _[ abs(ipl) > 0 ]$name

    urn_white_balls <- intersect(sel_ents, urn_balls)

    lapply(names(gmtl), defOvrByGmt, gmtl, urn_balls, urn_white_balls) |>
    rbindlist() |> _[order(fisher_pval)] |>
    _[, fisher_padj := p.adjust(fisher_pval, method='BH')] |>
    _[, pat := pat]
}

#' @importFrom stats fisher.test
defOvrByGmt <- function(goname, gmtl, urn_balls, urn_white_balls) {
    p.value <- NULL

    gmt_gns <- gmtl[[goname]]
    drawn_balls         <- intersect(gmt_gns, urn_balls)
    n_drawn_white_balls <- intersect(drawn_balls, urn_white_balls) |> length()
    n_drawn_black_balls <- length(drawn_balls) - n_drawn_white_balls
    n_not_drawn_balls   <- length(urn_balls) - length(drawn_balls)
    n_not_drawn_white_balls <- length(urn_white_balls) - n_drawn_white_balls
    n_not_drawn_black_balls <- n_not_drawn_balls - n_not_drawn_white_balls

    if ( n_drawn_white_balls > 0 ) {
        fisher_pval <- matrix(c(n_drawn_white_balls, n_not_drawn_white_balls,
            n_drawn_black_balls, n_not_drawn_black_balls),
            nrow=2, byrow=TRUE) |> fisher.test() |> _$p.value
    } else {
        fisher_pval <- NA
    }

    list( 
        goname = goname,
        n_drawn_white_gns     = n_drawn_white_balls,
        n_drawn_black_gns     = n_drawn_black_balls,
        n_not_drawn_white_gns = n_not_drawn_white_balls,
        n_not_drawn_black_gns = n_not_drawn_black_balls,
        fisher_pval = fisher_pval 
    )
}
