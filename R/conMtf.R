#' @title  Find consensus pathway motifs from a list of pathways
#'
#' @param subntwl  A list of igraph objects representing input pathways from
#'                 different samples. It is the output from `subNtw()`
#'
#' @param min_mtf_n_nodes  Number of minimum nodes in a motif. Default: 5
#'
#' @inheritParams ovrGMT
#'
#' @return  A list of igraph objects representing consensus pathway motifs
#'
#' @examples
#'
#' fsubntwl = system.file('extdata/conMtf/subntwl.rds', package='MPAC')
#' subntwl = readRDS(fsubntwl)
#'
#' fomic_gns = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
#' omic_gns = rownames(readRDS(fomic_gns))
#'
#' conMtf(subntwl, omic_gns, min_mtf_n_nodes=50)
#'
#' @export
#'
#' @import igraph
#'
conMtf <- function(subntwl, omic_genes=NULL, min_mtf_n_nodes=5) {
    . <- ent <- n_pats <- NULL

    n_subntws <- length(subntwl)
    names(subntwl) <- paste0('samp', seq_len(n_subntws))

    con_ents <- lapply(names(subntwl), function(pat) {
        data.table(pat=pat, ent=V(subntwl[[pat]])$name) %>% return()
    }) %>% rbindlist() %>%
    .[, list(n_pats = .N), by=ent] %>% .[ n_pats == n_subntws] %$% ent

    conl <- induced_subgraph(subntwl[[1]], con_ents) %>%
        decompose.graph(min.vertices=min_mtf_n_nodes)

    out_conl <- NULL
    if ( is.null(omic_genes) ) {
        out_conl <- conl 
    } else {
        out_conl <- Filter(function(grph) {
            length(intersect(V(grph)$name, omic_genes)) > 0
        }, conl) 
    }

    return(out_conl)
}
