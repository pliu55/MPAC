#' @title  Plot consensus pathway submodules
#'
#' @param  proteins  A vector of protein symbols to highlight in the plot.
#'                   Default: no protein will be highlighted.
#'
#' @inheritParams pltMtfPrtIPL
#'
#' @return  a plot of consensus pathway submodules
#'
#' @examples
#'
#' grphl <- system.file('extdata/pltMtfPrtIPL/grphl.rds',package='MPAC') |>
#'          readRDS()
#'
#' proteins <- system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC') |>
#'             readRDS() |> rownames() |> c('CD3G')
#'
#' pltConMtf(grphl, proteins) |> print()
#'
#' @export
#'
#' @import igraph
#' @import ggraph
#' @importFrom stringr str_wrap
#'
pltConMtf <- function(grphl, proteins=NULL){
    hasomic = lab = NULL

    ntw = do.call(igraph::union, grphl)
    V(ntw)$hasomic = (V(ntw)$name %in% proteins)
    V(ntw)$lab     = gsub('_', ' ', V(ntw)$name) |> stringr::str_wrap(width=20)

    set.seed(88888888)
    ggraph(ntw, layout='fr') +
    geom_edge_link(width=0.2, color='gray80',
                   arrow=arrow(length=unit(1, 'mm')),
                   start_cap=circle(1, 'mm'), end_cap=circle(1, 'mm')) +
    geom_node_point(aes(color=hasomic), size=1) +
    geom_node_text(aes(label=lab, color=hasomic), repel=TRUE, size=2.5) +
    theme_void() +
    theme( legend.position = 'none' ) +
    scale_color_manual( name   = 'protein with\nomic data',
                        values = c( 'TRUE'  = 'firebrick1',
                                    'FALSE' = 'gray10' ) )
}
