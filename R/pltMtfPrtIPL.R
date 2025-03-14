#' @title  Plot a heatmap of IPLs on proteins from consensus pathway submodules
#'
#' @param grphl  A list of igraph objects representing consensus pathway
#'               submodules. It is the output from `conMtf()`.
#'
#' @param proteins  A vector of proteins, of which IPLs to plot. Default: all
#'                  proteins that in both `grphl` and `fltmat`.
#'
#' @inheritParams pltNeiStt
#' @inheritParams pltOvrHm
#'
#' @return  A Kaplan-Meier plot
#'
#' @examples
#'
#' fltmat <- system.file('extdata/pltSttKM/ipl.rds',package='MPAC') |> readRDS()
#' cldt <- system.file('extdata/pltMtfPrtIPL/cl.rds',package='MPAC')|> readRDS()
#' grphl <- system.file('extdata/pltMtfPrtIPL/grphl.rds',package='MPAC') |>
#'          readRDS()
#'
#' pltMtfPrtIPL(fltmat, cldt, grphl, proteins=c('CD247', 'FASLG'))
#'
#' @export
#'
#' @import igraph
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#'
pltMtfPrtIPL <- function(fltmat, cldt, grphl, proteins=NULL){
    protdt = findJntProts(fltmat, grphl, proteins)

    nreps = NULL
    setnames(cldt, 1, 'nreps')
    srt_cldt = cldt[order(-nreps)][1, ] |>
               melt(id='nreps', variable='brc', value='icl')

    jnt_brcs = intersect(srt_cldt$brc, colnames(fltmat))
    brc = NULL
    sel_srt_cldt = srt_cldt[ brc %in% jnt_brcs ]
    pltmat = fltmat[protdt$prot, sel_srt_cldt$brc, drop=FALSE]

    hmMtfIplByCl(pltmat, protdt, sel_srt_cldt)
}

hmMtfIplByCl <- function(pltmat, protdt, cldt) {
    icl_lab = icl = NULL
    cldt[, icl_lab := paste0('c', icl, "\nn=", .N), by=icl]

    FONT_SIZE = 9
    CYAN    = '#56B4E9'
    MAGENTA = '#CC79A7'
    Heatmap( pltmat,
        col = circlize::colorRamp2(c(-0.5, 0, 0.5), c(CYAN, 'gray95', MAGENTA)),
        na_col = 'black',
        border            = 'black',
        rect_gp           = gpar(col='white'),
        column_split      = cldt$icl_lab,
        row_split         = protdt$mtf,
        show_row_names    = TRUE,
        show_column_names = FALSE,
        column_names_side = 'top',
        row_names_gp      = gpar(fontsize=FONT_SIZE),
        column_title_gp   = gpar(fontsize=FONT_SIZE),
        row_title_rot     = 0,
        row_title_gp      = gpar(fontsize=FONT_SIZE),
        show_row_dend     = FALSE,
        show_column_dend  = FALSE,
        cluster_column_slices = FALSE,
        cluster_row_slices    = FALSE,
        cluster_columns   = TRUE,
        cluster_rows      = FALSE,
        use_raster        = FALSE,
        heatmap_legend_param = list(
            title       = 'IPL',
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            grid_height = grid::unit(3, 'mm'),
            grid_width  = grid::unit(3, 'mm'),
            border      = 'black',
            title_position = 'lefttop',
            direction   = 'horizontal'
        )
    ) |> draw(
        heatmap_legend_side = 'bottom',
        column_title = paste0(ncol(pltmat), ' samples'),
        column_title_gp = gpar(fontsize=FONT_SIZE),
        background = 'transparent'
    )
}

findJntProts <- function(fltmat, grphl, proteins) {
    if ( is.null(proteins) ) {
        proteins = rownames(fltmat)
    }

    grph_nodes = do.call(igraph::union, grphl) |> V() |> _[]$name
    prots = rownames(fltmat) |> intersect(grph_nodes) |> intersect(proteins)
    lapply(seq_along(grphl), function(isub) {
        data.table(mtf = paste0('submodule ', isub),
                   prot = intersect(prots, V(grphl[[isub]])$name))
    }) |> rbindlist()
}
