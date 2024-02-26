#' @title  Plot a heatmap of pathway and omic states of a protein and its
#'         pathway neighbors
#'
#' @param  protein  Name of the protein to plot. It requires to have CN and RNA
#'                  state data, as well as pathway data from the input.
#'
#' @inheritParams runPrd
#' @inheritParams subNtw
#'
#' @return  A heatmap of pathway and omic states of a protein and its
#'          pathway neighbors
#'
#' @examples
#'
#' fpth = system.file('extdata/Pth/tiny_pth.txt', package='MPAC')
#'
#' fcn  = system.file('extdata/pltNeiStt/inp_focal.rds',       package='MPAC')
#' frna = system.file('extdata/pltNeiStt/inp_log10fpkmP1.rds', package='MPAC')
#' fflt = system.file('extdata/pltNeiStt/fltmat.rds',          package='MPAC')
#'
#' cn_state_mat  = readRDS(fcn)
#' rna_state_mat = readRDS(frna)
#' fltmat = readRDS(fflt)
#' protein = 'CD86'
#'
#' pltNeiStt(cn_state_mat, rna_state_mat, fltmat, fpth, protein)
#'
#' @export
#'
pltNeiStt <- function(cn_state_mat, rna_state_mat, fltmat, fpth, protein) {
    nodedt <- ppNode4PltNeiStt(fpth, protein)
    pltmat <- ppMat4PltNeiStt(cn_state_mat, rna_state_mat, fltmat, nodedt)

    id <- NULL
    plt_nodedt <- nodedt[ id %in% rownames(pltmat) ]

    makeHmNeiStt(pltmat, plt_nodedt)
}

#' @import ComplexHeatmap
#' @importFrom grid  gpar unit
#'
makeHmNeiStt <- function(pltmat, nodedt) {
    level <- type <- id <- row_grp <- NULL
    nodedt <- nodedt[order(level, type, id)]
    nodedt[, row_grp := factor(row_grp, levels=unique(nodedt$row_grp))]
    pltmat <- pltmat[ nodedt$id, ]
    Heatmap(pltmat,
        col = structure(c(MAGENTA, 'gray85', CYAN), names=c(1, 0, -1)),
        rect_gp            = gpar(col='white'),
        border             = 'black',
        row_split          = nodedt$row_grp,
        cluster_row_slices = FALSE,
        row_gap            = unit(0.2, 'cm'),
        row_title_rot      = 0,
        row_title_gp       = gpar(fontsize=FONT_SIZE),
        row_title_side     = 'right',
        cluster_rows       = TRUE,
        cluster_columns    = TRUE,
        clustering_distance_rows    = 'manhattan',
        clustering_distance_columns = 'manhattan',
        show_row_dend      = FALSE,
        show_column_dend   = FALSE,
        show_row_names     = TRUE,
        row_names_side     = 'left',
        show_column_names  = TRUE,
        row_names_gp       = gpar(fontsize=8),
        column_names_gp    = gpar(fontsize=7),
        column_names_side  = 'top',
        column_names_rot   = 30,
        width              = unit(0.45 * ncol(pltmat), 'cm'),
        use_raster         = FALSE,
        heatmap_legend_param = list(
            title       = "omic or\npathway state",
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            grid_height = unit(3.5, 'mm'),
            grid_width  = unit(3.5, 'mm'),
            at          = c(1, 0, -1),
            labels      = c('activated', 'normal', 'repressed'),
            border      = 'black'
        )
    ) %>% draw( heatmap_legend_side = 'bottom',
                background = 'transparent' )
}

ppMat4PltNeiStt <- function(cnmat, rnamat, iplmat, nodedt) {
    pats <- intersect(colnames(cnmat), colnames(rnamat)) |>
            intersect(colnames(iplmat)) |> sort()

    type <- NULL
    ents <- intersect(nodedt[ type != 'omic' ]$id, rownames(iplmat))
    pltmat <- rbind(cnmat[, pats], rnamat[, pats], sign(iplmat)[ents, pats])
    rownames(pltmat) <- c('CNA', 'RNA-seq', ents)

    return(pltmat)
}

ppNode4PltNeiStt <- function(fpth, gnname) {
    to <- from <- nei <- entity <- . <- desc <- row_grp <- loc <- NULL
    pthlist <- getNodeEdge(fpth)
    in_nodedt <- pthlist$nodedt
    in_edgedt <- pthlist$edgedt

    omic_edgedt <- copy(OMIC_EDGEDT) |> _[, to := gnname]
    edgedt <- in_edgedt[(from == gnname) | (to %in% gnname)] |>
        rbind(omic_edgedt) |> merge(DASHDT, by='title', all.x=TRUE) |>
        _[, nei := ifelse(from == gnname, to, ifelse(to == gnname, from, NA))]

    up_ents <- edgedt[(to == gnname) & (! from %in% OMIC_NODEDT$entity)]$from
    dn_ents <- edgedt[ from == gnname ]$to

    omic_nodedt <- copy(OMIC_NODEDT)
    outdt <- in_nodedt[entity %in% c(edgedt$from, edgedt$to)] |>
        rbind(omic_nodedt) |> setnames('entity', 'id') |>
        orderNodes(up_ents, dn_ents, gnname) |>
        merge(edgedt[, .(nei, desc)], by.x='id', by.y='nei', all.x=TRUE) |>
        _[, row_grp := ifelse(loc %in% c('omic', gnname), loc,
            paste0(loc, ":\n", desc))]

    return(outdt)
}

orderNodes <- function(nodedt, up_ents, dn_ents, gnname) {
    id <- label <- level <- loc <- NULL
    nodedt[, id := factor(id, 
        levels=unique(c(up_ents, dn_ents, 'CNA', gnname, 'RNA-seq')))] |>
    _[, label := gsub('_', ' ', id)] |>
    _[, level := ifelse(id %in% up_ents, 1, ifelse(id %in% dn_ents, 3, 2))] |>
    _[, loc := ifelse(id %in% up_ents, 'upstream',
        ifelse(id %in% dn_ents, 'downstream',
            ifelse(id %in% c('CNA', 'RNA-seq'), 'omic',
                ifelse(id == gnname, gnname, '')))) ]
}
