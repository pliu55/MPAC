#' @title  Plot a heatmap of over-represented gene sets for clsutered samples
#'
#' @param  ovrmat   A matrix containing over-representation adjusted P with
#'                  rows as gene set names and columns as sample IDs. It is the
#'                  output of the `ovrGMT()` function.
#'
#' @param  cldt     A data table with each row representing one clustering
#'                  result, and the first column denotes the number of
#'                  occurrences of a clustering result and the rest of columns
#'                  indicating each sample's cluster index. It is the output
#'                  of the `clSamp()` function. Only the most frequent
#'                  clustering result will be used to plot.
#'
#' @param  min_frc  A minimum fraction of samples in a cluster that have a gene
#'                  set significantly over-represented (adjusted P < 0.05).
#'                  This is used to select gene sets to plot. Default: 0.8
#'
#' @return  A heatmap with rows as over-represented gene sets and columns as
#'          samples splited by clusters.
#'
#' @examples
#'
#' ovrmat <- system.file('extdata/pltOvrHm/ovr.rds',package='MPAC') |> readRDS()
#' cldt   <- system.file('extdata/pltOvrHm/cl.rds', package='MPAC') |> readRDS()
#'
#' pltOvrHm(ovrmat, cldt)
#'
#' @export
#'
pltOvrHm <- function(ovrmat, cldt, min_frc=0.8) {
    nreps <- nsamps <- icl <- padj <- is_signif <- . <- goname <- frc <- NULL

    setnames(cldt, 1, 'nreps')
    srt_cldt <- cldt[order(-nreps)][1, ] |>
        melt(id='nreps', variable='brc', value='icl') |>
        _[, nsamps := .N, by=icl] |> _[, nreps := NULL]

    gonames <- as.data.table(ovrmat, keep.rownames='goname') |>
        melt(id='goname', variable='brc', value='padj') |>
        _[, padj := ifelse(is.na(padj), 1.0, padj)] |>
        _[, is_signif := ifelse(padj < 0.05, TRUE, FALSE)] |>
        merge(srt_cldt, by='brc', all.x=TRUE) |>
        _[, .(frc = sum(is_signif)/nsamps), by=.(goname, icl)] |>
        _[ frc >= min_frc ]$goname |> unique()

    pltmat <- ovrmat[ gonames, srt_cldt$brc ] |> log10()
    rownames(pltmat) <- gsub('_', ' ', gonames)

    makeOvrHm(pltmat, srt_cldt, min_frc)
}

#' @import ComplexHeatmap
#' @importFrom grid  gpar unit
#' @importFrom stringr  str_wrap
#' @importFrom scales   percent
#' @importFrom circlize colorRamp2
#' @importFrom viridis  cividis
#'
makeOvrHm <- function(pltmat, cldt, min_frc) {
    icl_lab <- icl <- nsamps <- NULL

    OVR_CLRS <- colorRamp2(seq(-4, 0, 0.1), rev(cividis(41)))
    FONT_SIZE <- 9
    row_title <- paste0( nrow(pltmat),
        ' gene sets significantly over-represented in >= ', percent(min_frc),
        ' samples in a group') |> str_wrap(width=45)
    cldt[, icl_lab := paste0('c', icl, "\nn=", nsamps)]

    Heatmap( pltmat,
        col = OVR_CLRS,
        na_col = 'black',
        border            = 'black',
        column_split      = cldt$icl_lab,
        column_title_gp   = gpar(fontsize=FONT_SIZE),
        cluster_column_slices = FALSE,
        show_row_names    = TRUE,
        row_names_gp      = gpar(fontsize=7),
        show_column_names = FALSE,
        column_names_side = 'top',
        column_names_gp   = gpar(fontsize=FONT_SIZE),
        show_row_dend     = FALSE,
        show_column_dend  = FALSE,
        row_title         = row_title,
        row_title_gp      = gpar(fontsize=FONT_SIZE),
        cluster_columns   = TRUE,
        cluster_rows      = TRUE,
        use_raster        = ifelse( (nrow(pltmat) > 20) | (ncol(pltmat) > 20),
            TRUE, FALSE),
        raster_quality    = 5,
        heatmap_legend_param = list(
            title       = expression(paste('log'[10], '(adjusted P)')),
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            legend_width = unit(0.8, 'inch'),
            title_position = 'lefttop',
            direction   = 'horizontal'
        )
    ) |> draw(
        background      = 'transparent',
        heatmap_legend_side = 'bottom',
        column_title    = paste0(ncol(pltmat), ' samples in total'),
        column_title_gp = gpar(fontsize=FONT_SIZE),
        padding = unit(c(2, 2, 2, 35), 'mm')
    )
}
