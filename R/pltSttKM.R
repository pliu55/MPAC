#' @title  Plot a Kaplan-Meier curve for samples stratified by given
#'         protein(s)' pathway states
#'
#' @param  cdrmat   A matrix containing survival data with rows as patient
#'                  samples and columns as survival event and time.
#'
#' @param  event  The column name in `cdrmat` to indicate survival event.
#'                Default: 'OS'.
#'
#' @param  time   The column name in `cdrmat` to indicate survival time.
#'                Default: 'OS_days'.
#'
#' @param  proteins  Rowname(s) in `fltmat`. Its/their pathway states will be
#'                   used to stratify patient samples. Default: all proteins
#'                   in `fltmat` will be used.
#'
#' @param strat_func  A function applied on protein(s) pathway states to
#'                    stratify patient samples. Available options: '>0', '<0',
#'                    Default: '>0', i.e., IPL >0 vs. the rest.
#'
#' @inheritParams pltNeiStt
#'
#' @return  A Kaplan-Meier plot
#'
#' @examples
#'
#' cdrmat <- system.file('extdata/pltSttKM/cdr.rds',package='MPAC') |> readRDS()
#' fltmat <- system.file('extdata/pltSttKM/ipl.rds',package='MPAC') |> readRDS()
#'
#' pltSttKM(cdrmat, fltmat, event='OS', time='OS_days',
#'          proteins=c('CD247', 'FASLG'))
#'
#' @export
#'
#'
#' @import ggplot2
#' @importFrom survival  Surv survfit
#' @importFrom survminer ggsurvplot
#' @importFrom stringr   str_wrap
#'
pltSttKM <- function(cdrmat, fltmat, event='OS', time='OS_days', proteins=NULL,
                     strat_func='>0'){

    if ( is.null(proteins) ) {
        proteins = rownames(fltmat)
    }

    pats = intersect(rownames(cdrmat), colnames(fltmat))
    cdrdt = cdrmat[pats, c(event, time)] |>
            as.data.table(keep.rownames='pat') |>
            setnames(c(event, time), c('event', 'time'))

    sttmat = fltmat[proteins, pats, drop=FALSE] |> sign()

    stt_sums = NULL
    if (strat_func == '>0') {
        stt_sums = colSums(sttmat)
    } else if ( strat_func == '<0') {
        stt_sums = colSums(sttmat) * (-1)
    } else {
        paste0("strat_func: ", strat_func, " not recognized\n",
               "It can only be >0 or <0\n") |> stop()
    }
    sel_pats = names(stt_sums[ stt_sums == length(proteins)])
    if ( length(sel_pats) == 0 ) {
        paste0('no sample for ', strat_func, " group\n") |> stop()
    }

    sel_lab = paste0(proteins, strat_func) |> paste(collapse=' & ') |>
              paste0(' (n=', length(sel_pats), ')') |> str_wrap(width=30)
    other_lab = paste0('others (n=', ncol(fltmat) - length(sel_pats), ')')

    grp = pat = NULL
    cdrdt[, grp := ifelse(pat %in% sel_pats, sel_lab, other_lab) |>
                   factor(levels=c(sel_lab, other_lab))]

    srByStt(cdrdt, event, time)
}

srByStt <- function(cdrdt, event, time) {
    frm = Surv(time, event) ~ grp
    MAGENTA = '#CC79A7'
    GREY    = '#999999'
    p = ggsurvplot(
        fit = do.call(survfit, args=list(formula=frm, data=cdrdt)),
        data = cdrdt,
        palette = c(MAGENTA, GREY),
        pval = TRUE,
        pval.size = 3,
        pval.coord = c(0.1, 0.05),
        xlab = time,
        ylab = paste0(event, ' probability'),
        legend = 'top',
        legend.title = '',
        title = '',
        ggtheme = getSrTM()
    )$plot + guides(color=guide_legend(nrow=2))

    return(p)
}

getSrTM <- function() {
    FONT_SIZE = 9

    PAPER_THEME = ggplot2::theme(
        legend.title      = ggplot2::element_text(size=FONT_SIZE),
        legend.text       = ggplot2::element_text(size=FONT_SIZE),
        axis.title        = ggplot2::element_text(size=FONT_SIZE),
        axis.text.x       = ggplot2::element_text(size=FONT_SIZE),
        axis.text.y       = ggplot2::element_text(size=FONT_SIZE),
        plot.title        = ggplot2::element_text(size=FONT_SIZE, hjust=0.5),
        legend.margin     = ggplot2::margin(b=-0.3, l=-0.3, unit='cm'),
        legend.key.size   = ggplot2::unit(0.35, 'cm'),
        panel.background  = ggplot2::element_rect(fill='transparent', color=NA),
        plot.background   = ggplot2::element_rect(fill='transparent', color=NA),
        legend.background = ggplot2::element_rect(fill='transparent', color=NA)
    )

    tm = theme_classic() + theme(aspect.ratio = 1) + PAPER_THEME

    return(tm)
}
