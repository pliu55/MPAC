#' @title  Filter IPLs from real data by distribution from permuted data
#'
#' @param realdt  A data.table object containing entities and their IPLs from
#'                real data. It is the output from `colRealIPL()`.
#'
#' @param permdt  A data.table object containing permutation index, entities
#'                and their IPLs from permuted data. It is the output from
#'                `colPermIPL()`.
#'
#' @param threads  Number of threads to run in parallel. Default: 1
#'
#' @return A matrix of filtered IPLs with rows as entities and columns as
#'         samples. Entities with IPLs observed by chance are set to NA.
#'
#' @export
#'
#' @examples
#'
#' freal = system.file('extdata/fltByPerm/real.rds', package='MPAC')
#' fperm = system.file('extdata/fltByPerm/perm.rds', package='MPAC')
#' realdt = readRDS(freal)
#' permdt = readRDS(fperm)
#'
#' fltByPerm(realdt, permdt)
#'
#' @importFrom BiocParallel  SnowParam bplapply
#'
fltByPerm <- function(realdt, permdt, threads=1) {
    pats <- intersect(names(realdt), names(permdt)) |> setdiff('entity')

    bp <- getBPPARAM(threads)
    bplapply(pats, fltByPat, realdt, permdt, BPPARAM=bp) |> rbindlist() |>
    dcast(entity ~ pat, value.var='flt_real_ipl') |>
    as.matrix(rownames='entity')
}

#' @importFrom stats  mad median
fltByPat <- function(pat, in_realdt, in_permdt) {
    realdt <- in_realdt[, c('entity', pat), with=FALSE] |>
        setnames(pat, 'real_ipl')

    perm_mat <- in_permdt[, c('entity', 'iperm', pat), with=FALSE] |>
        dcast(iperm ~ entity, value.var=pat) |> as.matrix(rownames='iperm')

    median_vec <- apply(perm_mat, 2, median)
    mad_vec    <- apply(perm_mat, 2, mad)

    perm_median <- perm_mad <- is_real <- med_m_3mad <- med_p_3mad <- NULL
    . <- real_ipl <- flt_real_ipl <- entity <- NULL

    maddt <- data.table(
        entity      = names(mad_vec),
        perm_median = median_vec,
        perm_mad    = mad_vec ) |>
    merge(realdt, y=_, by='entity', all=TRUE) |>
    _[, `:=`(
        med_m_3mad = perm_median - 3 * perm_mad,
        med_p_3mad = perm_median + 3 * perm_mad )] |>
    _[, is_real := ifelse( (real_ipl < med_m_3mad) |
                            (real_ipl > med_p_3mad), TRUE, FALSE )]

    maddt[, flt_real_ipl := ifelse(is_real == TRUE, real_ipl, NA)] |>
    _[, pat := pat] |>
    _[, .(pat, entity, flt_real_ipl)]
}
