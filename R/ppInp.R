#' @title Prepare input copy-number (CN) alteration and RNA data to run PARADIGM
#'
#' @inheritParams ppCnInp
#' @inheritParams ppRnaInp
#'
#' @return  A SummarizedExperiment object of CN and RNA state for PARADIGM
#'
#' @export
#'
#' @examples
#'
#' fcn = system.file('extdata/TcgaInp/focal_tumor.rds', package='MPAC')
#' ftumor = system.file('extdata/TcgaInp/log10fpkmP1_tumor.rds', package='MPAC')
#' fnorm = system.file('extdata/TcgaInp/log10fpkmP1_normal.rds', package='MPAC')
#'
#' cn_tumor_mat = readRDS(fcn)
#' rna_tumor_mat = readRDS(ftumor)
#' rna_norm_mat  = readRDS(fnorm)
#'
#' ppRealInp(cn_tumor_mat, rna_tumor_mat, rna_norm_mat)
#'
#' @importFrom SummarizedExperiment assays assays<-
#'
ppRealInp <- function(cn_tumor_mat, rna_tumor_mat, rna_normal_mat, rna_n_sd=2,
    threads=1) {
    cn_se <- ppCnInp(cn_tumor_mat)
    rna_se <- ppRnaInp(rna_tumor_mat, rna_normal_mat, threads=threads)

    out_rownames <- intersect(rownames(cn_se), rownames(rna_se))
    out_colnames <- intersect(colnames(cn_se), colnames(rna_se))
    real_se <- cn_se[out_rownames, out_colnames] |> copy()
    assays(real_se)$RNA_state <-
        assays(rna_se)$RNA_state[out_rownames, out_colnames]

    return(real_se)
}

#' @title Prepare input copy-number (CN) alteration data to run PARADIGM
#'
#' @param cn_tumor_mat  A matrix of tumor CN focal data with rows as genes
#'                      and columns as samples. A value of 0 means normal CN,
#'                      > 0 means amplification, and < 0 means deletion.
#'
#' @return  A SummarizedExperiment object of CN state for PARADIGM
#'
#' @export
#'
#' @examples
#'
#' fcn = system.file('extdata/TcgaInp/focal_tumor.rds', package='MPAC')
#' cn_tumor_mat = readRDS(fcn)
#'
#' ppCnInp(cn_tumor_mat)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#'
ppCnInp <- function(cn_tumor_mat) {
    cn_state_mat <- sign(cn_tumor_mat)
    SummarizedExperiment(assays=list(CN_state=cn_state_mat))
}


#' @title Prepare input RNA data to run PARADIGM
#'
#' @param rna_tumor_mat  A matrix of RNA data from tumor samples with rows as
#'                       genes and columns as samples
#'
#' @param rna_normal_mat  A matrix of RNA data from normal samples with rows
#'                        as genes and columns as samples
#'
#' @param rna_n_sd  Standard deviation range from fitted normal samples to
#'                  define RNA state. Default: 2, i.e. 2*sd
#'
#' @param threads  Number of threads to run in parallel. Default: 1
#'
#' @return A SummarizedExperiment of RNA state for PARADIGM
#'
#' @export
#'
#' @examples
#'
#' ftumor = system.file('extdata/TcgaInp/log10fpkmP1_tumor.rds', package='MPAC')
#' fnorm = system.file('extdata/TcgaInp/log10fpkmP1_normal.rds', package='MPAC')
#' rna_tumor_mat = readRDS(ftumor)
#' rna_norm_mat  = readRDS(fnorm)
#'
#' ppRnaInp(rna_tumor_mat, rna_norm_mat, threads=2)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#' @importFrom BiocParallel  SnowParam bplapply
#'
ppRnaInp <- function(rna_tumor_mat, rna_normal_mat, rna_n_sd=2, threads=1) {
    out_mat <- getBPPARAM(threads) |>
        bplapply(rownames(rna_normal_mat), fitByGn, rna_normal_mat, BPPARAM=_)|>
        rbindlist() |>
        defState(rna_tumor_mat, rna_n_sd)

    SummarizedExperiment(assays=list(RNA_state=out_mat))
}

#' @import data.table
#'
defState <- function(fitdt, tumor_mat, rna_n_sd=2) {
    norm_mean <- norm_sd <- m_m_nsd <- m_p_nsd <- val <- state <- NULL
    fitdt[, `:=`(
        m_m_nsd = norm_mean - rna_n_sd * norm_sd,
        m_p_nsd = norm_mean + rna_n_sd * norm_sd)]

    t(tumor_mat) |> as.data.table(keep.rownames='brc') |>
    melt(id='brc', variable='gnname', value='val') |>
    merge(fitdt, by='gnname', all.x=TRUE) |>
    _[, state := ifelse(val < m_m_nsd, -1, ifelse(val > m_p_nsd, 1, 0))] |>
    dcast(brc ~ gnname, value.var='state') |>
    _[, c('brc', fitdt$gnname), with=FALSE] |>
    as.matrix(rownames='brc') |> t()
}

#' @importFrom fitdistrplus mgedist
#'
fitByGn <- function(gnname, mat) {
    vec <- mat[gnname, ]
    outl <- list(gnname = gnname)
    scaling_dummy <- 100 ## for fitting on very small numbers
    if ( all(vec == vec[1]) ) {
        outl$norm_mean <- vec[1]
        outl$norm_sd   <- 0
    } else {
        fit <- fitdistrplus::mgedist(vec*scaling_dummy, distr='norm')
        outl$norm_mean <- fit$estimate['mean']/scaling_dummy
        outl$norm_sd   <- fit$estimate['sd']/scaling_dummy
    }

    return(outl)
}


#' @title Permute input genomic state data between genes in the same sample
#'
#' @param real_se  A SummarizedExperiment object of CN and RNA states from
#'                 real samples with rows as genes and columns as samples.
#'                 It is the output from `ppRealInp()`.
#'
#' @param n_perms  Number of permutations. Default: 100
#'
#' @inheritParams ppRnaInp
#'
#' @usage ppPermInp(real_se, n_perms=100, threads=1)
#'
#' @return  A list of SummarizedExperiment objects of permuted CN and RNA
#'          states. The metadata `i` in each obbect denotes its permutation
#'          index.
#'
#' @export
#'
#' @examples
#'
#' freal = system.file('extdata/TcgaInp/inp_real.rds', package='MPAC')
#' real_se = readRDS(freal)
#'
#' ppPermInp(real_se, n_perms=3)
#'
ppPermInp <- function(real_se, n_perms=100, threads=1) {
    ngns <- nrow(real_se)
    permlist <- lapply(seq_len(n_perms),
        function(x) sample(ngns, ngns, replace=FALSE))

    getBPPARAM(threads) |>
    bplapply(seq_len(n_perms), ppByIPerm, permlist, real_se, BPPARAM=_)
}

#'
#' @importFrom S4Vectors metadata metadata<-
#'
ppByIPerm <- function(iperm, permlist, real_se) {
    perms <- permlist[[iperm]]

    perm_se <- real_se[perms, ]
    rownames(perm_se) <- rownames(real_se)
    metadata(perm_se)$i <- iperm

    return(perm_se)
}
