#' @title Prepare input copy-number (CN) alteration data to run PARADIGM
#'
#' @param cn_tumor_mat  A matrix of tumor CN focal data with rows as genes 
#'                      and columns as samples 
#'
#' @return  A matrix of CN state for PARADIGM
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
#' @importFrom magrittr %>%
#'
ppCnInp <- function(cn_tumor_mat) {
    sign(cn_tumor_mat) %>% return()
}


#' @title Prepare input RNA data to run PARADIGM
#'
#' @param rna_tumor_mat  A matrix of RNA data from tumor samples with rows as
#'                       genes and columns as samples
#'
#' @param rna_normal_mat  A matrix of RNA data from normal samples with rows
#'                        as genes and columns as samples
#'
#' @param threads  Number of threads to run in parallel. Default: 1
#'
#' @return A matrix of RNA state for PARADIGM
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
#' @import data.table
#' @import magrittr
#' @importFrom BiocParallel  SnowParam bplapply
#'
ppRnaInp <- function(rna_tumor_mat, rna_normal_mat, threads=1) {
    . = NULL

    getBPPARAM(threads) %>%
    bplapply(rownames(rna_normal_mat), fitByGn, rna_normal_mat, BPPARAM=.) %>% 
    rbindlist() %>%
    defState(rna_tumor_mat) %>% 
    return()
}

defState <- function(fitdt, tumor_mat) {
    norm_mean = norm_sd = m_m_2sd = m_p_2sd = val = state = . = NULL
    fitdt[, `:=`(
        m_m_2sd = norm_mean - 2 * norm_sd,
        m_p_2sd = norm_mean + 2 * norm_sd)]

    t(tumor_mat) %>% as.data.table(keep.rownames='brc') %>%
    melt(id='brc', variable='gnname', value='val') %>%
    merge(fitdt, by='gnname', all.x=TRUE) %>%
    .[, state := ifelse(val < m_m_2sd, -1, ifelse(val > m_p_2sd, 1, 0))] %>%
    dcast(brc ~ gnname, value.var='state') %>%
    .[, c('brc', fitdt$gnname), with=FALSE] %>%
    as.matrix(rownames='brc') %>% t() %>%
    return()
}

#' @importFrom fitdistrplus mgedist
#'
fitByGn <- function(gnname, mat) {
    vec = mat[gnname, ]
    outl = list(gnname = gnname)
    scaling_dummy = 10 ## for fitting on very small numbers
    if ( all(vec == vec[1]) ) {
        outl$norm_mean = vec[1]
        outl$norm_sd   = 0
    } else {
        fit = fitdistrplus::mgedist(vec*scaling_dummy, distr='norm')
        outl$norm_mean = fit$estimate['mean']/scaling_dummy
        outl$norm_sd   = fit$estimate['sd']/scaling_dummy
    }

    return(outl)
}


#' @title Permute input genomic state data between genes in the same sample
#'
#' @param real_cn_mat  A matrix of CNA states from real samples with rows as
#'                     genes and columns as samples. It is the output from 
#'                     `ppCnInp()`.
#'
#' @param real_rna_mat  A matrix of RNA state from real samples with rows as 
#'                      genes and columns as samples. It is the output from 
#'                      `ppRnaInp()`.
#'
#' @param n_perms  Number of permutations. Default: 3
#'
#' @inheritParams ppRnaInp
#'
#' @usage ppPermInp(real_cn_mat, real_rna_mat, n_perms=3, threads=1)
#'
#' @return  A list of list of matrix. The top level is by permutation index and
#'          the next level saves permutation index, CNA matrix, RNA matrix. A 
#'          matrix contains permuted CNA or RNA states as the input for 
#'          PARADIGM.
#'
#' @export
#' 
#' @examples
#'
#' fcn  = system.file('extdata/TcgaInp/inp_focal.rds',       package='MPAC')
#' frna = system.file('extdata/TcgaInp/inp_log10fpkmP1.rds', package='MPAC')
#' real_cn_mat  = readRDS(fcn)
#' real_rna_mat = readRDS(frna)
#'
#' ppPermInp(real_cn_mat, real_rna_mat, n_perms=3)
#'
ppPermInp <- function(real_cn_mat, real_rna_mat, n_perms=3, threads=1) {
    . = NULL

    in_cnvmat = t(real_cn_mat)
    in_rnamat = t(real_rna_mat)

    ngns = ncol(in_cnvmat)

    permlist = lapply(seq_len(n_perms), 
        function(x) sample(ngns, ngns, replace=FALSE))

    getBPPARAM(threads) %>%
    bplapply(seq_len(n_perms), ppByIPerm, permlist, in_cnvmat, in_rnamat,
        BPPARAM=.) %>%
    return()
}

#' @importFrom magrittr %>%
#'
ppByIPerm <- function(iperm, permlist, in_cnvmat, in_rnamat) {
    perms = permlist[[iperm]]

    cnvmat = in_cnvmat[, perms]
    rnamat = in_rnamat[, perms]

    colnames(cnvmat) = colnames(in_cnvmat)
    colnames(rnamat) = colnames(in_rnamat)
    list(iperm=iperm, CN=t(cnvmat), RNA=t(rnamat)) %>% return()
}
