suppressMessages(require(SummarizedExperiment))

main <- function() {
    testPpCnInp()
    testPpRnaInp()
    testPpRealInp()

    testPpPermInp()
}

testPpRealInp <- function() {
    fcn = system.file('extdata/TcgaInp/focal_tumor.rds', package='MPAC')
    ftumor = system.file('extdata/TcgaInp/log10fpkmP1_tumor.rds',package='MPAC')
    fnorm = system.file('extdata/TcgaInp/log10fpkmP1_normal.rds',package='MPAC')
    fcmp = system.file('extdata/TcgaInp/inp_real.rds', package='MPAC') 
   
    cn_tumor_mat = readRDS(fcn)
    rna_tumor_mat = readRDS(ftumor)
    rna_norm_mat  = readRDS(fnorm)
   
    out_se = ppRealInp(cn_tumor_mat, rna_tumor_mat, rna_norm_mat)
    cmp_se = readRDS(fcmp)

    test_that('ppRealInp', {expect_identical(out_se, cmp_se)})
}

testPpCnInp <- function() {
    fin  = system.file('extdata/TcgaInp/focal_tumor.rds', package='MPAC')
    fcmp = system.file('extdata/TcgaInp/inp_focal.rds',   package='MPAC')
    out_se = readRDS(fin) |> ppCnInp()

    cmp_se = readRDS(fcmp)

    test_that('ppCnInp', {expect_identical(out_se, cmp_se)})
}

testPpRnaInp <- function() {
    ftumor= system.file('extdata/TcgaInp/log10fpkmP1_tumor.rds', package='MPAC')
    fnorm = system.file('extdata/TcgaInp/log10fpkmP1_normal.rds',package='MPAC')
    fcmp  = system.file('extdata/TcgaInp/inp_log10fpkmP1.rds',   package='MPAC')

    tumor_mat = readRDS(ftumor)
    norm_mat  = readRDS(fnorm)
    out_se = ppRnaInp(tumor_mat, norm_mat)

    cmp_se = readRDS(fcmp)

    test_that('ppRnaInp', {expect_identical(out_se, cmp_se)})
}

testPpPermInp <- function() {
    freal = system.file('extdata/TcgaInp/inp_real.rds', package='MPAC') 
    real_se = readRDS(freal)
    n_perms = 3
    set.seed(123456)
    outll = ppPermInp(real_se, n_perms)

    fcmp = system.file('extdata/TcgaInp/inp_perm.rds', package='MPAC')
    cmpll = readRDS(fcmp)

    lapply(seq_len(n_perms), testPpPermInpByIPerm, outll, cmpll)
}

testPpPermInpByIPerm <- function(iperm, outll, cmpll) {
    func = function(x) x$iperm == iperm
    outlist = Filter(func, outll)
    cmplist = Filter(func, cmpll)

    test_that(paste0('ppPermInp: ', iperm), {
        expect_identical(outlist, cmplist )
    })
}

main()
