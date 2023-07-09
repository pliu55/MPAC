require(magrittr)

main <- function() {
    testPpCnInp()
    testPpRnaInp()

    testPpPermInp()
}

testPpCnInp <- function() {
    fin  = system.file('extdata/TcgaInp/focal_tumor.rds', package='MPAC')
    fcmp = system.file('extdata/TcgaInp/inp_focal.rds',   package='MPAC')
    inmat = readRDS(fin)

    outmat = ppCnInp(inmat)

    cmpmat = readRDS(fcmp)

    test_that('ppCnInp', {expect_identical(outmat, cmpmat)})
}

testPpRnaInp <- function() {
    ftumor= system.file('extdata/TcgaInp/log10fpkmP1_tumor.rds', package='MPAC')
    fnorm = system.file('extdata/TcgaInp/log10fpkmP1_normal.rds',package='MPAC')
    fcmp  = system.file('extdata/TcgaInp/inp_log10fpkmP1.rds',   package='MPAC')

    tumor_mat = readRDS(ftumor)
    norm_mat  = readRDS(fnorm)
    outmat = ppRnaInp(tumor_mat, norm_mat)

    cmpmat = readRDS(fcmp)

    test_that('ppRnaInp', {expect_identical(outmat, cmpmat)})
}

testPpPermInp <- function() {
    fin_cn  = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
    fin_rna = system.file('extdata/TcgaInp/inp_log10fpkmP1.rds',
        package='MPAC')
   
    real_cn_mat  = readRDS(fin_cn)
    real_rna_mat = readRDS(fin_rna)
    n_perms = 3
    set.seed(123456)
    outll = ppPermInp(real_cn_mat, real_rna_mat, n_perms)

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
