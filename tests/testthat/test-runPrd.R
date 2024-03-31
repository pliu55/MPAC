main <- function() {
    # paradigm_bin = dlParadigmBin()

    testRunRealPrd()

    testRunPermPrd()
}

testRunPermPrd <- function() {
    fpermll = system.file('extdata/TcgaInp/inp_perm.rds', package='MPAC')
    fpth    = system.file('extdata/Pth/tiny_pth.txt',     package='MPAC')

    permll = readRDS(fpermll)
    pat = 'TCGA-CV-7100'
    outdir = tempdir()
    
    runPermPrd(permll, fpth, outdir, PARADIGM_bin=NULL, sampleids=c(pat))

    lapply(permll, testRunPermPrdByIPerm, pat, outdir)
}

testRunPermPrdByIPerm <- function(perm_se, pat, outdir) {
    iperm = S4Vectors::metadata(perm_se)$i
    out_prefix = paste0(outdir, '/p', iperm, '/', pat, '_')
    cmp_relpre = paste0('extdata/runPrd/p', iperm, '/', pat, '_')
    testRunPrdByPrefix('testRunPermPrd', out_prefix, cmp_relpre)
}

testRunRealPrd <- function() {
    freal = system.file('extdata/TcgaInp/inp_real.rds', package='MPAC')
    fpth  = system.file('extdata/Pth/tiny_pth.txt',     package='MPAC')

    real_se = readRDS(freal)
    pat = 'TCGA-CV-7100'

    outdir = tempdir()
    runPrd(real_se, fpth, outdir, PARADIGM_bin=NULL, sampleids=c(pat))

    out_prefix = paste0(outdir, '/', pat, '_')
    cmp_relpre = paste0('extdata/runPrd/', pat, '_')

    testRunPrdByPrefix('testRunRealPrd', out_prefix, cmp_relpre)
}

testRunPrdByPrefix <- function(testname, out_prefix, cmp_relpre) {
    # suffixes = c('cfg.txt', 'inp_focal.tsv', 'inp_log10fpkmP1.tsv', 'cpt.txt',
    #     'ipl.txt')
    suffixes = c('cfg.txt', 'inp_focal.tsv', 'inp_log10fpkmP1.tsv')
    lapply(suffixes, testFileContent, testname, out_prefix, cmp_relpre)
}

testFileContent <- function(suffix, testname, out_prefix, cmp_relpre) {
    fout = paste0(out_prefix, suffix)
    fcmp = paste0(cmp_relpre, suffix) %>% system.file(package='MPAC')
    outlines = readLines(fout)
    cmplines = readLines(fcmp)

    test_that(paste0(testname, ':', suffix), {
        expect_identical(outlines, cmplines)
    })
}

main()
