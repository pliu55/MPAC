require(magrittr)

main <- function() {
    paradigm_bin = dlParadigmBin()

    testRunRealPrd(paradigm_bin)

    testRunPermPrd(paradigm_bin)
}

testRunPermPrd <- function(paradigm_bin) {
    fpermll = system.file('extdata/TcgaInp/inp_perm.rds', package='MPAC')
    fpth    = system.file('extdata/Pth/tiny_pth.txt',     package='MPAC')

    permll = readRDS(fpermll)
    pat = 'TCGA-CV-7100'
    outdir = tempdir()
    
    runPermPrd(permll, fpth, outdir, paradigm_bin, sampleids=c(pat))

    lapply(permll, testRunPermPrdByIPerm, pat, outdir)
}

testRunPermPrdByIPerm <- function(permlist, pat, outdir) {
    iperm = permlist$iperm
    out_prefix = paste0(outdir, '/p', iperm, '/', pat, '_')
    cmp_relpre = paste0('extdata/runPrd/p', iperm, '/', pat, '_')
    testRunPrdByPrefix('testRunPermPrd', out_prefix, cmp_relpre)
}

testRunRealPrd <- function(paradigm_bin) {
    fcn  = system.file('extdata/TcgaInp/inp_focal.rds',       package='MPAC')
    frna = system.file('extdata/TcgaInp/inp_log10fpkmP1.rds', package='MPAC')
    fpth = system.file('extdata/Pth/tiny_pth.txt',            package='MPAC')

    cnmat  = readRDS(fcn)
    rnamat = readRDS(frna)
    pat = 'TCGA-CV-7100'

    outdir = tempdir()
    runPrd(cnmat, rnamat, fpth, outdir, paradigm_bin, sampleids=c(pat))

    out_prefix = paste0(outdir, '/', pat, '_')
    cmp_relpre = paste0('extdata/runPrd/', pat, '_')

    testRunPrdByPrefix('testRunRealPrd', out_prefix, cmp_relpre)
}

testRunPrdByPrefix <- function(testname, out_prefix, cmp_relpre) {
    suffixes = c('cfg.txt', 'inp_focal.tsv', 'inp_log10fpkmP1.tsv', 'cpt.txt',
        'ipl.txt')
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
