main <- function() {
    testColRealIPL()

    testColPermIPL()
}

testColRealIPL <- function() {
    indir = system.file('/extdata/runPrd/', package='MPAC')
    pat = 'TCGA-CV-7100'
    
    out_alldt = colRealIPL(indir, sampleids=NULL)
    out_patdt = colRealIPL(indir, sampleids=c(pat))

    cmpdt = system.file('/extdata/colIPL/real.rds', package='MPAC') |> readRDS()
    testColIPL('testColRealIPLAll', out_alldt, cmpdt)
    testColIPL('testColRealIPLPat', out_patdt, cmpdt)
}

testColPermIPL <- function() {
    indir = system.file('/extdata/runPrd/', package='MPAC')
    n_perms = 3
    pat = 'TCGA-CV-7100'

    out_alldt = colPermIPL(indir, n_perms, sampleids=NULL)
    out_patdt = colPermIPL(indir, n_perms, sampleids=c(pat))

    cmpdt = system.file('/extdata/colIPL/perm.rds', package='MPAC') |> readRDS()
    testColIPL('testColRealIPLAll', out_alldt, cmpdt)
    testColIPL('testColRealIPLPat', out_patdt, cmpdt)
}

testColIPL <- function(tag, outdt, cmpdt) {
    test_that(tag, {
        expect_identical(outdt, cmpdt)
    })
}

main()
