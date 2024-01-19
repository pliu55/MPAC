require(magrittr)

main <- function() {
    testClSamp()
}

testClSamp <- function() {
    fovr = system.file('extdata/clSamp/ovrmat.rds', package='MPAC')
    fcmp = system.file('extdata/clSamp/clmat.rds',  package='MPAC')

    ovrmat = readRDS(fovr)
    set.seed(123456)
    outmat = clSamp(ovrmat)

    cmpmat = readRDS(fcmp)

    test_that('testClSamp', {
        expect_identical(outmat, cmpmat)
    })
}

main()
