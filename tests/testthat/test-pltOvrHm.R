main <- function() {
    testPltOvrHm()
}

testPltOvrHm <- function() {
    fovr = system.file('extdata/pltOvrHm/ovr.rds', package='MPAC')
    fcl  = system.file('extdata/pltOvrHm/cl.rds',   package='MPAC')

    ovrmat = readRDS(fovr)
    cldt   = readRDS(fcl)

    test_that('testPltOvrHm', {
        pltOvrHm(ovrmat, cldt) |>
        expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
