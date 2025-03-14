main <- function() {
    testPltMtfPrtIPL()
}

testPltMtfPrtIPL <- function() {
    fltmat = system.file('extdata/pltSttKM/ipl.rds',package='MPAC') |> readRDS()
    cldt = system.file('extdata/pltMtfPrtIPL/cl.rds',package='MPAC')|> readRDS()
    grphl = system.file('extdata/pltMtfPrtIPL/grphl.rds',package='MPAC') |>
            readRDS()

    test_that('testPltMtfPrtIPL', {
        pltMtfPrtIPL(fltmat, cldt, grphl, proteins=c('CD247', 'FASLG')) |>
        print() |>
        expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
