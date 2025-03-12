main <- function() {
    testPltSttKM()
}

testPltSttKM <- function() {
    cdrmat = system.file('extdata/pltSttKM/cdr.rds', package='MPAC') |>readRDS()
    fltmat = system.file('extdata/pltSttKM/ipl.rds', package='MPAC') |>readRDS()

    prts = c('CD247', 'FASLG')

    test_that('testPltSttKM', {
        pltSttKM(cdrmat, fltmat, event='OS', time='OS_days', proteins=prts) |>
        print() |>
        expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
