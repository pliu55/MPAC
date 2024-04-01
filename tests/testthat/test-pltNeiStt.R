main <- function() {
    testPltNeiStt()
}

testPltNeiStt <- function() {
    fpth = system.file('extdata/Pth/tiny_pth.txt', package='MPAC')

    freal = system.file('extdata/pltNeiStt/inp_real.rds', package='MPAC')
    fflt  = system.file('extdata/pltNeiStt/fltmat.rds',   package='MPAC')

    real_se = readRDS(freal)
    fltmat = readRDS(fflt)
    protein = 'CD86'

    test_that('testPltNeiStt', {
        pltNeiStt(real_se, fltmat, fpth, protein) |>
        expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
