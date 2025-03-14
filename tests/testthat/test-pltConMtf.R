main <- function() {
    testPltConMtf()
}

testPltConMtf <- function() {
    grphl = system.file('extdata/pltMtfPrtIPL/grphl.rds',package='MPAC') |>
            readRDS()

    prots = c('CD3G', 'CD86')

    test_that('testPltConMtf', {
        pltConMtf(grphl, prots) |> print() |> expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
