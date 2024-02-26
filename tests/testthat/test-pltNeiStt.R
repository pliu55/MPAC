main <- function() {
    testPltNeiStt()
}

testPltNeiStt <- function() {
    fpth = system.file('extdata/Pth/tiny_pth.txt', package='MPAC')

    fcn  = system.file('extdata/pltNeiStt/inp_focal.rds',       package='MPAC')
    frna = system.file('extdata/pltNeiStt/inp_log10fpkmP1.rds', package='MPAC')
    fflt = system.file('extdata/pltNeiStt/fltmat.rds',          package='MPAC')

    cn_state_mat  = readRDS(fcn)
    rna_state_mat = readRDS(frna)
    fltmat = readRDS(fflt)
    protein = 'CD86'

    test_that('testPltNeiStt', {
        pltNeiStt(cn_state_mat, rna_state_mat, fltmat, fpth, protein) |>
        expect_no_error()
    })

    fpdf = 'Rplots.pdf'
    if ( file.exists(fpdf) ) file.remove(fpdf)
}

main()
