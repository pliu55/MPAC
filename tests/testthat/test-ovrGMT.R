main <- function() {
    testOvrGMT()
}

testOvrGMT <- function() {
    fntw   = system.file('extdata/subNtw/subntwl.rds',    package='MPAC')
    fgmt   = system.file('extdata/ovrGMT/fake.gmt',       package='MPAC')
    ffocal = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
    fcmp   = system.file('extdata/ovrGMT/ovr.rds',        package='MPAC')

    subntwl = readRDS(fntw) |> lapply(function(s) upgrade_graph(s))
    omic_gns = readRDS(ffocal) |> rownames()

    outmat = ovrGMT(subntwl, fgmt, omic_gns)
    cmpmat = readRDS(fcmp)

    test_that('testOvrGMT', {
        expect_equal(outmat, cmpmat)
    })
}

main()
