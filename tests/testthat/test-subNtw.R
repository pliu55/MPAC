suppressMessages(require(igraph))

main <- function() {
    testSubNtw()
}

testSubNtw <- function() {
    fflt = system.file('extdata/fltByPerm/flt_real.rds', package='MPAC')
    fpth = system.file('extdata/Pth/tiny_pth.txt',       package='MPAC')
    fgmt = system.file('extdata/ovrGMT/fake.gmt',        package='MPAC')

    fltdt = readRDS(fflt)
    outl = subNtw(fltdt, fpth, fgmt, min_n_gmt_gns=1)
    cmpl = system.file('extdata/subNtw/subntwl.rds', package='MPAC') |> 
           readRDS()

    lapply(names(outl), cmpSubNtwByPat, outl, cmpl)
}

cmpSubNtwByPat <- function(pat, outl, cmpl) {
    test_that('testSubNtw', {
        identical_graphs(outl[[pat]], cmpl[[pat]]) |> 
        expect_identical(TRUE)
    })
}

main()
