require(data.table)
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
    out = outl[[pat]]
    cmp = cmpl[[pat]] |> upgrade_graph()

    test_that('testSubNtw: vertex', {
        expect_identical( sort(V(out)$name), sort(V(cmp)$name) )
    })

    test_that('testSubNtw: edge', {
        expect_identical(
            as_edgelist(out) |> as.data.table() |> _[order(V1, V2)],
            as_edgelist(cmp) |> as.data.table() |> _[order(V1, V2)]
        )
    })
}

main()
