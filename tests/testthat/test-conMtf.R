require(data.table)
suppressMessages(require(igraph))

main <- function() {
    testConMtf()
}

testConMtf <- function() {
    fin = system.file('extdata/conMtf/subntwl.rds', package='MPAC')
    subntwl = readRDS(fin)

    ffocal = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
    omic_gns = readRDS(ffocal) |> rownames()

    outl = conMtf(subntwl, omic_gns, min_mtf_n_nodes=50)

    cmpl = system.file('extdata/conMtf/con_grph.rds', package='MPAC') |>
           readRDS()

    lapply(seq_len(length(outl)), testConMtfByIndi, outl, cmpl)
}

testConMtfByIndi <- function(ind, outl, cmpl) {
    out = outl[[ind]]
    cmp = cmpl[[ind]]

    test_that('testConMtf: vertex', {
        expect_identical( sort(V(out)$name), sort(V(cmp)$name) )
    })

    test_that('testConMtf: edge', {
        expect_identical(
            as_edgelist(out) |> as.data.table() |> _[order(V1, V2)],
            as_edgelist(cmp) |> as.data.table() |> _[order(V1, V2)]
        )
    })
}

main()
