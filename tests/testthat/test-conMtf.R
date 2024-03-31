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
    test_that(paste0('testConMtf: ', ind), {
        igraph::identical_graphs(outl[[ind]], cmpl[[ind]]) |> 
        expect_identical(TRUE)
    })
}

main()
