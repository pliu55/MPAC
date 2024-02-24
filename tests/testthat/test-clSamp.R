main <- function() {
    testClSamp()
}

testClSamp <- function() {
    fovr = system.file('extdata/clSamp/ovrmat.rds', package='MPAC')
    fcmp = system.file('extdata/clSamp/clmat.rds',  package='MPAC')

    ovrmat = readRDS(fovr)
    outdt = clSamp(ovrmat, n_neighbors=10, n_random_runs=5) %>%
            suppressWarnings()
    
    ## take the most frequent one
    freqmat = outdt[order(-nreps)] %>% as.matrix(rownames='nreps') %>%
        .[1, ] %>% t() %>% t()
    colnames(freqmat) = 'icl'
    cmpmat = readRDS(fcmp)

    test_that('testClSamp', {
        expect_identical(freqmat, cmpmat)
    })
}

main()
