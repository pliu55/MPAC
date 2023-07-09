require(magrittr)

main <- function() {
    testFltByPerm()
}

testFltByPerm <- function() {
    realdt = system.file('extdata/fltByPerm/real.rds', package='MPAC') %>% 
             readRDS()
    permdt = system.file('extdata/fltByPerm/perm.rds', package='MPAC') %>%
             readRDS()

    outmat = fltByPerm(realdt, permdt)

    cmpmat = system.file('extdata/fltByPerm/flt_real.rds', package='MPAC') %>%
             readRDS()

    test_that('testFltByPerm', {
        expect_identical(outmat, cmpmat)
    })
}

main()
