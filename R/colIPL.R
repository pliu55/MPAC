#' @title  Collect Inferred Pathway Levels (IPLs) from PARADIGM runs on real 
#'         data
#'
#' @param indir  Input folder that saves PARADIGM results. It should be set as 
#'               the same as `outdir` as in `runPrd()`.
#'
#' @param sampleids  Sample IDs for which IPLs to be collected. If not provided,
#'                   all files with suffix '_ipl.txt' in `indir` will be 
#'                   collected. Default: NULL.
#'
#' @return  A data.table object with columns of pathway entities and their IPLs.
#'
#' @export
#'
#' @examples
#'
#' indir = system.file('/extdata/runPrd/', package='MPAC')
#'
#' colRealIPL(indir)
#'
colRealIPL <- function(indir, sampleids=NULL) {
    colIPL(indir, sampleids) %>% return()
}

#' @title  Collect Inferred Pathway Levels (IPLs) from PARADIGM runs on permuted
#'         data
#'
#' @inheritParams colRealIPL
#'
#' @param  n_perms  Number of permutations to collect.
#'
#' @return  A data.table object with columns of permutation index, pathway 
#'          entities and their IPLs.
#'
#' @export
#'
#' @examples
#'
#' indir = system.file('/extdata/runPrd/', package='MPAC')
#' n_perms = 3
#'
#' colPermIPL(indir, n_perms)
#'
colPermIPL <- function(indir, n_perms, sampleids=NULL) {
    . <- NULL

    lapply(seq_len(n_perms), function(iperm) {
        ipldt <- paste0(indir, '/p', iperm, '/') %>% colIPL(sampleids)
        brcs <- names(ipldt) %>% setdiff('entity')
        ipldt[, iperm := iperm] %>%
        .[, c('entity', 'iperm', brcs), with=FALSE] %>%
        return()
    }) %>% rbindlist() %>%
    return()
}

colIPL <- function(indir, sampleids) {
    fipls <- NULL
    if ( is.null(sampleids) ) {
        fipls <- list.files(path=indir, pattern="*_ipl.txt", full.names=TRUE,
            recursive=FALSE)
    } else {
        fipls <- paste0(indir, '/', sampleids, '_ipl.txt')
    }

    . <- sampleid <- fipl <- NULL
    fdt <- data.table(fipl = fipls) %>%
        .[, sampleid := basename(fipl) %>% tstrsplit('_ipl.txt') %$% .[[1]] ]

    Map(function(fin, sampleid) {
        readIPL(fin) %>% .[, sampleid := sampleid] %>% return()
    }, fdt$fipl, fdt$sampleid) %>% rbindlist() %>%
    dcast(entity ~ sampleid, value.var='ipl') %>%
    return()
}

readIPL <- function(fin) {
    readLines(fin) %>%
    lapply(function(line) {
        if ((! grepl('^> ',     line, perl=TRUE)) &
            (! grepl('__\\d\t', line, perl=TRUE)) ) {
            words <- strsplit(line, "\t")[[1]]
            list(entity=words[1], ipl=words[2]) %>%
            return()
        }
    }) %>% rbindlist() %>%
    return()
}
