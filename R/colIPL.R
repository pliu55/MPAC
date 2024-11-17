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
#' @param file_tag  A string of output file name tag. Default: NULL
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
colRealIPL <- function(indir, sampleids=NULL, file_tag=NULL) {
    colIPL(indir, sampleids, file_tag)
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
    lapply(seq_len(n_perms), function(iperm) {
        ipldt <- paste0(indir, '/p', iperm, '/') |> colIPL(sampleids, iperm)
        brcs <- names(ipldt) |> setdiff('entity')
        ipldt[, iperm := iperm] |>
        _[, c('entity', 'iperm', brcs), with=FALSE]
    }) |> rbindlist()
}

colIPL <- function(indir, sampleids, file_tag=NULL) {
    fipls <- NULL
    suffix <- ifelse(is.null(file_tag), '_ipl.txt',
        paste0('_', file_tag, '_ipl.txt'))

    if ( is.null(sampleids) ) {
        fipls <- list.files(path=indir, pattern="*_ipl.txt", full.names=TRUE,
            recursive=FALSE)
    } else {
        fipls <- paste0(indir, '/', sampleids, suffix)
    }

    sampleid <- fipl <- NULL
    fdt <- data.table(fipl = fipls) |>
        _[, sampleid := basename(fipl) |> tstrsplit(suffix) |> _[[1]] ]

    Map(function(fin, sampleid) {
        readIPL(fin) |> _[, sampleid := sampleid]
    }, fdt$fipl, fdt$sampleid) |> rbindlist() |>
    dcast(entity ~ sampleid, value.var='ipl')
}

readIPL <- function(fin) {
    readLines(fin) |>
    lapply(function(line) {
        if ((! grepl('^> ',     line, perl=TRUE)) &
            (! grepl('__\\d\t', line, perl=TRUE)) ) {
            words <- strsplit(line, "\t")[[1]]
            list(entity=words[1], ipl=as.numeric(words[2]))
        }
    }) |> rbindlist()
}
