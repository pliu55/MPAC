#' @importFrom BiocParallel SnowParam MulticoreParam
#'
getBPPARAM <- function(threads) {
    bp <- NULL
    if ( getOS() == 'Windows' ) {
        bp <- SnowParam(workers=threads, type='SOCK')
    } else {
        bp <- MulticoreParam(workers=threads)
    }
    return(bp)
}

getOS <- function() {
    os <- toupper(.Platform$OS.type)
    sysinf <- Sys.info()
    if ( !is.null(sysinf) ) {
        os <- sysinf['sysname']
        if ( os == 'Darwin' ) {
            os <- "OSX"
        } else if ( os == 'Linux' ) {
            os <- "LINUX"
        } else {
            os <- "WINDOWS"
        }
    } else {
        if ( grepl("^darwin", R.version$os, perl=TRUE) ) {
            os <- "OSX"
        } else if ( grepl("linux-gnu", R.version$os, perl=TRUE) ) {
            os <- "LINUX"
        } else {
            os <- "WINDOWS"
        }
    }

    return(os)
}

getNodeEdge <- function(fpth) {
    . <- entity <- type <- from <- to <- title <- NULL

    wordslist <- readLines(fpth) |> strsplit("\t")

    nodedt <- lapply(wordslist, function(words) {
        if (length(words) == 2) as.list(words)
    }) |> rbindlist() |> setnames( c('type', 'entity') ) |>
    _[, .(entity, type)]

    edgedt <- lapply(wordslist, function(words) {
        if (length(words) == 3) as.list(words)
    }) |> rbindlist() |> setnames( c('from', 'to', 'title') ) |>
    _[, .(from, to, title)]

    list(nodedt=nodedt, edgedt=edgedt)
}
