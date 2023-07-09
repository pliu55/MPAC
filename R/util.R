#' @import magrittr
#' @importFrom BiocParallel SnowParam
#'
getBPPARAM <- function(threads) {
    . = NULL
    ifelse(Sys.info()[["sysname"]]=="Windows", 'SOCK', 'FORK') %>%
    SnowParam(workers=threads, type=.) %>%
    return()
}

getOS <- function() {
    os = toupper(.Platform$OS.type)
    sysinf = Sys.info()
    if ( !is.null(sysinf) ) {
        os = sysinf['sysname']
        if ( os == 'Darwin' ) {
            os = "OSX"
        } else if ( os == 'Linux' ) {
            os = "LINUX"
        } else {
            os = "WINDOWS"
        }
    } else {
        if ( grepl("^darwin", R.version$os, perl=TRUE) ) {
            os = "OSX"
        } else if ( grepl("linux-gnu", R.version$os, perl=TRUE) ) {
            os = "LINUX"
        } else {
            os = "WINDOWS"
        }
    }

    return(os)
}

#' @title  Download PARADIGM binary and return its location
#'
#' @return full path of downloaded PARADIGM binary
#'
#' @importFrom utils download.file
#'
dlParadigmBin <- function() {
    os = getOS()
    
    ## to get $url, right-click the 'raw' button on Github file page
    url = 'https://github.com/sng87/paradigm-scripts/raw/master/public/exe/'
    if ( os == 'LINUX' ) {
        url = paste0(url, 'LINUX/paradigm')
    } else if ( os == 'OSX' ) {
        url = paste0(url, 'MACOS/paradigm')
    }
    flocal = paste0(tempdir(), '/paradigm')

    ## wget --no-check-certificate
    ## curl -LJ $url
    if ( ! file.exists(flocal) ) {
        download.file(url, flocal, method='curl', extra='-LJ')
    }
    Sys.chmod(flocal, mode='755')

    return(flocal)
}

## file.exists cannot give T/F for character(0)
fileExists <- function(f) {
    is_existed = ifelse( identical(f, character(0)), FALSE,
        ifelse( ! file.exists(f), FALSE, TRUE))
    return(is_existed)
}


getNodeEdge <- function(fpth) {
    . = entity = type = from = to = title = NULL

    wordslist = readLines(fpth) %>% strsplit("\t")

    nodedt = lapply(wordslist, function(words) {
        if (length(words) == 2) as.list(words) %>% return()
    }) %>% rbindlist() %>% setnames( c('type', 'entity') ) %>% 
    .[, .(entity, type)]

    edgedt = lapply(wordslist, function(words) {
        if (length(words) == 3) as.list(words) %>% return()
    }) %>% rbindlist() %>% setnames( c('from', 'to', 'title') ) %>%
    .[, .(from, to, title)]

    list(nodedt=nodedt, edgedt=edgedt) %>% return()
}
