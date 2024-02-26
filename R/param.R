FONT_SIZE = 11
MAGENTA = '#CC79A7'
CYAN    = '#56B4E9'

DASHDT <- data.table( rbind(
    c( '-a|',        'protein repression'         ),
    c( '-t|',        'transcriptional repression' ),
    c( 'component>', 'is a component of'          ),
    c( 'member>',    'is a member of'             ),
    c( '-a>',        'protein activation'         ),
    c( '-t>',        'transcriptional activation' ),
    c( '-omic>',     'omic'                       )
)) |> setnames( c('title', 'desc') )

OMIC_NODEDT <- data.table( rbind(
    c( 'CNA',     'omic' ),
    c( 'RNA-seq', 'omic' )
)) |> setnames( c('entity', 'type') )

OMIC_EDGEDT <- data.table( rbind(
    c( 'CNA',     '-omic>' ),
    c( 'RNA-seq', '-omic>' )
)) |> setnames( c('from', 'title') )
