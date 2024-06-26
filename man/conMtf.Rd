% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conMtf.R
\name{conMtf}
\alias{conMtf}
\title{Find consensus pathway motifs from a list of pathways}
\usage{
conMtf(subntwl, omic_genes = NULL, min_mtf_n_nodes = 5)
}
\arguments{
\item{subntwl}{A list of igraph objects representing input pathways from
different samples. It is the output from \code{subNtw()}}

\item{omic_genes}{A vector of gene symbols to narrow down
over-representation calculation to only those with input
genomic data. If not provided, all genes in the GMT file
will be considered. Default: NULL.}

\item{min_mtf_n_nodes}{Number of minimum nodes in a motif. Default: 5}
}
\value{
A list of igraph objects representing consensus pathway motifs
}
\description{
Find consensus pathway motifs from a list of pathways
}
\examples{

fsubntwl = system.file('extdata/conMtf/subntwl.rds', package='MPAC')
subntwl = readRDS(fsubntwl)

fomic_gns = system.file('extdata/TcgaInp/inp_focal.rds', package='MPAC')
omic_gns = rownames(readRDS(fomic_gns))

conMtf(subntwl, omic_gns, min_mtf_n_nodes=50)

}
