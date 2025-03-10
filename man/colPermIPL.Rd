% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colIPL.R
\name{colPermIPL}
\alias{colPermIPL}
\title{Collect Inferred Pathway Levels (IPLs) from PARADIGM runs on permuted
data}
\usage{
colPermIPL(indir, n_perms, sampleids = NULL, threads = 1)
}
\arguments{
\item{indir}{Input folder that saves PARADIGM results. It should be set as
the same as \code{outdir} as in \code{runPrd()}.}

\item{n_perms}{Number of permutations to collect.}

\item{sampleids}{Sample IDs for which IPLs to be collected. If not provided,
all files with suffix '_ipl.txt' in \code{indir} will be
collected. Default: NULL.}

\item{threads}{Number of threads to run in parallel. Default: 1}
}
\value{
A data.table object with columns of permutation index, pathway
entities and their IPLs.
}
\description{
Collect Inferred Pathway Levels (IPLs) from PARADIGM runs on permuted
data
}
\examples{

indir = system.file('/extdata/runPrd/', package='MPAC')
n_perms = 3

colPermIPL(indir, n_perms)

}
