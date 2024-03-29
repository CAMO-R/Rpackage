% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdsGlobal.R
\name{mdsGlobal}
\alias{mdsGlobal}
\title{Genome-wide MDS plot for multiple pairs}
\usage{
mdsGlobal(acs, model.name, sep = "_", file)
}
\arguments{
\item{acs:}{pairwise c-scores or d-scores matrix (c-scores is recommended)}

\item{model.name:}{a vector of dataset names.}

\item{sep:}{a character string to separate the data terms. Please avoid characters appeared in model.name.}

\item{file:}{the output file name.}
}
\value{
A MDS plot of c-scores or d-scores for all study pairs to visualize global distance between each data pair.
}
\description{
The \code{mdsGlobal} is function to generate : multidimensional scaling plot for genome-wide c-scores or d-scores.
}
\examples{
\dontrun{
#mcmc.merge.list from the merge step (see the example in function 'merge')
#ACS_ADS_global from the multi_ACS_ADS_global step (see the example in 'multi_ACS_ADS_global')
dataset.names = c("hb","hs","ht","ha","hi","hl",
                  "mb","ms","mt","ma","mi","ml")
mdsGlobal(ACS_ADS_global$ACS,dataset.names,sep="_",file="~/globalMDS.pdf")
}
}
