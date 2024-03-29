% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readData.R
\name{read.rawData}
\alias{read.rawData}
\title{Read raw data: microarray or RNAseq.}
\usage{
read.rawData(data.file, sep = ",", quote = "\\"", header = T)
}
\arguments{
\item{data.file:}{the file name of the data to be read.}

\item{sep:}{the field separator character. Values on each line of
the file are separated by this character. If \code{sep = ""}, the
separator is 'white space', that is one or more spaces, tabs,
newlines or carriage returns. The default is \code{,}. See
\link{read.csv}.}

\item{quote:}{the set of quoting characters. To disable quoting
altogether, use \code{quote = ""}.  See \link{scan} for the
behaviour on quotes embedded in quotes. Quoting is only
considered for columns read as character, which is all of them
unless colClasses is specified.  The default is \code{\"}. See
\link{read.csv}}

\item{header:}{a logical value indicating whether the file contains
the names of the variables as its first line. If missing, the
value is determined from the file format: \code{header} is set to
\code{TRUE} if and only if the first row contains one fewer field
than the number of columns. The default is \code{TRUE}. See
\link{read.csv}.}
}
\value{
the data matrix.
}
\description{
Read raw data: microarray & RNAseq
The \code{read.rawData} is function to read raw data
}
\examples{
\dontrun{
data(hb)
write.csv(data,file="rawData.csv")
data <- read.rawData(data.file="rawData.csv")
}
}
