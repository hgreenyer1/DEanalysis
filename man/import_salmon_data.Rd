% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importData.R
\name{import_salmon_data}
\alias{import_salmon_data}
\title{Import Salmon Data}
\usage{
import_salmon_data(salmon_dir, tx2gene, out_dir = NULL)
}
\arguments{
\item{salmon_dir}{Directory where salmon output file folders are located}

\item{tx2gene}{A txdb object}

\item{out_dir}{option to specify output directory}
}
\value{
a list containing two objects: an annotations df and corresponding salmon data extracted via txi_salmon
}
\description{
Import Salmon Data
}
\examples{
import_salmon_data('~/salmonDirectory', tx2gene)
}
