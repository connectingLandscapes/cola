% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cola_fun.R
\name{runCDPOP}
\alias{runCDPOP}
\title{Run CDPOP model}
\usage{
runCDPOP(
  py = Sys.getenv("COLA_PYTHON_PATH"),
  inputvars = NULL,
  agevars = NULL,
  cdmat = NULL,
  xy = NULL,
  tempFolder
)
}
\arguments{
\item{py}{Python location}
}
\value{
Path with CDPOP results
}
\description{
Run CDPOP model
}
\examples{
runCDPOP( )
}
\author{
Ivan Gonzalez <ig299@nau.edu>

Patrick Jantz <Patrick.Jantz@gmail.com>
}
