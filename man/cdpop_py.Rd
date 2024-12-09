% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cola_fun.R
\name{cdpop_py}
\alias{cdpop_py}
\title{Runs CDOPOP2}
\usage{
cdpop_py(
  outshp,
  outtif,
  param5 = 0.5,
  param6 = 1000,
  tif,
  incrk,
  inlcc,
  maskedcsname = paste0(tempfile(), ".tif"),
  py = Sys.getenv("COLA_PYTHON_PATH")
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
