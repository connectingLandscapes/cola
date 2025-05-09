% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cola_fun.R
\name{randPtsFun}
\alias{randPtsFun}
\title{Create random points}
\usage{
randPtsFun(rvect, npts, rmin, rmax)
}
\arguments{
\item{rvect}{Raster. Raster object co be sampled}

\item{npts}{Numeric. Number of points}

\item{rmin}{Numeric. Raster minimum values o¿to consider}

\item{rmax}{Numeric. Raster maximun value to consider}
}
\value{
Path with CDPOP results
}
\description{
Create random points
}
\examples{
hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
newPoints <- randPtsFun(terra::rast(hs), 10, 0.2, 0.8)
}
\author{
Ivan Gonzalez <ig299@nau.edu>

Patrick Jantz <Patrick.Jantz@gmail.com>
}
