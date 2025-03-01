% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Olink_plate_randomizer.R
\name{olink_plate_randomizer}
\alias{olink_plate_randomizer}
\title{Randomly assign samples to plates}
\usage{
olink_plate_randomizer(
  Manifest,
  PlateSize = 96,
  SubjectColumn,
  iterations = 500,
  available.spots,
  seed
)
}
\arguments{
\item{Manifest}{tibble/data frame in long format containing all sample ID's. Sample ID column must be named SampleID.}

\item{PlateSize}{Integer. Either 96 or 48. 96 is default.}

\item{SubjectColumn}{(Optional) Column name of the subject ID column. Cannot contain missings. If provided, subjects are kept on the same plate.}

\item{iterations}{Number of iterations for fitting subjects on the same plate.}

\item{available.spots}{Integer. Number of wells available on each plate. Maximum 40 for T48 and 88 for T96. Can also take a vector equal to the number of plates to be used indicating the number of wells available on each plate.}

\item{seed}{Seed to set. Highly recommend setting this for reproducibility.}
}
\value{
Tibble including SampleID, SubjectID etc. assigned to well positions.
}
\description{
Samples can be randomly assigned to plates using base::sample with an option to keep Subjects on the same plate.  Olink Data Science no longer recommends forced balanced randomization considering other clinical variables.
}
\examples{
\donttest{randomized.manifest <- olink_plate_randomizer(manifest, seed=12345)}
\donttest{randomized.manifest <- olink_plate_randomizer(manifest,SubjectColumn="SubjectID",
                                                        available.spots=c(88,88), seed=12345)}
}
\keyword{plates}
\keyword{randomized}
