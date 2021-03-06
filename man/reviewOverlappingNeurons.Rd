% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCALPEL_reviewNeurons.R
\name{reviewOverlappingNeurons}
\alias{reviewOverlappingNeurons}
\title{Save additional frames for overlapping neurons from SCALPEL.}
\usage{
reviewOverlappingNeurons(scalpelOutput, neuronSet, numFrames = 10)
}
\arguments{
\item{scalpelOutput}{An object returned by one of the SCALPEL functions:
\code{\link{scalpel}}, \code{\link{scalpelStep2}}, or \code{\link{scalpelStep3}}.}

\item{neuronSet}{The set of neurons that should be reviewed:
use \code{"A"} for those resulting from \code{\link{scalpelStep2}} and saved as \code{scalpelOutput$A}, or use \code{"Afilter"} for those resulting from
\code{\link{scalpelStep3}} and saved as \code{scalpelOutput$Afilter}. This argument is ignored if the class of \code{scalpelOutput} is \code{scalpelStep2}.}

\item{numFrames}{The maximum number of frames that should be saved for each neuron being considered. Each neuron has a number of frames
equal to the number of members in that neuron's cluster that can be plotted. All frames will be saved when the total number of available frames for the neuron
is less than \code{numFrames}. The default value is 10.}
}
\value{
None
}
\description{
We use this function after running \code{\link{reviewNeurons}} or \code{\link{reviewNeuronsInteractive}} to plot additional frames for neurons that overlap with others.
These frames are saved, and the classification for the neurons can then be updated using \code{\link{updateNeurons}} or \code{\link{updateNeuronsInteractive}}.
}
\examples{
\dontrun{
### many of the functions in this package are interconnected so the
### easiest way to learn to use the package is by working through the vignette,
### which is available at ajpete.com/software

#assumes you have run the example for the "reviewNeuronsInteractive" function

reviewOverlappingNeurons(scalpelOutput = scalpelOutput, neuronSet = "A")
}
}
\seealso{
\code{\link{reviewNeurons}}, \code{\link{updateNeurons}}, \code{\link{reviewNeuronsInteractive}}, \code{\link{updateNeuronsInteractive}}
}
