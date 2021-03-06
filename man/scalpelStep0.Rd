% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCALPEL_stepFunctions.R
\name{scalpelStep0}
\alias{scalpelStep0}
\title{Perform Step 0 of SCALPEL.}
\usage{
scalpelStep0(
  outputFolder,
  rawDataFolder,
  videoHeight,
  fileType = "R",
  processSeparately = TRUE
)
}
\arguments{
\item{outputFolder}{The existing directory where the results should be saved.}

\item{rawDataFolder}{The directory where the raw data version of Y is saved. The data should be a
PxT matrix, where P is the total number of pixels per image frame and T
the number of frames of the video, for which the (i,j)th element contains the
fluorescence of the ith pixel in the jth frame. To create Y, you should
vectorize each 2-dimensional image frame by concatenating the columns of the image frame.  If the data is
saved in a single file, it should be named "Y_1.mat", "Y_1.rds", "Y_1.txt", or "Y_1.txt.gz" (depending on \code{fileType}),
and if the data is split over multiple files, they should be split into chunks of the columns
and named consecutively ("Y_1.mat", "Y_2.mat", etc.; "Y_1.rds", "Y_2.rds", etc.; "Y_1.txt", "Y_2.txt", etc.; or "Y_1.txt.gz", "Y_2.txt.gz", etc.).}

\item{videoHeight}{The height of the video (in pixels).}

\item{fileType}{Indicates whether raw data is an .rds (default value; \code{fileType="R"}), .mat (\code{fileType="matlab"}), .txt (\code{fileType="text"}), or
.txt.gz (\code{fileType="zippedText"}) file. Any text files should not have row or column names.}

\item{processSeparately}{Logical scalar giving whether the multiple raw data files should be
processed individually, versus all at once. Processing the files separately may be preferable for larger videos.
The default value is \code{TRUE}; this argument is ignored if the raw data is saved in a single file.}
}
\value{
An object of class \code{scalpelStep0}, which can be summarized using \code{\link{summary}}, used to run SCALPEL Step 1 using \code{\link{scalpelStep1}},
or can be used with the plotting functions \code{\link{plotFrame}}, \code{\link{plotThresholdedFrame}}, and \code{\link{plotVideoVariance}}.
\itemize{
\item{\code{minRaw, maxRaw, minDeltaf, maxDeltaf}: }{Minimum and maximum values for the raw and processed videos.}
\item{\code{partsRaw, partsDeltaf}: }{Vectors indicating the indices of the raw and processed data files, respectively.}
\item{\code{nFramesRaw, nFramesDeltaf}: }{The number of frames in each part of the raw and processed data.}
\item{\code{lowThreshold, highThreshold}: }{The default lowest and highest threshold values for image segmentation that may be used in Step 1.}
\item{Other elements: }{As specified by the user.}
}
}
\description{
This step involves data pre-processing. We read in the raw data version of Y and perform standard pre-processing
techniques in order to smooth the data both temporally and spatially, remove the bleaching effect,
and calculate a standardized fluorescence.
}
\details{
Several files containing data from this step and a summary of the step are saved in "outputFolder".
}
\examples{
\dontrun{
### many of the functions in this package are interconnected so the
### easiest way to learn to use the package is by working through the vignette,
### which is available at ajpete.com/software

#existing folder to save results (update this to an existing folder on your computer)
outputFolder = "scalpelResultsStepByStep"
#location on computer of raw data in R package to use
rawDataFolder = gsub("Y_1.rds", "", system.file("extdata", "Y_1.rds", package = "scalpel"))
#video height of raw data in R package
videoHeight = 30
#run Step 0 of SCALPEL
Step0Out = scalpelStep0(outputFolder = outputFolder,
                       rawDataFolder = rawDataFolder, videoHeight = videoHeight)
summary(Step0Out)
}
}
\seealso{
The entire SCALPEL pipeline can be implemented using the \code{\link{scalpel}} function. The
other steps in the pipeline can be run using the \code{\link{scalpelStep1}}, \code{\link{scalpelStep2}}, \code{\link{scalpelStep3}} functions.
Results from this step can be summarized using \code{\link{summary}}, loaded at a later time using \code{\link{getScalpelStep0}}, and plotted using
\code{\link{plotFrame}}, \code{\link{plotThresholdedFrame}}, and \code{\link{plotVideoVariance}}.
}
