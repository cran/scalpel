% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCALPEL_stepFunctions.R
\name{scalpelStep3}
\alias{scalpelStep3}
\title{Perform Step 3 of SCALPEL.}
\usage{
scalpelStep3(
  step2Output,
  lambdaMethod = "trainval",
  lambda = NULL,
  minClusterSize = 1,
  alpha = 0.9,
  removeBorder = FALSE,
  excludeReps = NULL
)
}
\arguments{
\item{step2Output}{An object of class \code{scalpel} or \code{scalpelStep2}, which result from running the
\code{\link{scalpel}} or \code{\link{scalpelStep2}} functions, respectively.}

\item{lambdaMethod}{A description of how lambda should be chosen: either \code{"trainval"} (default),
\code{"distn"}, or \code{"user"}. A value of \code{"trainval"} means lambda will be chosen using a training/validation
set approach. A value of \code{"distn"} means lambda will be chosen as the negative of the 0.1\% quantile
of elements of active pixels (i.e., those contained in at least one dictionary element) of Y.
Using \code{"distn"} is computationally faster than \code{"trainval"}. Alternatively with \code{"user"},
the value of lambda can be directly specified using \code{lambda}.}

\item{lambda}{The value of lambda to use when fitting the sparse group lasso. By default, the value is automatically
chosen using the approach specified by \code{lambdaMethod}. If a value is provided for \code{lambda}, then \code{lambdaMethod}
will be ignored.}

\item{minClusterSize}{The minimum number of preliminary dictionary elements that a cluster must contain in order to be included
in the sparse group lasso. The default value is 1 (i.e., all possible dictionary elements are included).}

\item{alpha}{The value of alpha to use when fitting the sparse group lasso. The default value is 0.9.}

\item{removeBorder}{A logical scalar indicating whether the dictionary elements containing pixels in the 10-pixel
border of the video should be removed prior to fitting the sparse group lasso. The default value is \code{FALSE}.}

\item{excludeReps}{A vector giving the indices of which dictionary elements to exclude, where the indices refer
to the columns of \code{step2Output$A}. The default value is \code{NULL} and no dictionary elements are excluded.
Users may also specify \code{"discarded"}, which will exclude all dictionary elements discarded
using a previous call to \code{\link{reviewNeurons}} or \code{\link{reviewNeuronsInteractive}}.}
}
\value{
An object of class \code{scalpelStep3}, which can be summarized using \code{\link{summary}} and used with the plotting functions
\code{\link{plotResults}}, \code{\link{plotResultsAllLambda}}, \code{\link{plotSpatial}}, \code{\link{plotTemporal}}, and \code{\link{plotBrightest}}.
\itemize{
\item{\code{Afilter}: }{A matrix containing the filtered dictionary elements, where the ith column of \code{Afilter}
is a vector of 1's and 0's, indicating whether each pixel is contained in the ith filtered dictionary element.
Note that \code{Afilter} is equivalent to \code{A} after removing the components excluded due to being on the border
(if \code{removeBorder=TRUE}) or having fewer preliminary dictionary elements in their cluster than \code{minClusterSize}.}
\item{\code{Zhat}: }{A matrix containing the estimated temporal components, where the ith row of \code{Zhat}
is the estimated calcium trace corresponding to the ith spatial component (i.e., the ith column of \code{Afilter}).}
\item{\code{lambda}: }{The value of lambda used in fitting the sparse group lasso.}
\item{\code{ZhatList}: }{A list of matrices containing the estimated temporal components for alternative values of \code{lambda}
specified in \code{lambdaSeq}. These can be plotted using \code{\link{plotResultsAllLambda}}}.
\item{\code{lambdaSeq}: }{A vector with length equaling the length of \code{ZhatList}, where the ith element indicates the value of lambda
corresponding to the temporal components in \code{ZhatList[[i]]}.}
\item{\code{clustersUse}: }{A vector with length equaling the number of columns of \code{Afilter}, where the ith element indicates which column of
\code{step2Output$A} the ith column of \code{Afilter} equals.}
\item{Other elements: }{As specified by the user or returned from a previous step.}
}
}
\description{
This step involves spatial component selection and temporal component estimation. We estimate the temporal
components corresponding to the dictionary elements from Step 2 by solving a sparse group lasso problem
with a non-negativity constraint.
}
\details{
To solve the sparse group lasso problem in this step, we minimize the following over \code{Z} with all non-negative elements:

\code{0.5*sum((Y - AfilterTilde \%*\% Z)^2) + lambda*alpha*sum(Z)} \cr
\code{ + lambda*(1-alpha)*sum(sqrt(rowSums(Z^2)))}

where \code{AfilterTilde} is a scaled version of \code{Afilter}.

Several files containing data from this step and a summary of the step are saved in the folder
"outputFolder/Step1_version/Step2_omega_omega_cutoff_cutoff/Step3_lambdaMethod_lambdaMethod_
minClusterSize_minClusterSize_alpha_alpha_removeBorder_removeBorder" where \code{version} is a 5-digit unique
 identifier that is automatically generated in Step 1, \code{omega} and \code{cutoff} are the user-supplied parameters
 from Step 2, and \code{lambdaMethod}, \code{minClusterSize}, \code{alpha}, and \code{removeBorder} are
 the user-supplied parameters from this step. If dictionary elements were manually excluded using \code{excludeReps},
 this is appended to the folder name.
}
\examples{
\dontrun{
### many of the functions in this package are interconnected so the
### easiest way to learn to use the package is by working through the vignette,
### which is available at ajpete.com/software

#assumes you have run the example for the "scalpelStep2" function

#run Step 3 of SCALPEL
Step3Out = scalpelStep3(step2Output = Step2Out)
summary(Step3Out)
}
}
\seealso{
The entire SCALPEL pipeline can be implemented using the \code{\link{scalpel}} function. The
other steps in the pipeline can be run using the \code{\link{scalpelStep0}}, \code{\link{scalpelStep1}}, \code{\link{scalpelStep2}} functions.
Results from this step can be summarized using \code{\link{summary}}, loaded at a later time using \code{\link{getScalpelStep3}}, and plotted using
\code{\link{plotSpatial}}, \code{\link{plotTemporal}}, \code{\link{plotResults}}, and \code{\link{plotBrightest}}.
}
