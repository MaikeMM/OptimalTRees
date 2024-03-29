% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finalization.R
\name{tune}
\alias{tune}
\title{Find optimal parameters (alpha and depth) for OCT}
\usage{
tune(trainingdata, validationdata, trainingweights, validationweights,
  features, dep_var, maxdepth, minleafsize, numbertries,
  misclassification_weights)
}
\arguments{
\item{trainingdata}{Dataframe that will be used to train the OCT}

\item{validationdata}{Dataframe that will be used to validate the OCT}

\item{trainingweights}{Training weights}

\item{validationweights}{Validation wights}

\item{features}{Character array of the features that will be used in training the OCT}

\item{dep_var}{Character that indicated the dependent variable (class name)}

\item{maxdepth}{Set the maximum depth of any node of the final tree,
with the root node counted as depth 0.}

\item{minleafsize}{The minimum number of observations in any terminal <leaf> node.}

\item{numbertries}{The number of searches that will be performed at each possible depth,
to tune the parameters and find the optimal tree.}
}
\value{
A list of two - $vbest gives the misclassification on the validation set at
optimal alpha $bestalpha
}
\description{
Create multiple batches of trees of varying depths, for which then the optimal parameter alpha and corresponding
validation error is found. The global smallest validation error decides the optimal depth, and optimal alpha corresponding
to that depth.
}
