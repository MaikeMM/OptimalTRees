% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finalization.R
\name{batch_tuneCP}
\alias{batch_tuneCP}
\title{Find optimal parameter alpha for OCT}
\usage{
batch_tuneCP(trees, losses, trainingdata, validationdata, trainingweights,
  validationweights, misclassification_weights)
}
\arguments{
\item{trees}{A list of dtree objects}

\item{losses}{A numeric array of corresponding values of loss functions of the dtree objects in trees,
found with alpha = 0}

\item{trainingdata}{Dataframe that will be used to train the OCT}

\item{validationdata}{Dataframe that will be used to validate the OCT}
}
\value{
A list of two, vbest is the best possible misclassification rate on the validation data and alphabest
is the best corresponding value for parameter alpha.
}
\description{
For a batch of tree, construct a mean curve of validation error as a function of complexity parameter
such that the optimal parameter alpha (smallest validation error) can be found.
}
