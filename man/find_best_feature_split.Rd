% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_search_functions.R
\name{find_best_feature_split}
\alias{find_best_feature_split}
\title{Find best split value for certain feature with weighted observations}
\usage{
find_best_feature_split(feature, object.subdtree, datasub, weightssub,
  minleafsize, misclassification_weights)
}
\arguments{
\item{feature}{Character name of feature that will be split on.}

\item{object.subdtree}{Object of class dtree. The split of the rootnode will be optimized.}

\item{datasub}{A dataframe: all observations that reach the rootnode of the subtree}

\item{minleafsize}{The minimum number of observations in any terminal <leaf> node.}
}
\value{
A list of two, where object.bestsubdtree is an object of class dtree, misclass is the total
number of misclassified observations, and splitpossible is a logical value that indicated whether a
split on this feature is possible.
}
\description{
For a specific branchnode (root of inputted tree) and feature, find the optimal split.
}
