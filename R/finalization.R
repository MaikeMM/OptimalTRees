#' Create the optimal classification tree (OCT)
#'
#' This function created an interpretable decision tree to classify data. It created batches of trees of specific
#' depths to optimize the parameters alpha (complexity parameter) and maximum depth of the tree, then creates
#' a new batch of trees with these optimal parameters. The best performing tree from this batch is returned.
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data A dataframe: all character variables must be factors, and the dependent
#' variable must be a factor. NA rows will be deleted.
#' @param weights weights
#' @param maxdepth Set the maximum depth of any node of the final tree,
#' with the root node counted as depth 0.
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @param numbertries The number of searches that will be performed at each possible depth,
#' to tune the parameters and find the optimal tree.
#' @param misclassification_weights
#' @return An object of class dtree that is the optimal classification tree
#' @export
optimaltrees <- function(formula, data, weights = NULL, maxdepth = 6, minleafsize = 1, numbertries = 10, misclassification_weights = NULL){
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (nrow(data) != length(weights)) stop("dataframe row dimension and weights length must be the same")
  if (minleafsize < 0) stop('minleafsize must be nonnegative')
  if (maxdepth < 1) stop('maxdepth must be positive')
  if (numbertries < 2) stop('numbertries must be greater than 1')
  # Obtain features and dependent variable from formula input
  mod <- as.formula(formula)
  dep_var <- all.vars(mod)[1]
  features <- tail(all.vars(mod), -1)
  if (features[1] == "."){
    features <- (names(data))[names(data) != dep_var]
  }

  # Prepare the data (normalize numeric features and let character features be factors)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  originaldata <- data
  res <- prepare_dataset(data, weights, features, dep_var)
  data <- res$data; weights <- res$weights; minmax <- res$minmax

  dep_values <- unlist(select(data, dep_var), use.names = F); classes <- unique(dep_values);
  if (is.null(misclassification_weights)) {
    misclassification_weights <- matrix(1, nrow = length(classes), ncol = length(classes)) - diag(length(classes))
    rownames(misclassification_weights) <- classes; colnames(misclassification_weights) <- classes
  } else {
    if (dim(misclassification_weights)[1] != dim(misclassification_weights)[2]) stop('misclassification_weights must be square')
    if (any(!(classes %in% colnames(misclassification_weights))) | any(!(classes %in% rownames(misclassification_weights)))) {
      stop('column and rownames of misclassification_weights must contain all possible classes')
    }
  }

  # Split data in training and validation set
  indexset <- sample(1:(dim(data)[1]), (dim(data)[1])*0.7)
  trainingdata <- data[indexset, ]; validationdata <- data[-indexset, ]
  trainingweights <- weights[indexset]; validationweights <- weights[-indexset]

  # Find optimal parameters for alpha and depth
  optimalparams <- tune(trainingdata, validationdata, trainingweights, validationweights, features, dep_var, maxdepth, minleafsize, numbertries, misclassification_weights)

  # Find new batch of trees with optimal parameters
  cat('final batch with alpha = ', optimalparams$alphabest, ' and D = ', optimalparams$Dbest, 'at ', format(Sys.time(), "%X"), '\n')
  results <- mclapply(X = 1:numbertries, FUN = function(i){
    CART <- makeCART(data = trainingdata, trainingweights, maxdepth = optimalparams$Dbest, minleafsize = minleafsize, features = features, dep_var = dep_var)
    newtree <- local_search(object.dtree = CART, data = trainingdata, trainingweights, alpha = optimalparams$alphabest, maxdepth = optimalparams$Dbest,
                            minleafsize = minleafsize, misclassification_weights = misclassification_weights)
  })

  trees <- vector('list', numbertries); losses <- rep(Inf, numbertries)
  for (try in 1:numbertries){
    trees[[try]] <- results[[try]]$object.dtree
    losses[try] <- results[[try]]$loss
  }

  # For each tree in batch, find misclassification on full dataset
  misclasspertree <- function(object.dtree, data, weights){
    find_where_and_leafs(object.dtree, data, weights, misclassification_weights)$weighted_misclass
  }
  mis_fulldata <- sapply(trees, misclasspertree, data = data, weights = weights)
  mis_trainingdata <- sapply(trees, misclasspertree, data = trainingdata, weights = trainingweights)
  mis_validationdata <- sapply(trees, misclasspertree, data = validationdata, weights = validationweights)
  # Select tree with lowest misclassification
  #index <- which(mis == min(mis))[1]
  #object.besttree <- trees[[index]]

  # Denormalize split values of optimal tree
  denorm <- function(object.besttree) {
    for (i in 1:length(object.besttree@bvec)){
      b <- object.besttree@bvec[i]
      if (b != 0){
        feature <- features[which(object.besttree@A[, i] == 1)]
        index <- which(names(data) == feature)
        newb <- b * minmax[2, index] + (1 - b) * minmax[1, index]
        object.besttree@bvec[i] <- newb
      }
    }
    return(object.besttree)
  }
  res <- lapply(trees, denorm)
  original_trainingdata <- originaldata[indexset, ]; original_validationdata <- originaldata[-indexset, ]
  original_trainingdata$weights <- trainingweights; original_validationdata$weights <- validationweights
  return(list(trees = res, mis_fulldata = mis_fulldata, mis_trainingdata = mis_trainingdata,
              mis_validationdata = mis_validationdata, trainingdata = original_trainingdata, validationdata = original_validationdata))
}

#' Find optimal parameters (alpha and depth) for OCT
#'
#' Create multiple batches of trees of varying depths, for which then the optimal parameter alpha and corresponding
#' validation error is found. The global smallest validation error decides the optimal depth, and optimal alpha corresponding
#' to that depth.
#' @param trainingdata Dataframe that will be used to train the OCT
#' @param validationdata Dataframe that will be used to validate the OCT
#' @param trainingweights Training weights
#' @param validationweights Validation wights
#' @param features Character array of the features that will be used in training the OCT
#' @param dep_var Character that indicated the dependent variable (class name)
#' @param maxdepth Set the maximum depth of any node of the final tree,
#' with the root node counted as depth 0.
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @param numbertries The number of searches that will be performed at each possible depth,
#' to tune the parameters and find the optimal tree.
#' @param misclassification_weights
#' @return A list of two - $vbest gives the misclassification on the validation set at
#' optimal alpha $bestalpha
#' @export
tune <- function(trainingdata, validationdata, trainingweights, validationweights, features, dep_var, maxdepth, minleafsize, numbertries, misclassification_weights){
  # Decide how many trees will be selected for pruning
  batchsize <- max(2, floor(numbertries * 0.1))

  # For each depth, make a batch of 'numbertries' trees and select 'batchsize' of these with the lowest loss.
  # With these selected trees, find optimal parameter alpha and corresponding misclassification on
  # validation set for the specific depth. Save alpha and depth with lowest misclassification rate.
  vbest <- Inf
  for (D in 1:maxdepth){
    cat('at depth ', D, 'at ', format(Sys.time(), "%X"), '\n')
    # Create batch of 'numbertries' trees and select 'batchsize' of these with the lowest loss.
    treesatdepth <- vector('list', batchsize); lossesatdepth <- rep(Inf, batchsize)
    if (D == 1){
      # If depth == 1, the same tree will always be found
      CART <- makeCART(data = trainingdata, trainingweights, maxdepth = D, minleafsize = minleafsize, features = features, dep_var = dep_var)
      newtree <- local_search(object.dtree = CART, data = trainingdata, weights =  trainingweights, alpha = 0, minleafsize = minleafsize,
                              maxdepth = D, misclassification_weights = misclassification_weights)
      for (i in 1:batchsize) {
        treesatdepth[[i]] <- newtree$object.dtree
        lossesatdepth[[i]] <- newtree$loss
      }
    } else {
      results <- mclapply(X = 1:numbertries, FUN = function(i){
        CART <- makeCART(data = trainingdata, trainingweights, maxdepth = D, minleafsize = minleafsize, features = features, dep_var = dep_var)
        newtree <- local_search(object.dtree = CART, data = trainingdata, trainingweights, alpha = 0, maxdepth = D,
                                minleafsize = minleafsize, misclassification_weights = misclassification_weights)
        })
      indexset <- order(sapply(results, function(ls){ls$loss}))[1:batchsize]
      for (i in 1:batchsize){
        treesatdepth[[i]] <- results[[indexset[i]]]$object.dtree
        lossesatdepth[i] <- results[[indexset[i]]]$loss
      }
      #if (newtree$loss < max(lossesatdepth)){
      #  index <- which(lossesatdepth == max(lossesatdepth))[1]
      #  treesatdepth[[index]] <- newtree$object.dtree
      #  lossesatdepth[index] <- newtree$loss
      #}
    }
    # Find optimal alpha and corresponding misclassficiation
    params <- batch_tuneCP(trees = treesatdepth, losses = lossesatdepth, trainingdata = trainingdata,
                          validationdata = validationdata, trainingweights, validationweights, misclassification_weights = misclassification_weights)

    # Update parameters
    if (params$vbest < vbest) {
      vbest <- params$vbest; alphabest <- params$alphabest; Dbest <- D
    }
  }
  return(list(alphabest = alphabest, Dbest = Dbest))
}

#' Find optimal parameter alpha for OCT
#'
#' For a batch of tree, construct a mean curve of validation error as a function of complexity parameter
#' such that the optimal parameter alpha (smallest validation error) can be found.
#'
#' @param trees A list of dtree objects
#' @param losses A numeric array of corresponding values of loss functions of the dtree objects in trees,
#' found with alpha = 0
#' @param trainingdata Dataframe that will be used to train the OCT
#' @param validationdata Dataframe that will be used to validate the OCT
#' @return A list of two, vbest is the best possible misclassification rate on the validation data and alphabest
#' is the best corresponding value for parameter alpha.
batch_tuneCP <- function(trees, losses, trainingdata, validationdata, trainingweights, validationweights, misclassification_weights){
  # Prepare necessary numeric arrays
  alpha_long <- seq(0, 1, length.out = 1001)
  fm <- matrix(0, nrow = length(losses), ncol = 1001)

  # Define necessary function to construct the curve of validation error as a function of complexity parameter
  valcurve_inside_sum <- function(i, alpha_longel, alpha, V) {
    if (i == 1) {
      psum <- V[i] * as.numeric(alpha_longel >= alpha[i])
    } else {
      psum <- (V[i] - V[i-1]) * as.numeric(alpha_longel >= alpha[i])
    }
    return(psum)
  }
  valcurve <- function(alpha_longel, alpha, V){
    endi <- length(alpha)
    psums <- sapply(1:endi, valcurve_inside_sum, alpha_longel = alpha_longel, alpha = alpha, V = V)
    return(sum(psums))
  }

  # For each tree in the batch, find the curve of validation error as a function of complexity parameter
  orderset <- order(losses)
  for (index in orderset){
    tree <- trees[[index]]
    result_prune <- prune(tree, trainingdata, validationdata, trainingweights, validationweights, misclassification_weights)
    alpha <- result_prune$alpha; V <- result_prune$V
    fm[index,] <- sapply(alpha_long, valcurve, alpha = alpha, V = V)
  }

  # Take the mean of all curves to find optimal alpha and corresponding validation
  f <- apply(fm, 2, mean)
  vbest <- min(f)
  alphabestvec <- alpha_long[which(f == vbest)]
  alphabest <- (min(alphabestvec) + max(alphabestvec))/2
  return(list(vbest = vbest, alphabest = alphabest))
}

#' Prune OCT
#'
#' For each tree, there is a corresponding interval for alpha from which any value will generate this tree.
#' This procedure finds the discrete values of alpha for which the tree changes, and for each range of alpha,
#' finds the validation error.
#'
#' @param object.dtree An object of class dtree
#' @param trainingdata Dataframe that will be used to train the OCT
#' @param validationdata Dataframe that will be used to validate the OCT
#' @return A list of two: vector alpha of critical values for complexity parameter and vector of corresponding
#' validation errors for each critical value
prune <- function(object.dtree, trainingdata, validationdata, trainingweights, validationweights, misclassification_weights){
  alpha <- 0
  V <- find_where_and_leafs(object.dtree, validationdata, validationweights, misclassification_weights)$weighted_misclass
  dep_values <- unlist(select(trainingdata, object.dtree@dep_var), use.names = F); classes <- unique(dep_values)
  where <- tapply(trainingweights, dep_values, sum, simplify = T)[classes]
  which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights))); W <- misclassification_weights[which_W, which_W];
  Lhat <- min(where %*% t(W))

  i <- 1
  while (complexity_tree_fun(object.dtree) > 0) {
    i <- i + 1
    alpha[i] <- Inf

    trainingZ <- find_pos_points(object.dtree, trainingdata)

    numberbranchnodes <- length(object.dtree@bvec)
    branchnodes<- (1:numberbranchnodes)[apply(object.dtree@A, 2, function(x) {any(x == 1)})]
    for (node in branchnodes){
      subtree <- find_subtree(rootnode = node, object.dtree); object.subdtree <- subtree$object.subdtree
      indexset <- which(trainingZ[, node] == 1); subtrainingdata <- trainingdata[indexset, ]; subtrainingweights <- trainingweights[indexset]

      dep_values <- unlist(select(subtrainingdata, object.dtree@dep_var), use.names = F)
      classes <- unique(dep_values); where <-  tapply(subtrainingweights, dep_values, sum, simplify = T)[classes]
      which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights))); W <- misclassification_weights[which_W, which_W];

      misclass_leaf <- min(where %*% t(W))
      misclass_subtree <- find_where_and_leafs(object.subdtree, data = subtrainingdata, subtrainingweights, misclassification_weights)$weighted_misclass
      complexity_subtree <- complexity_tree_fun(object.subdtree)

      alpha_local <- (misclass_leaf - misclass_subtree) / (complexity_subtree * Lhat)

      if (alpha_local < alpha[i]){
        nodehat <- node
        alpha[i] <- alpha_local
      }
    }

    subtree_to_be_removed <- find_subtree(rootnode = nodehat, object.dtree)
    for (branchnode in subtree_to_be_removed$branchnodes){
      object.dtree@A[, branchnode] <- 0
      object.dtree@bvec[branchnode] <- 0
      object.dtree@cats[branchnode] <- vector('list', 1)
    }
    object.dtree <- cut_tree_fun(object.dtree)
    V[i] <- find_where_and_leafs(object.dtree, data = validationdata, validationweights, misclassification_weights)$weighted_misclass
  }
  return(list(alpha = alpha, V = V))
}
