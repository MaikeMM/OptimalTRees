#' CART procedure - random subset of features considered
#'
#' @param data A dataframe: all character variables must be factors, and the dependent
#' variable must be a factor.
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @param maxdepth Set the maximum depth of any node of the final tree,
#' with the root node counted as depth 0.
#' @param features The features that should be considered. If empty the first n - 1 column
#' names of the dataframe will be taken as features.
#' @param dep_var The dependent variable. If empty the last column
#' name of the dataframe will be taken as dependent variable.
#' @return An object of class dtree
makeCART <- function(data, weights = NULL, minleafsize = 1, maxdepth = 6, features = NULL, dep_var = NULL){
  if (missing(data)) stop('data is missing with no default')

  if (is.null(features)){
    features <- head(names(data), -1)
    warning('No features were specified,')
  }
  if (is.null(dep_var)){
    dep_var <- tail(names(data), 1)
  }
  if (is.null(weights)){
    weights <- rep(1, nrow(data))
  }

  if (length(weights) != nrow(data)) stop('the number of data rows and the length of weights must be equal')
  smalldata <- select(data, c(features, dep_var))
  if (any(!complete.cases(smalldata))) stop('data cannot cannot contain incomplete cases')
  if (minleafsize < 1) stop('minleafsize must be positive')
  if (maxdepth < 1) stop('maxdepth must be positive')

  make_best_ginisplit <- function(data, weights, subsetfeatures, minleafsize, features, dep_var){
    numberfeatures <- length(features)

    splitexists <- F
    ginibest <- Inf

    A <- matrix(0, nrow = numberfeatures, ncol = 1); bvec <- 0; cats <- vector('list', 1)
    object.dtree <- new('dtree', A = A, bvec = bvec, cats = cats,
                        dep_var = dep_var, leaf_classes = character(), features = features)

    N <- dim(data)[1]
    dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F)
    if (length(unique(dep_values)) == 1) return(list(object.dtree = object.dtree, splitpossible = F))


    smallwhere <- tapply(weights, dep_values, sum, simplify = T); smallwhere[is.na(smallwhere)] <- 0
    bestadjustedgini <-  sum(smallwhere ** 2) / sum(smallwhere)

    splitpossible <- F
    for (feature in subsetfeatures){
      featindex <- which(feature == features)
      object.dtree@A[featindex, 1] <- 1

      # Order all data
      values <- unlist(select(data, feature), use.names = F)
      ord <- order(values); ordered_values <- values[ord]; ordered_dep_values <- dep_values[ord]
      ordered_table <- table(ordered_values); ordered_weights <- weights[ord]

      m <- length(unique(values))
      # If m equals 1, there cannot be a split
      if (m > 1) {
        if (!is.factor(values)){
          # If here, there are multiple values in unique_sorted_values and this is a numeric array

          # Initialize variables
          classes <- unique(dep_values)
          where <- matrix(0, nrow = 2, ncol = length(classes)); colnames(where) <- classes

          where[1, ] <- 0
          where[2, ] <- tapply(weights, dep_values, sum, simplify = T)[classes]

          # Now loop through all observations and consider the new 'where' matrix if the split was such that
          # a particular observation (and all observation with a lower value at specific feature) follow the
          # left branch. Count the misclassification and minimum leaf size to check eligibility.
          index <- 1
          for (j in 1:(length(ordered_table)-1)){

            indexset <- index:(index + ordered_table[j] - 1)
            localtable <- tapply(ordered_weights[indexset], ordered_dep_values[indexset], sum, simplify = T)[classes]
            localtable[is.na(localtable)] <- 0

            # Remove observation from one of the right leafs and add it to one of the left leafs
            where[1, ] <- where[1, ] + localtable
            where[2, ] <- where[2, ] - localtable

            # Count the adjusted gini gain
            adjustedgini <- sum(where[1, ] ** 2) / sum(where[1, ]) + sum(where[2, ] ** 2) / sum(where[2, ])

            # Count the minimum leaf size
            minleafsizelocal <- min(rowSums(where))
            if (adjustedgini > bestadjustedgini & minleafsizelocal >= minleafsize){
              # If here then the misclassification is better than the current and the minimum leaf size is greater
              # than the prescribed minimum
              bestadjustedgini <- adjustedgini
              b <- mean(as.numeric(names(ordered_table[j:(j+1)])))
              object.dtree@bvec[1] <- b; object.dtree@cats[1] <- vector('list', 1)
              bestfeatindex <- featindex

              splitpossible <- T
            }
            index <- index + ordered_table[j]
          }
        } else {
          # If here, the feature is a categorical feature and we must consider combinations of different levels
          # as possible splits
          levelsfeat <- unique(droplevels(ordered_values)); numberlevels <- length(levelsfeat)
          ordered_values <- droplevels(ordered_values)
          classes <- unique(droplevels(ordered_dep_values))
          crosstable <- t(sapply(X = levelsfeat, FUN = function(x) {
            indexset <- which(ordered_values == x)
            tab <- tapply(ordered_weights[indexset], ordered_dep_values[indexset], sum, simplify = T)[classes];
            tab[is.na(tab)] <- 0
            tab}, simplify = T))
          rownames(crosstable) <- levelsfeat;

          where <- matrix(0, nrow = 2, ncol = length(classes)); colnames(where) <- classes

          # Loop over all possible combinations with outer and inner loop
          outeroptions <- floor(numberlevels/2)
          for (n in 1:outeroptions){
            options <- combn(levelsfeat, n); numberoptions <- dim(options)[2]
            for (i in 1:numberoptions){
              catsplit <- options[, i]

              where[1, ] <- colSums(crosstable[catsplit, , drop = F])
              where[2, ] <- colSums(crosstable[levelsfeat[!(levelsfeat %in% catsplit)], , drop = F ])

              # Count the adjusted gini gain
              adjustedgini <- sum(where[1, ] ** 2) / sum(where[1, ]) + sum(where[2, ] ** 2) / sum(where[2, ])

              # Count the minimum leaf size
              minleafsizelocal <- min(rowSums(where))

              if (adjustedgini > bestadjustedgini & minleafsizelocal >= minleafsize){
                # If here then the misclassification is better than the current and the minimum leaf size is greater
                # than the prescribed minimum
                bestadjustedgini <- adjustedgini
                object.dtree@bvec[1] <- 0; object.dtree@cats[[1]] <- catsplit
                bestfeatindex <- featindex
                splitpossible <- T
              }
            }
          }
        }
      }
      object.dtree@A[featindex, 1] <- 0
    }
    if (splitpossible){
      object.dtree@A[bestfeatindex, 1] <- 1
    }
    return(list(object.dtree = object.dtree, splitpossible = splitpossible))
  }



  numberfeatures <- length(features)
  numberfeatures_split <- floor(sqrt(numberfeatures))

  subsetfeatures <- sample(features, numberfeatures_split)
  firstsplit <- make_best_ginisplit(data, weights, subsetfeatures, minleafsize, features = features, dep_var = dep_var)
  object.dtree <- firstsplit$object.dtree

  treeDone <- F
  while(!treeDone){
    Z <- find_pos_points(object.dtree, data)

    numbernodes <- dim(Z)[2];
    depthtree <- log(numbernodes + 1)/log(2) - 1

    if (depthtree == maxdepth) break

    leafnodes <- (2 ^ depthtree) : (2 ^ (depthtree + 1) - 1)

    colsA <- matrix(0, nrow = numberfeatures, ncol = length(leafnodes))
    elsbvec <- rep(0, length(leafnodes))
    elscats <- vector('list', length(leafnodes))

    object.dtree@A <- cbind(object.dtree@A, colsA);
    object.dtree@bvec <- c(object.dtree@bvec, elsbvec);
    object.dtree@cats <- c(object.dtree@cats, elscats)

    for (node in leafnodes){
      pnode <- find_parent_node(node)
      if (!all(object.dtree@A[, pnode] == 0)){
        subsetfeatures <- sample(features, numberfeatures_split)

        indexset <- which(Z[, node] == 1); datasub <- data[indexset, ]; weightssub <- weights[indexset]
        newsplit <- make_best_ginisplit(datasub, weightssub, subsetfeatures, minleafsize, features, dep_var)
        if (newsplit$splitpossible){
          object.dtree@A[, node] <- newsplit$object.dtree@A;
          object.dtree@bvec[node] <- newsplit$object.dtree@bvec;
          object.dtree@cats[node] <- newsplit$object.dtree@cats;
        }
      }
    }

    if (all(object.dtree@A[, leafnodes] == 0)) {
      treeDone <- T
      object.dtree <- cut_tree_fun(object.dtree)
    }
  }
  object.dtree@leaf_classes <- find_class_leafs_and_misclass_and_min_leaf_size_fun(object.dtree, data)$leaf_classes
  return(object.dtree)
}


