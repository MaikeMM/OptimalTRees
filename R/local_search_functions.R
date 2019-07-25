#' Perform local search procedure to find OCT
#'
#' The local search procedure designed by Dunn takes a CART decision tree as input, and then exhaustively
#' considers changes to each node to optimize the tree.
#'
#' @param object.dtree Object of class dtree, usually a tree created by CART procedure
#' where each split only a random subset of features is considered.
#' @param data  A dataframe: all character variables must be factors, and the dependent
#' variable must be a factor.
#' @param alpha A complexity parameter
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @param maxdepth  Set the maximum depth of any node of the final tree,
#' with the root node counted as depth 0.
#' @return A list of two, where the first element is an object of class dtree and the second
#' element is the loss associated with the dtree.
local_search <- function(object.dtree, data, weights = NULL, alpha = 0, minleafsize = 1, maxdepth = 6, misclassification_weights = NULL){
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (nrow(data) != length(weights)) stop("dataframe row dimension and weights length must be the same")
  if (class(object.dtree) != 'dtree') stop('object.dtree must be of class dtree')
  if (any(is.na(data))) stop('NA values in data not allowed')
  if (alpha < 0) stop('alpha must be nonnegative')
  if (minleafsize < 0) stop('minleafsize must be nonnegative')
  if (maxdepth < 1) stop('maxdepth must be positive')

  dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F);
  classes <- unique(dep_values); where <- tapply(weights, dep_values, sum, simplify = T)[classes]
  if (is.null(misclassification_weights)) {
    misclassification_weights <- matrix(1, nrow = length(classes), ncol = length(classes)) - diag(length(classes))
    rownames(misclassification_weights) <- classes; colnames(misclassification_weights) <- classes
  } else {
    if (dim(misclassification_weights)[1] != dim(misclassification_weights)[2]) stop('misclassification_weights must be square')
    if (any(!(classes %in% colnames(misclassification_weights))) | any(!(classes %in% rownames(misclassification_weights)))) {
      stop('column and rownames of misclassification_weights must contain all possible classes')
    }
  }
  # Define loss function
  which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights))); W <- misclassification_weights[which_W, which_W];

  Lhat <- min(where %*% t(W))

  lossfunction <- function(object.dtree){
    Z <- find_pos_points(object.dtree, data); leafs_subtree <- (length(object.dtree@bvec) + 1):(length(object.dtree@bvec) * 2 + 1)
    where <- t(apply(Z[, leafs_subtree, drop = F], 2, function(x) {
                                                        index <- which(x == 1);
                                                        tab <- tapply(weights[index], dep_values[index], sum, simplify = T)[classes]
                                                        tab[is.na(tab)] <- 0
                                                        tab
                                                      }))

    weighted_misclass <- sum(apply(where %*% t(W), 1, min));
    loss <- 1/Lhat * weighted_misclass + alpha * complexity_tree_fun(object.dtree)
  }

  # Define function that decides with nodes are eligible for consideration
  work_with_this_node <- function(node, object.dtree){
    if (node == 1) T else if (any(object.dtree@A[, find_parent_node(node)] == 1)) T else F
  }
  # Start while-loop to consider all nodes. The loop ends when every node is visited but nothing is changed.
  loopFinished <- F
  while (!loopFinished){
    # Assess which nodes are allowed to change, such that the structure of the tree is not nonsensical and
    # the maximum depth is not exceeded.
    #show(object.dtree)
    numbernodes <- (length(object.dtree@bvec) * 2) + 1
    l <- sapply(1:numbernodes, work_with_this_node, object.dtree = object.dtree)
    endbranch <- (numbernodes - 1) / 2; beginleaf <- (numbernodes + 1) / 2; endleaf <- numbernodes
    if (log(beginleaf)/log(2) == maxdepth) nodes <- (1:endbranch)[l[1:endbranch]] else nodes <- (1:endleaf)[l]

    # Find the position of all points
    Z <- find_pos_points(object.dtree, data)

    # Shuffle the nodes that can be considered in random order and loop through them.
    shufflednodes <- sample(nodes, size = length(nodes), replace = F)
    for (node in shufflednodes){
      #cat(node, '\n')
      # Find the subtree rooted at the considered node and the subset of data points in the subtree
      subtree <- find_subtree(rootnode = node, object.dtree); object.subdtree <- subtree$object.subdtree
      indexset <- which(Z[,node] == 1)
      datasub <- data[indexset, ]; weightssub <- weights[indexset]
      # Optimize node parallel
      best_subtree <- optimize_node_parallel(subtree, datasub, weightssub, data, weights, alpha, minleafsize, misclassification_weights)
      # Check if the subtree has changed
      if (best_subtree$newtree){
        # If here, the subtree needs to be replaced
        object.bestsubdtree <- best_subtree$object.bestsubdtree
        if (subtree$leafnode){
          # If here rootnode was previously a leaf node and now becomes a branch node
          depthNode <- floor(log(node)/log(2))

          object.dtree@A <- cbind(object.dtree@A, matrix(0, nrow = dim(object.dtree@A)[1], ncol = 2^depthNode))
          object.dtree@bvec <- c(object.dtree@bvec, rep(0, 2^depthNode))
          object.dtree@cats <- c(object.dtree@cats, vector('list', 2^depthNode))

          object.dtree@A[, node] <- object.bestsubdtree@A
          object.dtree@bvec[node] <- object.bestsubdtree@bvec
          object.dtree@cats[node] <- object.bestsubdtree@cats
        } else if (length(object.subdtree@bvec) == length(object.bestsubdtree@bvec)){
          # If here rootnode was previously a branch node and stays that, while subtree does not change size
          object.dtree@A[, subtree$branchnodes] <- object.bestsubdtree@A
          object.dtree@bvec[subtree$branchnodes] <- object.bestsubdtree@bvec
          object.dtree@cats[subtree$branchnodes] <- object.bestsubdtree@cats
        } else {
          # If here rootnode was previously a branch node and stays that, while subtree becomes smaller
          # Thus the rootnode is replaced by the lower or upper subtree.
          number_oldcols <- length(object.subdtree@bvec); number_newcols <- length(object.bestsubdtree@bvec)

          object.dtree@A[, subtree$branchnodes[1:number_newcols]] <- object.bestsubdtree@A
          object.dtree@bvec[subtree$branchnodes[1:number_newcols]] <- object.bestsubdtree@bvec
          object.dtree@cats[subtree$branchnodes[1:number_newcols]] <- object.bestsubdtree@cats

          object.dtree@A[, subtree$branchnodes[(number_newcols + 1):number_oldcols]] <- 0
          object.dtree@bvec[subtree$branchnodes[(number_newcols + 1):number_oldcols]] <- 0
          object.dtree@cats[subtree$branchnodes[(number_newcols + 1):number_oldcols]] <- vector('list', number_oldcols - number_newcols)
        }
        # Check if there is an unneccessary depth in the tree and if so, adjust.
        object.dtree <- cut_tree_fun(object.dtree)
        # Restart the loop
        break
      } else {
        # If here, the subtree stays the same.
        # Now check if we have assessed all nodes without a change
        if (node == tail(shufflednodes, 1)){
          # If here the local search procedure is done and the outer while-loop must be exited
          loopFinished <- T
        }

      }
    }
  }
  # Find the classification and add to the dtree object
  analysis <- find_where_and_leafs(object.dtree, data, weights, misclassification_weights)
  object.dtree@leaf_classes <- analysis$leaf_classes

  # Find the loss associated with the tree
  loss <- lossfunction(object.dtree)
  return(list(object.dtree = object.dtree, loss = loss))
}

#' Perform Optimize Node Parallel procedure
#'
#' The Optimize Node Parallel procedure checks to see if a change in the split (or leaf) at the rootnode
#' of a subtree is necessary to reduce the loss associated to the subtree. If the rootnode is a branchnode,
#' it is checked whether 1) there exists a better split, 2) it is better to replace with the lower subtree, or
#' 3) it is better to replace with upper subtree. If the rootnode is a leaf, it is checked wheter it is better
#' to apply a split.
#'
#'
#' @param subtree A list with an object of class dtree, the associated branchnodes
#' (result of function find_subtree) and a logical indicating if the subtree is a leaf node.
#' @param datasub  A dataframe: all observations that reach the rootnode of the subtree
#' @param data A dataframe: the full dataset
#' @param alpha A complexity parameter
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @return A list of two, where the first element is an object of class dtree and the second
#' element is a logical value indicating if the subtree has changed.
optimize_node_parallel <- function(subtree, datasub, weightssub, data, weights, alpha, minleafsize, misclassification_weights){

  # Define loss function
  dep_values_full <- unlist(select(data, subtree$object.subdtree@dep_var), use.names = F);
  classes_full <- unique(dep_values_full); where_full <- tapply(weights, dep_values_full, sum, simplify = T)[classes_full]
  which_W_full <- sapply(classes_full, function(x) which(x == rownames(misclassification_weights))); W_full <- misclassification_weights[which_W_full, which_W_full];
  Lhat <- min(where_full %*% t(W_full))

  dep_values_sub <- unlist(select(datasub, subtree$object.subdtree@dep_var), use.names = F); classes_sub <- unique(dep_values_sub)
  which_W_sub <- sapply(classes_sub, function(x) which(x == rownames(misclassification_weights))); W_sub <- misclassification_weights[which_W_sub, which_W_sub];

  lossfunction <- function(object.dtree){

    if (length(classes_sub) == 1){return(alpha * complexity_tree_fun(object.dtree))}

    Z <- find_pos_points(object.dtree, datasub); leafs_subtree <- (length(object.dtree@bvec) + 1):(length(object.dtree@bvec) * 2 + 1)

    where_sub <- t(apply(Z[, leafs_subtree, drop = F], 2, function(x) {
                                                        index <- which(x==1);
                                                        tab <- tapply(weightssub[index], dep_values_sub[index], sum, simplify = T)[classes_sub]
                                                        tab[is.na(tab)] <- 0
                                                        tab
                                                      }))

    weighted_misclass <- sum(apply(where_sub %*% t(W_sub), 1, min));
    loss <- 1/Lhat * weighted_misclass + alpha * complexity_tree_fun(object.dtree)
    return(loss)
  }

  # Initialize variables
  newtree <- F
  object.subdtree <- subtree$object.subdtree; object.bestsubdtree <- object.subdtree
  bestloss <- lossfunction(object.bestsubdtree)

  if (subtree$leafnode){
    # If here, then the root node is a leaf node.
    # Check if it is better to remain a leaf node or to apply a split

    # Find best split (if possible)
    best_parallel_tree <- find_best_parallel_split(object.subdtree, datasub, weightssub, minleafsize, misclassification_weights)

    if (best_parallel_tree$newtree){
      # If here, a split is possible
      loss <- lossfunction(best_parallel_tree$object.bestsubdtree)
      if (loss < bestloss){
        # If here, the loss associated with the new split is less than before
        # The leaf node becomes a branch node

        object.bestsubdtree <- best_parallel_tree$object.bestsubdtree
        newtree <- T
      }
    }
    # Return the best tree

    return(list(object.bestsubdtree = object.bestsubdtree, newtree = newtree))

  } else {
    # If here, the rootnode was a branchnode

    # 1) Find best split (if possible)
    best_parallel_tree <- find_best_parallel_split(object.subdtree, datasub, weightssub, minleafsize, misclassification_weights)


    if (best_parallel_tree$newtree){
      # If here, the found split is different from the current split

      loss <- lossfunction(best_parallel_tree$object.bestsubdtree)
      if (loss < bestloss){
        # If here, the loss associated with the new split is less than before
        object.bestsubdtree <- best_parallel_tree$object.bestsubdtree
        bestloss <- loss
        newtree <- T
      }
    }

    # 2) Replace subtree with lower subtree
    lower_subtree <- find_subtree(rootnode = 2, object.subdtree)
    object.lowersubdtree <- lower_subtree$object.subdtree

    if (subtree$branchnodes[1] != 1 | !all(object.lowersubdtree@A[, 1] == 0)){
      # If here, there are two possibilities: i) the rootnode is the rootnode of the whole tree (that is,
      # the subtree is the full tree) but the lower subtree is not a leaf node, or ii) the rootnode is not
      # the rootnode of the full tree. In either case, it may be replaces by the lower subtree.

      loss <- lossfunction(object.lowersubdtree)
      if (loss < bestloss){
        # If here, the subtree should be replaced by the lower subtree

        object.bestsubdtree <- object.lowersubdtree
        bestloss <- loss
        newtree <- T
      }
    }

    upper_subtree <- find_subtree(rootnode = 3, object.subdtree)
    object.uppersubdtree <- upper_subtree$object.subdtree
    if (subtree$branchnodes[1] != 1| !all(object.uppersubdtree@A[, 1] == 0)){
      # If here, there are two possibilities: i) the rootnode is the rootnode of the whole tree (that is,
      # the subtree is the full tree) but the upper subtree is not a leaf node, or ii) the rootnode is not
      # the rootnode of the full tree. In either case, it may be replaces by the upper subtree.

      loss <- lossfunction(object.uppersubdtree)
      if (loss < bestloss){
        # If here, the subtree should be replaced by the upper subtree

        object.bestsubdtree <- object.uppersubdtree
        bestloss <- loss
        newtree <- T
      }
    }

  }
  # Return the best tree
  return(list(object.bestsubdtree = object.bestsubdtree, newtree = newtree))
}

#' Find best parallel split
#'
#' For a specific branchnode (root of inputted tree), find the optimal split.
#'
#' @param object.subdtree Object of class dtree. The split of the rootnode will be optimized.
#' @param datasub  A dataframe: all observations that reach the rootnode of the subtree
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @return A list of two, where the first element is an object of class dtree and the second
#' element is a logical value indicating if the subtree has changed.
find_best_parallel_split <- function(object.subdtree, datasub, weightssub, minleafsize, misclassification_weights){
  # Initialize variables
  features <- object.subdtree@features; dep_values <- unlist(select(datasub, object.subdtree@dep_var), use.names = F)
  classes <- unique(dep_values);

  if (length(classes) == 1) {return(list(object.bestsubdtree = object.subdtree, newtree = F))}

  which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights)))
  W <- misclassification_weights[which_W, which_W]

  leafs_subtree <- (length(object.subdtree@bvec) + 1):(length(object.subdtree@bvec) * 2 + 1);
  Z <- find_pos_points(object.subdtree, datasub);
  where <- t(apply(Z[, leafs_subtree, drop = F], 2, function(x) { index <- which(x==1);
                                                          tab <- tapply(weightssub[index], dep_values[index], sum, simplify = T)[classes]
                                                          tab[is.na(tab)] <- 0
                                                          tab
                                                         }))

  bestmisclass <- sum(apply(where %*% t(W), 1, min));

  object.bestsubdtree <- object.subdtree
  newtree <- F

  # Loop through without parallel
  for (feature in features){
      #cat('At feature ', feature, ' with ', length(dep_values), ' observations and ', length(classes), ' classes \n')
      #tic()
      res <- find_best_feature_split(feature, object.subdtree, datasub, weightssub, minleafsize, misclassification_weights)
      #toc()

      if (res$splitpossible){
    # If here, there is a split possible for this feature
        if (res$misclass < bestmisclass){
          # If here, it is better than the current split
          # Replace the subtree
          object.bestsubdtree <- res$object.bestsubdtree
          bestmisclass <- res$misclass
          newtree <- T
        }
      }
  }
  return(list(object.bestsubdtree = object.bestsubdtree, newtree = newtree))
}


#' Find best split value for certain feature with weighted observations
#'
#' For a specific branchnode (root of inputted tree) and feature, find the optimal split.
#'
#' @param feature Character name of feature that will be split on.
#' @param object.subdtree Object of class dtree. The split of the rootnode will be optimized.
#' @param datasub  A dataframe: all observations that reach the rootnode of the subtree
#' @param minleafsize The minimum number of observations in any terminal <leaf> node.
#' @return A list of two, where object.bestsubdtree is an object of class dtree, misclass is the total
#' number of misclassified observations, and splitpossible is a logical value that indicated whether a
#' split on this feature is possible.
find_best_feature_split <- function(feature, object.subdtree, datasub, weightssub, minleafsize, misclassification_weights){
  # Initialize variables, clear all entries associated with rootnode of object.subdtree
  featindex <- which(feature == object.subdtree@features)
  object.bestsubdtree <- object.subdtree

  object.bestsubdtree@A[, 1] <- 0;  object.bestsubdtree@A[featindex, 1] <- 1
  object.bestsubdtree@bvec[1] <- 0; object.bestsubdtree@cats[1] <- vector('list', 1)

  unique_sorted_values <- sort(unique(unlist(select(datasub, feature), use.names = F)))
  dep_values <- unlist(select(datasub, object.subdtree@dep_var), use.names = F)
  m <- length(unique_sorted_values)
  n <- length(unique(dep_values))
  bestmisclass <- Inf; splitpossible <- F

  # If m equals 1, there cannot be a split
  if (m > 1 & n > 1) {
    if (!is.factor(unique_sorted_values)){
      # If here the feature is numeric
      values <- unlist(select(datasub, feature), use.names = F)
      ord <- order(values); ordered_values <- values[ord]; ordered_dep_values <- dep_values[ord]
      ordered_datasub <- datasub[ord, ]; ordered_table <- table(ordered_values)
      ordered_weightssub <- weightssub[ord]

      unique_ordered_values <- unique(ordered_values);
      intermediate_values <- sapply(1:(length(unique_ordered_values) - 1),
                                    function(i) (unique_ordered_values[i] + unique_ordered_values[i + 1]) / 2)

      classes <- unique(dep_values); which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights)))
      W <- misclassification_weights[which_W, which_W]


      # Find the upper and lower subtrees and the position of all points if they would follow either the
      # lower or upper tree
      if (length(object.subdtree@bvec) > 1) {
        leafs_subtrees <- ((length(object.subdtree@bvec) + 1) / 2):length(object.subdtree@bvec);
        number_leafs_subtrees <- length(leafs_subtrees)

        object.lower <- find_subtree(2, object.subdtree)$object.subdtree
        object.upper <- find_subtree(3, object.subdtree)$object.subdtree

        Zleft <- find_pos_points(object.lower, ordered_datasub)
        Zright <- find_pos_points(object.upper, ordered_datasub)
      } else {
        leafs_subtrees <- 1; number_leafs_subtrees <- 1

        Zleft <- matrix(1, nrow = nrow(datasub), ncol = 1)
        Zright <- matrix(1, nrow = nrow(datasub), ncol = 1)
      }


      # Find the leafs that should be taken into account for the calculation of the minimum leaf size
      shouldbesomething <- c(colSums(Zleft)[leafs_subtrees] > 0, colSums(Zright)[leafs_subtrees] > 0)

      where <- matrix(0, nrow = number_leafs_subtrees * 2, ncol = length(classes)); colnames(where) <- classes
      where[(number_leafs_subtrees + 1):(number_leafs_subtrees * 2), ] <- t(apply(Zright[, leafs_subtrees, drop = F], 2,
                                                                                  function(x) { index <- which(x==1);
                                                                                  tab <- tapply(ordered_weightssub[index], ordered_dep_values[index],sum, simplify = T)[classes]
                                                                                  tab[is.na(tab)] <- 0
                                                                                  tab
                                                                                  }))

      posZleft <- Zleft; posZright <- Zright; pos_values <- ordered_values; pos_dep_values <- ordered_dep_values; pos_weights <- ordered_weightssub
      for (i in 1:length(intermediate_values)){
        int_value <- intermediate_values[i]
        indices <- which(pos_values < int_value)

        where[1:number_leafs_subtrees, ] <- where[1:number_leafs_subtrees, ] +
          t(apply(posZleft[indices, leafs_subtrees, drop = F], 2,
                  function(x) { index <- which(x==1);
                  tab <- tapply(pos_weights[index], pos_dep_values[index], sum, simplify = T)[classes]; tab[is.na(tab)] <- 0
                  tab
                  }))

        where[(number_leafs_subtrees + 1):(number_leafs_subtrees * 2), ] <- where[(number_leafs_subtrees + 1):(number_leafs_subtrees * 2), ] -
          t(apply(posZright[indices, leafs_subtrees, drop = F], 2,
                  function(x) { index <- which(x==1);
                  tab <- tapply(pos_weights[index], pos_dep_values[index], sum, simplify = T)[classes]; tab[is.na(tab)] <- 0
                  tab
                  }))

        # Count the misclassification
        weighted_misclass <- sum(apply(where %*% t(W), 1, min))

        # Count the minimum leaf size
        minleafsizelocal <- min(rowSums(where)[shouldbesomething])
        if (weighted_misclass < bestmisclass & minleafsizelocal >= minleafsize){
          # If here then the misclassification is better than the current and the minimum leaf size is greater
          # than the prescribed minimum
          best_b <- int_value
          object.bestsubdtree@bvec[1] <- best_b
          bestmisclass <- weighted_misclass
          splitpossible <- T
        }
        posZleft <- posZleft[-indices, , drop = F]; posZright <- posZright[-indices, , drop = F]; pos_values <- pos_values[-indices];
        pos_dep_values <- pos_dep_values[-indices]; pos_weights <-pos_weights[-indices]
      }
    } else {
      # If here, the feature is a categorical feature and we must consider combinations of different levels
      # as possible splits

      # Initialize variables
      values <- unlist(select(datasub, feature), use.names = F)
      ord <- order(values); ordered_values <- values[ord]; ordered_dep_values <- dep_values[ord]
      ordered_datasub <- datasub[ord, ]; ordered_table <- table(droplevels(ordered_values))
      ordered_weightssub <- weightssub[ord]

      levelsfeat <- names(ordered_table); numberlevels <- length(levelsfeat)

      classes <- unique(dep_values); which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights)))
      W <- misclassification_weights[which_W, which_W]

      leafs_subtrees <- ((length(object.subdtree@bvec) + 1) / 2):length(object.subdtree@bvec);
      number_leafs_subtrees <- length(leafs_subtrees)

      # Find the upper and lower subtrees and the position of all points if they would follow either the
      # lower or upper tree
      object.lower <- find_subtree(2, object.subdtree)$object.subdtree
      object.upper <- find_subtree(3, object.subdtree)$object.subdtree

      Zleft <- find_pos_points(object.lower, ordered_datasub)
      Zright <- find_pos_points(object.upper, ordered_datasub)

      # Find the leafs that should be taken into account for the calculation of the minimum leaf size
      shouldbesomething <- c(colSums(Zleft)[leafs_subtrees] > 0, colSums(Zright)[leafs_subtrees] > 0)

      # where_list is a list of length 'number of levels in the feature' and described for each level
      # where the observation of that value would be if they would either follow the lower or the
      # upper branch
      where_list <- vector('list', length(ordered_table))

      listindex <- 1
      index <- 1
      for (lev in levelsfeat){
        where <- matrix(0, nrow = length(leafs_subtrees) * 2, ncol = length(classes))
        colnames(where) <- classes
        indexset <- index:(index + ordered_table[lev] - 1)

        where[1:number_leafs_subtrees,] <-
          t(apply(Zleft[indexset, leafs_subtrees, drop = F], 2,
                  function(x) {index <- which(x==1);
                  tab <- tapply(ordered_weightssub[index], ordered_dep_values[index], sum, simplify = T)[classes]; tab[is.na(tab)] <- 0
                  tab
                  }))
        where[(1 + number_leafs_subtrees):(2 * number_leafs_subtrees), ] <-
          t(apply(Zright[indexset, leafs_subtrees, drop = F], 2,
                  function(x) {index <- which(x==1);
                  tab <- tapply(ordered_weightssub[index], ordered_dep_values[index], sum, simplify = T)[classes]; tab[is.na(tab)] <- 0
                  tab
                  }))


        where_list[[listindex]] <- where
        index <- index + ordered_table[lev]
        listindex <- listindex + 1
      }

      # Loop over all possible combinations with outer and inner loop
      outeroptions <- floor(numberlevels/2)
      for (n in 1:outeroptions){
        options <- combn(1:length(ordered_table), n); numberoptions <- dim(options)[2]
        for (i in 1:numberoptions){
          list_indices <- options[, i]

          bigwhere <- matrix(0, nrow = length(leafs_subtrees) * 2, ncol = length(classes))
          for (i in 1:length(where_list)){
            if (is.element(i, list_indices)){
              bigwhere[1:number_leafs_subtrees, ] <- bigwhere[1:number_leafs_subtrees, ] + where_list[[i]][1:number_leafs_subtrees, ]
            } else {
              bigwhere[(1 + number_leafs_subtrees):(number_leafs_subtrees * 2), ] <-
                bigwhere[(1 + number_leafs_subtrees):(number_leafs_subtrees * 2), ] +
                where_list[[i]][(1 + number_leafs_subtrees):(number_leafs_subtrees * 2), ]
            }
          }
          # Count the misclassification
          weighted_misclass <- sum(apply(bigwhere %*% t(W), 1, min))


          # Count the minimum leaf size
          minleafsizelocal <- min(rowSums(bigwhere)[shouldbesomething])
          if (weighted_misclass < bestmisclass & minleafsizelocal >= minleafsize){
            # If here then the misclassification is better than the current and the minimum leaf size is greater
            # than the prescribed minimum

            object.bestsubdtree@cats[[1]] <- names(ordered_table)[list_indices]
            bestmisclass <- weighted_misclass
            splitpossible <- T
          }
        }
      }
    }
  }
  if (splitpossible){
    return(list(object.bestsubdtree = object.bestsubdtree, misclass = bestmisclass, splitpossible = T))
  } else {
    return(list(splitpossible = F))
  }
}

