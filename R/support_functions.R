#' Find complexity of a dtree object
#'
#' @param object.dtree Object of class dtree
#' @return Numeric value that is the complexity (number of branchnodes)
complexity_tree_fun <- function(object.dtree){
  sum(sum(object.dtree@A))
}

#' Reduce depth of dtree if possible
#'
#' If there is a whole depth layer of 'branchnodes' that contains no branch node, this dimension is removed
#' @param object.dtree Object of class dtree
#' @return Object of class dtree
cut_tree_fun <- function(object.dtree){
  A <- object.dtree@A; bvec <- object.dtree@bvec; cats <- object.dtree@cats
  if (all(A == 0)) return(object.dtree)
  number_branch_nodes <- length(bvec)
  branch_nodes <- ((number_branch_nodes + 1)/2):number_branch_nodes
  if (all(bvec[branch_nodes] == 0) & all(compare.list(cats[branch_nodes], vector('list', length(branch_nodes))))){
    A <- A[, -branch_nodes, drop = F]
    bvec <- bvec[-branch_nodes]
    cats <- cats[-branch_nodes]
  }
  object.dtree@A <- A; object.dtree@bvec <- bvec; object.dtree@cats <- cats
  return(object.dtree)
}

#' Basic analysis of performance of dtree ob dataset
#'
#' Finds the classes associated with leaf nodes, the misclassification per leaf node and the minimum
#' leaf size.
#'
#' @param object.dtree Object of class dtree
#' @param data A dataframe with features and dependent variable class.
#' @return A list with leaf_classes an factor array of classes prescribed by each leaf node,
#' misclass a numeric array of the misclassification on the data associated with each leaf node, and
#' minleafsize the minimum leaf size in a leaf node.
find_class_leafs_and_misclass_and_min_leaf_size_fun <- function(object.dtree, data){
  if (length(object.dtree@bvec) == 0){
    dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F)
    cl <- get_class_node(dep_values)
    misclass <-  sum(dep_values != cl)
    return(list(leaf_classes = cl,
                misclass = misclass,
                minleafsize = length(dep_values)))
  }
  Z <- find_pos_points(object.dtree, data)
  dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F)
  number_nodes <- dim(Z)[2]
  leaf_nodes <- ((number_nodes + 1) / 2):number_nodes
  leaf_classes <- rep(NA, length(leaf_nodes))
  misclass <- rep(0, length(leaf_nodes))
  index <- 1
  for (node in leaf_nodes){
    indexset <- which(Z[, node] == 1)
    if (length(indexset) != 0){
      node_classes <- dep_values[indexset]
      leaf_classes[index] <- as.character(get_class_node(node_classes))
      misclass[index] <- sum(node_classes != leaf_classes[index])
    }
    index <- index + 1
  }
  work_with_this_node <- function(node, object.dtree){
    if (node == 1) {
      T
    } else {
      A <- object.dtree@A
      pnode <- find_parent_node(node)
      if (any(A[,pnode] == 1)) {T} else {F}
    }
  }
  indexset <- sapply(X = 1:(dim(Z)[2]), FUN = work_with_this_node, object.dtree = object.dtree)
  minleafsize <- min(colSums(Z)[indexset])
  return(list(leaf_classes = leaf_classes, misclass = misclass, minleafsize = minleafsize))
}

#' Basic analysis of performance of dtree ob dataset
#'
#' Finds the classes associated with leaf nodes, the misclassification per leaf node and the minimum
#' leaf size.
#'
#' @param object.dtree Object of class dtree
#' @param data A dataframe with features and dependent variable class.
#' @return A list with leaf_classes an factor array of classes prescribed by each leaf node,
#' misclass a numeric array of the misclassification on the data associated with each leaf node, and
#' minleafsize the minimum leaf size in a leaf node.
find_where_and_leafs <- function(object.dtree, data, weights, misclassification_weights) {
  dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F); classes <- unique(dep_values)
  leafs <- (length(object.dtree@bvec) + 1):(length(object.dtree@bvec) * 2 + 1)
  Z <- find_pos_points(object.dtree, data)

  where <- t(apply(Z[, leafs, drop = F], 2, function(x) {
    index <- which(x == 1);
    tab <- tapply(weights[index], dep_values[index], sum, simplify = T)[classes]
    tab[is.na(tab)] <- 0
    tab}))
  if (length(classes) == 1) where <- t(where)
  colnames(where) <- classes

  which_W <- sapply(classes, function(x) which(x == rownames(misclassification_weights))); W <- misclassification_weights[which_W, which_W];
  leaf_classes <- apply(where %*% t(W), 1, function(x) {
    if (all(x == 0)) {NA} else {as.character(classes[which(x == min(x))[1]])}
  })
  weighted_misclass <- sum(apply(where %*% t(W), 1, min))
  return(list(leaf_classes = leaf_classes, where = where, weighted_misclass = weighted_misclass))
}

#' Find the classification from a set of observations
#'
#' @param obs A numeric array with classes
#' @return The most prevalent class
get_class_node <- function(obs) {
  uniqobs <- unique(obs)
  as.character(uniqobs[which.max(tabulate(match(obs, uniqobs)))])
}



#' Find position of all data points in tree
#'
#' @param object.dtree Object of class dtree
#' @param data A dataframe with features and dependent variable class.
#' @return A matrix in which each row indicated a observation and each column indicates a node. A 1 in
#' position (i, j) indicated that observation i is in node j.
find_pos_points <- function(object.dtree, data){
  number_branch_nodes <- length(object.dtree@bvec)
  number_nodes <- number_branch_nodes * 2 + 1
  number_observations <- dim(data)[1]
  Z <- matrix(0, nrow = number_observations, ncol = number_nodes)
  Z[, 1] <- 1
  for (branchnode in 1:number_branch_nodes){
    dnodes <- find_daughter_nodes(branchnode)

    indexset <- which(Z[, branchnode] == 1)
    if (length(indexset) != 0){
      feat <- which(object.dtree@A[, branchnode] == 1); obs <- unlist(select(data[indexset, ], object.dtree@features[feat]), use.names = F)
      if (object.dtree@bvec[branchnode] != 0){
        newz <- obs < object.dtree@bvec[branchnode]
        Z[indexset, dnodes[1]] <- as.numeric(newz)
        Z[indexset, dnodes[2]] <- as.numeric(!newz)

      } else if (!is.null(object.dtree@cats[[branchnode]])){
        newz <- is.element(obs, object.dtree@cats[[branchnode]])
        Z[indexset, dnodes[1]] <- as.numeric(newz)
        Z[indexset, dnodes[2]] <- as.numeric(!newz)

      } else {
        Z[indexset, dnodes[1]] <- 0
        Z[indexset, dnodes[2]] <- 1
      }
    } else {
      Z[, dnodes] <- 0
    }
  }
  return(Z)
}


#' Find parent nodes
#'
#' @param node Numeric value
#' @return A numeric array of all parent of node
find_parent_nodes <- function(node){
  if (node == 1){
    return('Error: root_node')
  }
  current_parents <- find_parent_node(node)

  if (findAll){
    while (tail(current_parents, n = 1) > 1){
      next_parent <- find_parent_node(tail(current_parents, n = 1))
      current_parents <- c(current_parents, next_parent)
    }
  }
  return(current_parents)
}

#' Find parent node
#'
#' @param node Numeric value
#' @return A numeric value indicating the direct parent of node
find_parent_node <- function(node){
  if (node%%2 == 0){
    pnode <- node/2
  } else {
    pnode <- (node - 1)/2
  }
  return(pnode)
}

#' Find daughter nodes
#'
#' @param node Numeric value
#' @return A numeric array of 2 indicating the direct daughters of node
find_daughter_nodes <- function(node){
  depth <- floor(log(node)/log(2))
  side <- node - 2^depth
  dnodes <- c(2 ^ (depth + 1) + side * 2, 2 ^ (depth + 1) + side * 2 + 1)
  return(dnodes)
}

#' Find subtree rooted at node
#'
#' @param rootnode Rootnode of subtree
#' @param object.dtree Object of class dtree
#' @return A list with object.subdtree being an object of class dtree and the subtree rooted at the rootnode
#' inputted, branchnodes is the a numeric array with all branchnodes of the subtree, and leafnode is a logical
#' indicating wheter the rootnode is a leaf node
find_subtree <- function(rootnode, object.dtree){
  number_branchnodes <- dim(object.dtree@A)[2]
  depth_tree <- log(number_branchnodes + 1)/log(2)
  depth_node <- floor(log(rootnode)/log(2))
  depth_subtree <- depth_tree - depth_node
  if (depth_subtree == 0){
    leafnode <- T; branchnodes_subtree <- numeric()

    A_subtree <- matrix(0, nrow = dim(object.dtree@A)[1], ncol = 1)
    bvec_subtree <- 0
    cats_subtree <- vector('list', 1)
  } else {
    leafnode <- F
    number_branchnodes_subtree <- 2 ^ (depth_subtree) - 1
    number_nodes_subtree <- 2 ^ (depth_subtree + 1) - 1
    nodes_subtree <- rootnode
    for (index in 1:number_branchnodes_subtree){
      daughters <- find_daughter_nodes(nodes_subtree[index])
      nodes_subtree[(index*2):(index*2+1)] <- daughters
    }
    branchnodes_subtree <- nodes_subtree[1:number_branchnodes_subtree]
    A_subtree <- object.dtree@A[, branchnodes_subtree, drop = F]
    bvec_subtree <- object.dtree@bvec[branchnodes_subtree]
    cats_subtree <- object.dtree@cats[branchnodes_subtree]
  }
  subtree <- new("dtree",
                 A = A_subtree,
                 bvec = bvec_subtree,
                 cats = cats_subtree,
                 dep_var = object.dtree@dep_var,
                 leaf_classes = character(),
                 features = object.dtree@features)
  res <- list(object.subdtree = subtree, branchnodes = branchnodes_subtree, leafnode = leafnode)
  return(res)
}

#' Prepare dataset for local search procedure
#'
#' This function normalized all numeric features from 0 to 1 and forces all character features to be factors,
#' as well as forcing the dependent variable (class) to be a factor
#'
#' @param data A dataframe
#' @param weights Weights
#' @param features Character array of the features that will be used in training the OCT
#' @param dep_var Character that indicated the dependent variable (class name)
#' @return A list with data being a dataframe with normalized numeric features, and all other variables factors,
#' and a matrix with the minimum and maximum values of all numeric features.
#' @export
#' @importFrom dplyr select
prepare_dataset <- function(data, weights, features, dep_var){
  norm_feat <- function(value, minvec, maxvec){
    (value - minvec) / (maxvec - minvec)
  }
  minmax <- matrix(0, nrow = 2, ncol = length(features))

  smalldata <- dplyr::select(data, c(features, dep_var))
  na_cases <- which(!complete.cases(data))
  if (length(na_cases) > 0) {data <- smalldata[-na_cases, ]; weights <- weights[-na_cases]}
  for (feature in features){
    index <- which(names(data) == feature)
    if (is.character(data[, index])){
      data[, index] <- as.factor(data[, index])
    } else if (!is.factor(data[, index])){
      vec <- data[, index]; minvec <- min(vec); maxvec <- max(vec)
      minmax[, index] <- rbind(minvec, maxvec)
      data[, index] <- sapply(data[, index], norm_feat, minvec = minvec, maxvec = maxvec)
    }
  }
  index <- which(names(data) == dep_var)
  data[, index] <- as.factor(data[, index])
  return(list(data = data, weights = weights, minmax = minmax))
}

