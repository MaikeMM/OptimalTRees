#' # NOT USED
#' find_best_feature_split2 <- function(feature, object.subdtree, datasub, minleafsize){
#'   featindex <- which(feature == object.subdtree@features)
#'   object.subdtree@A[, 1] <- 0;  object.subdtree@A[featindex, 1] <- 1
#'   object.subdtree@bvec[1] <- 0; object.subdtree@cats[1] <- vector('list', 1)
#'
#'   unique_sorted_values <- sort(unique(unlist(select(datasub, feature), use.names = F)))
#'   m <- length(unique_sorted_values)
#'
#'   bestmisclass <- Inf; splitpossible <- F
#'   if (m > 1) {
#'     if (!is.factor(unique_sorted_values)){
#'       while (m > 20){
#'         localbestmisclass <- Inf
#'         splitvalues <- split(unique_sorted_values, ceiling(seq_along(unique_sorted_values) / (m / 10)))
#'         for (index in 1:length(splitvalues)){
#'           values <- splitvalues[[index]]; b <- (values[1] + tail(values, 1))/2
#'           object.subdtree@bvec[1] <- b
#'           analysis <- find_class_leafs_and_misclass_and_min_leaf_size_fun(object.subdtree, datasub)
#'           misclass <- sum(analysis$misclass); minleafsizelocal <- analysis$minleafsize
#'           if (misclass < localbestmisclass & minleafsizelocal >= minleafsize){
#'             bestindex <- index
#'             localbestmisclass <- misclass
#'           }
#'         }
#'         if (is.infinite(localbestmisclass)) return(list(splitpossible = F))
#'         if (bestindex == 1) {
#'           unique_sorted_values <- c(splitvalues[[1]], splitvalues[[2]])
#'         } else if (bestindex == length(splitvalues)) {
#'           unique_sorted_values <- c(splitvalues[[bestindex - 1]], splitvalues[[bestindex]])
#'         } else {unique_sorted_values <- c(splitvalues[[bestindex - 1]], splitvalues[[bestindex]], splitvalues[[bestindex + 1]])}
#'         m <- length(unique_sorted_values)
#'       }
#'       for (i in 1:(m-1)){
#'         b <- (unique_sorted_values[i] + unique_sorted_values[i + 1]) / 2
#'         object.subdtree@bvec[1] <- b
#'         analysis <- find_class_leafs_and_misclass_and_min_leaf_size_fun(object.subdtree, datasub)
#'         misclass <- sum(analysis$misclass); minleafsizelocal <- analysis$minleafsize
#'         if (misclass < bestmisclass & minleafsizelocal >= minleafsize){
#'           object.bestsubdtree <- object.subdtree
#'           bestmisclass <- misclass
#'           splitpossible <- T
#'         }
#'       }
#'     } else {
#'       levelsfeat <- levels(unique_sorted_values); numberlevels <- length(levelsfeat); ceil <- floor(numberlevels/2)
#'       for (n in 1:ceil){
#'         options <- combn(levelsfeat, n); numberoptions <- dim(options)[2]
#'         for (i in 1:numberoptions){
#'           object.subdtree@cats[[1]] <- options[, i]
#'           analysis <- find_class_leafs_and_misclass_and_min_leaf_size_fun(object.subdtree, datasub)
#'           misclass <- sum(analysis$misclass); minleafsizelocal <- analysis$minleafsize
#'           if (misclass < bestmisclass & minleafsizelocal >= minleafsize){
#'             object.bestsubdtree <- object.subdtree
#'             bestmisclass <- misclass
#'             splitpossible <- T
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   if (splitpossible){
#'     return(list(object.bestsubdtree = object.bestsubdtree, misclass = bestmisclass, splitpossible = T))
#'   } else {
#'     return(list(splitpossible = F))
#'   }
#' }
#'
#'
#' #for (feature in features){
#' #  res <- find_best_feature_split(feature, object.subdtree, datasub, minleafsize)
#' #  if (res$splitpossible){
#'     # If here, there is a split possible for this feature
#' #    if (res$misclass < bestmisclass){
#'       # If here, it is better than the current split
#'       # Replace the subtree
#' #      object.bestsubdtree <- res$object.bestsubdtree
#' #      bestmisclass <- res$misclass
#' #      newtree <- T
#' #    }
#' #  }
#' #}
#'
#' #' Find position of all data points in tree
#' #'
#' #' @param object.dtree Object of class dtree
#' #' @param data A dataframe with features and dependent variable class.
#' #' @return A matrix in which each row indicated a observation and each column indicates a node. A 1 in
#' #' position (i, j) indicated that observation i is in node j.
#' find_pos_points <- function(object.dtree, data){
#'   number_branch_nodes <- dim(object.dtree@A)[2]
#'   depth_tree <- log(number_branch_nodes + 1)/log(2)
#'   number_observations <- dim(data)[1]
#'   number_nodes <- 2^(depth_tree+1) - 1
#'   Z <- matrix(0, nrow = number_observations, ncol = number_nodes)
#'   for (node in 1:number_nodes){
#'     if (node == 1){
#'       Z[, node] <- 1
#'     } else {
#'       newz <- points_satisfy_split(node, object.dtree, data)
#'       pnode <- find_parent_node(node)
#'       Z[, node] <- newz * Z[, pnode]
#'     }
#'   }
#'   return(Z)
#' }
#'
#' #' Find which points satisfy a split
#' #'
#' #' @param node Numeric value, must be branchnode in object.dtree
#' #' @param object.dtree Object of class dtree
#' #' @param data A dataframe with features and dependent variable class.
#' #' @return An numeric array in which each value indicated whether an observations satisfies the split at
#' #' the branch node (1) or not (0)
#' points_satisfy_split <-function(node, object.dtree, data){
#'   pnode <- find_parent_node(node)
#'   avec <- object.dtree@A[, pnode];
#'
#'   if (all(avec == 0)){
#'     if (node%%2 == 0){
#'       z <- rep.int(x = 0, times = dim(data)[1])
#'       return(z)
#'     } else {
#'       z <- rep.int(x = 1, times = dim(data)[1])
#'       return(z)
#'     }
#'   } else {
#'     feat <- which(avec == 1)
#'     b <- object.dtree@bvec[pnode]; cats <- object.dtree@cats
#'     if (b == 0){
#'       split <- cats[[pnode]]
#'       if (node%%2 == 0){  #yes branch of parent node was followed
#'         applySplit <- function(x) {if (is.element(x, split)) {z <- 1} else {z <- 0}}
#'       } else { #no branch of parent node was followed
#'         applySplit <- function(x) {if (!is.element(x, split)) {z <- 1} else {z <- 0}}
#'       }
#'     } else {
#'       if (node%%2 == 0){ #yes branch of parent node was followed
#'         applySplit <- function(x) {if (x < b) {z <- 1} else {z <- 0}}
#'       } else { #no branch of parent node was followed
#'         applySplit <- function(x) {if (x >= b) {z <- 1} else {z <- 0}}
#'       }
#'     }
#'   }
#'   z <- apply(select(data, object.dtree@features[feat]), 1, applySplit)
#'   return(z)
#' }

