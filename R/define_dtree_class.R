# define class tree and methods -----


#' An S4 class to represent a decision tree.
#'
#' @slot A A numeric matrix with p rows and q columns, where each column represents a branch node.
#' Each column contains at most one '1', which corresponds to the feature on which the split is made.
#' @slot bvec A length q numeric vector where each entry represent a branch node and the split value
#' (if the split is on a numeric feature).
#' @slot cats A list of length q, where each list entry represents a branch and is a character array
#' (if the split is on a categorical feature).
#' @slot dep_var A character string for the dependent variable
#' @slot leaf_classes The predicted class for point in certain leaf node
#' @slot features A character vector of all p feature on which the tree is built.
#' @export dtree
dtree <- setClass("dtree",
         slots=list(A = "matrix",
                    bvec = "numeric",
                    cats = "list",
                    dep_var = "character",
                    leaf_classes = "character",
                    features = "character"))

# define class tree and methods -----


#' Displays object of 'dtree' class
#' @docType methods
#' @rdname dtree-method
#' @export

setMethod("show",
          "dtree",
          function(object) {
            if (length(object@bvec) == 0 | (object@bvec[1] == 0 & is.null(object@cats[[1]]))){
              cat('No branch nodes')
            } else {
              numbernodes <- length(object@bvec) * 2 + 1
              rules <- rep(NA, numbernodes)
              for (node in 1:length(object@bvec)){
                if (any(object@A[, node] == 1)){
                  depth <- floor(log(node)/log(2))
                  indent <- paste(rep('\t', depth), collapse = " ")
                  param <- object@features[which(object@A[, node] == 1)]
                  if (object@bvec[node] != 0){
                    split <- round(object@bvec[node], 2)
                    rule <- paste(param, '<', split)
                  } else {
                    split <- object@cats[[node]]
                    split <- paste(split, collapse = " or ")
                    rule <- paste(param, '=', split)
                  }
                  rules[node] <- paste0(indent, ' Node ', node, ' : ', rule, '\n')
                }
              }
              leafs <- (length(object@bvec) + 1):((length(object@bvec) + 1) * 2 - 1)
              if (length(object@leaf_classes) == length(leafs)){
                indent <- paste(rep('\t', depth + 1), collapse = " ")
                index <- 1
                for (node in leafs) {
                  if (!is.na(object@leaf_classes[index])){
                    rules[node] <- paste0(indent, ' Leaf: ', object@dep_var, ' = ', object@leaf_classes[index], '\n')
                  }
                  index <- index + 1
                }
              }
              depthtree <- log(numbernodes + 1) / log(2) - 1
              ord <- make_order(depthtree)
              for (lin in ord){
                if (!is.na(rules[lin])){
                  cat(rules[lin])
                }
              }
              cat('Features are ', paste(object@features), '\n', 'Dependent variable is ', object@dep_var)
              visualize_tree(object)
            }
          }
)
#' Predict Method for dtree Fits
#'
# #' @param object An object of class dtree
# #' @param newdata A data frame
# #' @return An array of predicted classes for each observationin the newdata data frame.
#' @export
setMethod('predict',
          signature(object="dtree"),
          function(object, newdata, ...){
            if (missing(newdata)){
              return('Error: needs new data')
            }
            leafs <- (length(object@bvec) + 1):(length(object@bvec) * 2 + 1)
            leaf_classes <- object@leaf_classes
            Z <- find_pos_points(object, newdata)
            y <- sapply(1:dim(Z)[1], function(rowindex){
              whichleaf <- which(Z[rowindex, leafs] == 1)
              yel <- leaf_classes[whichleaf]
            })
            as.factor(y)
          })



make_order <- function(depthtree){
  ord <- 1
  depth <- 1
  while (depth <= depthtree + 1){
    for (node in ord){
      depthnode <- floor(log(node)/log(2))
      if (depthnode == depth - 1){
        dnodes <- find_daughter_nodes(node)
        index <- which(node == ord)
        ord <- append(ord, dnodes, index)
      }
    }
    depth <- depth + 1
  }
  return(ord)
}


# visualize_tree -----


#' Visualize decision tree
#'
#' \code{visualize_tree} plots the decision tree including splits and
#' leaf assignments in the plot window and prints splits and leafs in command window.
#'
#' @param object.tree an object of class \code{dtree}
#'
#' @examples A <- matrix(0, nrow = 4, ncol = 3); A[3,1] <- 1; A[4,3] <- 1
#' tree1 <- new("dtree", A = A, bvec = c(2.5, 0, 1.8), cats = vector('list', 3),
#'              dep_var = 'Species', leaf_classes = c(NA, 'setosa', 'versicolor', 'virginica'),
#'              features = c("Sepal.Length", "Sepal.Width",  "Petal.Length", "Petal.Width"))
#' visualize_tree(tree1)
#' @importFrom diagram coordinates
#' @importFrom diagram openplotmat
#' @importFrom diagram bentarrow
#' @importFrom diagram textrect
#' @importFrom diagram textempty
#' @importFrom diagram textellipse
visualize_tree <- function(object.dtree){
  A <- object.dtree@A; bvec <- object.dtree@bvec; cats <- object.dtree@cats;
  leaf_classes <- object.dtree@leaf_classes; features <- object.dtree@features
  numberbranchnodes <- dim(A)[2]
  depthtree <- log(numberbranchnodes + 1)/log(2)
  numberleafnodes <- 2^(depthtree + 1) - 1 - numberbranchnodes
  elpos <- coordinates (2 ^ (0:depthtree))

  coorarrows <- double()
  for (node in 1:numberbranchnodes){
    daughters <- find_daughter_nodes(node)
    if (bvec[node] != 0 | !is.null(cats[[node]])){ # not simple branch node
      if (daughters[1] <= numberbranchnodes){ # node connected to branch node
        for (dnode in daughters){
          if (bvec[dnode] != 0 | !is.null(cats[[dnode]]) ){
            coorarrows <- c(coorarrows, c(node, dnode))
          } else {
            depthnode <- floor(log(dnode)/log(2))
            for (i in 1:(depthtree-depthnode)){
              dnode <- find_daughter_nodes(dnode)[2]
            }
            coorarrows <- c(coorarrows, c(node, dnode))
          }
        }
      } else {
        coorarrows <- c(coorarrows, c(node, daughters[1], node, daughters[2]))
      }
    }
  }
  fromto <- matrix(ncol = 2, byrow = TRUE,
                   data = coorarrows)

  nr     <-nrow(fromto)

  arrpos <- matrix(ncol = 2, nrow = nr)

  titl <- paste0("Decision Tree")
  openplotmat(main = titl)

  for (i in 1:nr) {
    arrpos[i, ] <-  bentarrow (
      to = elpos[fromto[i, 2], ],
      from = elpos[fromto[i, 1], ],
      lwd = 2, arr.pos = 0.6,
      arr.length = 0
    )
  }

  for (node in 1:(numberbranchnodes)){
    if (bvec[node] != 0 | !is.null(cats[[node]])){
      midnode <- elpos[node,]
      decision_parameter <- features[which(A[, node] == 1)]
      if (bvec[node] != 0){
        decision_value <- bvec[node]
        textrect   (midnode ,radx = .055, rady = .04, cex = 0.8, shadow.size = 0,
                    lab = paste(decision_parameter, "\n", '<', round(decision_value, 2)))
      } else {
        decision_value <- cats[[node]]
        textrect   (midnode, radx = .055, rady = .04, cex = 0.8, shadow.size = 0,
                    lab = paste(decision_parameter, "\n", 'is', paste(decision_value, collapse = " or ")))
      }
    }
  }
  index <- 1
  for (node in (numberbranchnodes + 1):(numberbranchnodes + numberleafnodes)){
    if (!is.na(leaf_classes[index])){
      class_node <- leaf_classes[index]
      textellipse   (elpos[node,], radx = .025, rady = .025, cex = 0.8, shadow.size = 0,
                  lab = paste(class_node))
    }
    index <- index + 1
  }
}

#' Visualize decision tree with data
#'
#' \code{visualize_tree} plots the decision tree including splits and
#' leaf assignments in the plot window and prints splits and leafs in command window.
#'
#' @param object.tree an object of class \code{dtree}
#' @param data dataframe
#' @param weights weights
#' @export
#' @importFrom dplyr select
#' @importFrom diagram coordinates
#' @importFrom diagram openplotmat
#' @importFrom diagram bentarrow
#' @importFrom diagram textrect
#' @importFrom diagram textempty
#' @importFrom diagram textellipse
visualize_tree_data <- function(object.dtree, data, weights){
  A <- object.dtree@A; bvec <- object.dtree@bvec; cats <- object.dtree@cats;
  leaf_classes <- object.dtree@leaf_classes; features <- object.dtree@features
  numberbranchnodes <- dim(A)[2]
  depthtree <- log(numberbranchnodes + 1)/log(2)
  numberleafnodes <- 2^(depthtree + 1) - 1 - numberbranchnodes

  dep_values <- unlist(select(data, object.dtree@dep_var), use.names = F); classes <- unique(dep_values)
  nodes <- 1:(2^(depthtree + 1) - 1)
  Z <- find_pos_points(object.dtree, data)

  where <- t(apply(Z[, nodes, drop = F], 2, function(x) {
    index <- which(x == 1);
    tab <- tapply(weights[index], dep_values[index], sum, simplify = T)[classes]
    tab[is.na(tab)] <- 0
    tab}))
  if (length(classes) == 1) where <- t(where)
  colnames(where) <- classes

  sumtotal <- sum(weights)

  titl <- paste0("Decision Tree")
  openplotmat(main = titl)

  textsize <- 1 - (depthtree ** 1.2)*0.1

  rec_width <- strwidth(paste0('(', paste0(round(where[1,]/sum(where[1,]), 2), collapse = ", "), ')'), cex = textsize)*1.5
  rec_height <- strheight(paste0('line1 \n line2'), cex = textsize)*1.5
  ellip_width <- strwidth(paste0('(', paste0(round(where[1,]/sum(where[1,]), 2), collapse = ", "), ')'), cex = textsize)*1.1
  ellip_height <- strheight(paste0('line1 \n line2 \n line3'), cex = textsize)*1.5

  wraptext <- function(rule){
    wrappedrule <- rule
    i <- 1
    while (i < nchar(wrappedrule)){
      substr_rule <- substr(wrappedrule, 1, i)
      if (strwidth(substr_rule, cex = textsize, font = 2) > rec_width){
        wrappedrule <- paste0(prevsubstr_rule, '\n', substr(wrappedrule, i, nchar(wrappedrule)))
        i <- i + 1
      }
      i <- i + 1
      prevsubstr_rule <- substr_rule
    }
    return(wrappedrule)
  }



  elpos <- coordinates (2 ^ (0:depthtree))

  coorarrows <- double()
  for (node in 1:numberbranchnodes){
    daughters <- find_daughter_nodes(node)
    if (bvec[node] != 0 | !is.null(cats[[node]])){ # not simple branch node
      if (daughters[1] <= numberbranchnodes){ # node connected to branch node
        for (dnode in daughters){
          if (bvec[dnode] != 0 | !is.null(cats[[dnode]]) ){
            coorarrows <- c(coorarrows, c(node, dnode))
          } else {
            depthnode <- floor(log(dnode)/log(2))
            for (i in 1:(depthtree-depthnode)){
              dnode <- find_daughter_nodes(dnode)[2]
            }
            coorarrows <- c(coorarrows, c(node, dnode))
          }
        }
      } else {
        coorarrows <- c(coorarrows, c(node, daughters[1], node, daughters[2]))
      }
    }
  }
  fromto <- matrix(ncol = 2, byrow = TRUE,
                   data = coorarrows)

  nr     <- nrow(fromto)

  arrpos <- matrix(ncol = 2, nrow = nr)

  for (i in 1:nr) {
    arrpos[i, ] <-  bentarrow (
      to = elpos[fromto[i, 2], ],
      from = elpos[fromto[i, 1], ],
      lwd = 2, arr.pos = 0.6,
      arr.length = 0
    )
  }

  for (node in 1:(numberbranchnodes)){
    if (bvec[node] != 0 | !is.null(cats[[node]])){
      midnode <- elpos[node,]
      decision_parameter <- features[which(A[, node] == 1)]
      in_this_node <- where[node, ]
      prop <- paste0('(', paste0(round(in_this_node/sum(in_this_node), 2), collapse = ", "), ')')
      tot <- paste0(round(sum(in_this_node)/sumtotal*100, 2), '%')
      textrect (midnode ,radx = rec_width/2, rady = rec_height/2, cex = textsize, shadow.size = 0,
                lab = paste(prop, "\n", tot))
      if (bvec[node] != 0){
        decision_value <- bvec[node]
        rule <- paste(decision_parameter, '<', round(decision_value, 2))
        wrapped_rule <- wraptext(rule)
        textempty (midnode - c(0, rec_height/2 + strheight(wrapped_rule, cex = textsize, font = 2)/2)*1.1, cex = textsize, font = 2, lab = wrapped_rule)
        #textrect   (midnode ,radx = .055, rady = .04, cex = textsize, shadow.size = 0,
        #            lab = paste(decision_parameter, "\n", '<', round(decision_value, 2)))
      } else {
        decision_value <- cats[[node]]
        rule <- paste(decision_parameter, ' is ', paste(decision_value, collapse = " or "))
        wrapped_rule <- wraptext(rule)
        textempty (midnode - c(0, rec_height/2 + strheight(wrapped_rule, cex = textsize, font = 2)/2)*1.1, cex = textsize, font = 2, lab = wrapped_rule)
        #            lab = paste(decision_parameter, "\n", 'is', paste(decision_value, collapse = " or ")))
      }
    }
  }
  index <- 1
  for (node in (numberbranchnodes + 1):(numberbranchnodes + numberleafnodes)){
    if (!is.na(leaf_classes[index])){
      class_node <- leaf_classes[index]
      in_this_node <- where[node, ]
      if (depthtree > 3 | (depthtree > 2 & ncol(where) > 2)){
        if (node %% 2 == 0){
          midnode <- elpos[node, ] + c(0, ellip_height/4)
        } else {
          midnode <- elpos[node, ] - c(0, ellip_height/4)
        }
      } else {
        midnode <- elpos[node, ]
      }

      prop <- paste0('(', paste0(round(in_this_node/sum(in_this_node), 2), collapse = ", "), ')')
      tot <- paste0(round(sum(in_this_node)/sumtotal*100, 2), '%')
      textellipse (midnode ,radx = ellip_width/2, rady = ellip_height/2, cex = textsize, shadow.size = 0,
                lab = paste(class_node, "\n", prop, "\n", tot))

      #textellipse   (elpos[node,], radx = .025, rady = .025, cex = textsize, shadow.size = 0,
      #               lab = paste(class_node))
    }
    index <- index + 1
  }
}
