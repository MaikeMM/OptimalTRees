make_artificial_dataset_and_tree <- function(numberfeatures, numberclasses, numberobs, depthtree,
                                             densitysplits = 1, probfeatcat = 0, numberlevelscat = 4){
  isCat <- runif(numberfeatures) < probfeatcat

  features <- character()
  df <- data.frame(matrix(0, nrow = numberobs, ncol = numberfeatures + 1))
  for (i in 1:numberfeatures){
    feat <- paste0('X', as.character(i))
    features[i] <- feat
    if (isCat[i]){
      df[, i] <- as.factor(sample(LETTERS[1:numberlevelscat], numberobs, replace = T))
    } else {
      df[, i] <- runif(numberobs)
    }
  }

  numberbranchnodes <- (2^depthtree) - 1

  A <- matrix(0, nrow = numberfeatures, ncol = 1)
  bvec <- rep(0, 1)
  cats <- vector('list', 1)
  featindex <- sample(1:numberfeatures, 1)
  A[featindex, 1] <- 1
  if (isCat[featindex]){
    poslevels <- LETTERS[1:numberlevelscat]; ceil <- floor(numberlevelscat/2);
    n <- sample(1:ceil, 1); options <- combn(poslevels, n); optionindex <- sample(1:(dim(options)[2]), 1)
    split <- options[, optionindex]
    cats[[1]] <- as.character(split)
  } else {
    orderedvalues <- sort(df[, featindex])
    minsplit <- mean(head(orderedvalues, 2)); maxsplit <- mean(tail(orderedvalues, 2))
    split <- runif(1, min = minsplit, max = maxsplit)
    bvec[1] <- split
  }


  object.dtree <- new('dtree',
                      A = A,
                      bvec = bvec,
                      cats = cats,
                      dep_var = 'Y',
                      leaf_classes = character(),
                      features = features)

  for (depth in 1:(depthtree-1)){
    Z <- find_pos_points(object.dtree, df)
    branchnodes_depth <- (2 ^ depth):(2 ^ (depth + 1) - 1)
    A <- cbind(A, matrix(0, nrow = numberfeatures, ncol = length(branchnodes_depth)))
    bvec <- c(bvec, rep(0, length(branchnodes_depth)))
    cats <- c(cats, vector('list', length(branchnodes_depth)))
    for (node in branchnodes_depth){
      indexset <- which(Z[, node] == 1); dfsub <- df[indexset,]
      isbranch <- runif(1) < densitysplits
      pnode <- find_parent_node(node)
      cansplit <- any(A[, pnode] == 1)
      if (isbranch & cansplit){
        featindex <- sample(1:numberfeatures, 1); tryindex <- 1
        while (length(unique(dfsub[, featindex])) < 2){
          tryindex <- tryindex + 1; if (tryindex > numberfeatures*2) {featindex <- 0; break}
          featindex <- sample(1:numberfeatures, 1)
        }
      } else {
        featindex <- 0
      }
      if (featindex != 0){
        A[featindex, node] <- 1
        if (isCat[featindex]){
          poslevels <- unique(dfsub[, featindex]); ceil <- floor(numberlevelscat/2);
          n <- sample(1:ceil, 1); options <- combn(poslevels, n); optionindex <- sample(1:(dim(options)[2]), 1)
          split <- options[, optionindex]
          cats[[node]] <- as.character(split)
        } else {
          orderedvalues <- sort(dfsub[, featindex])
          minsplit <- mean(head(orderedvalues, 2)); maxsplit <- mean(tail(orderedvalues, 2))
          split <- runif(1, min = minsplit, max = maxsplit)
          bvec[node] <- split
        }
      }
    }
    object.dtree@A <- A; object.dtree@bvec <- bvec; object.dtree@cats <- cats
  }
  classes <- character()
  for (i in 1:numberclasses){
    cl <- paste0('c', i)
    classes[i] <- cl
  }


  Z <- find_pos_points(object.dtree, df)
  leafs <- (numberbranchnodes + 1):(numberbranchnodes * 2 + 1)
  l <- colSums(Z[, leafs]) > 0

  leaf_classes <- rep(NA, length(leafs))
  index <- 1
  for (index in 1:length(leafs)){
    if (l[index]){
      node <- leafs[index]
      cl <- sample(classes, 1)
      leaf_classes[index] <- cl
      indexset <- which(Z[,node] == 1)
      df[indexset, numberfeatures + 1] <- cl
    }
    index <- index + 1
  }
  names(df) <- c(features, 'Y')
  df$Y <- as.factor(df$Y)
  object.dtree@leaf_classes <- leaf_classes
  return(list(data = df, object.dtree = object.dtree))
}
