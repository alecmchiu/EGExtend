
#' Remove class-specific mean from each class in matrix
#'
#' Remove the class-specific mean for each feature/column.
#'
#' @param mat A matrix with observations (e.g. individuals) as rows and features as columns.
#' @param class_df A data frame with two columns. The first columns are observation identifiers and the second column denotes class/group.
#'
#' @return A matrix with the class-specific mean removed from each observation in its appropriate group.
#' @export
#'
#' @examples
#'
#' # Generate 200 observation by 10 features matrix
#' # 100 samples of class one with features means with 1-10
#' # 100 samples of class two with feature means 10-1
#' obs <- rbind(
#'   sapply(1:10, function(x) {
#'     rnorm(100, mean = x)
#'   }),
#'   sapply(10:1, function(x) {
#'     rnorm(100, mean = x)
#'   })
#' )
#'
#' # Check feature means for each class
#' colMeans(obs[1:100, ])
#' colMeans(obs[101:200, ])
#'
#' # Create class data frame
#' classes <- cbind(1:200, c(rep(0, 100), rep(1, 100)))
#' colnames(classes) <- c("ID", "Class")
#'
#' # Remove class-specific means from each group
#' demean_obs <- demean(obs, classes)
#'
#' # Check if means are near 0 for each class
#' colMeans(demean_obs[1:100, ])
#' colMeans(demean_obs[101:200, ])
demean <- function(mat, class_df) {
  mat2 <- mat
  for (each in unique(class_df[, 2])) {
    ids <- subset(class_df, class_df[, 2] == each)[, 1]
    mat2[ids, ] <- scale(mat2[ids, ], scale = FALSE)
  }
  return(mat2)
}

#' Generate pairwise interaction matrix
#'
#' This function generates a matrix where the columns/features are the pairwise
#' products of the input matrix features. The new matrix will retain the same
#' number of rows/observations, but have \emph{m} choose 2 columns for each pair
#' of features, where \emph{m} is the number of features in the input matrix.
#' When using this function where the mean has been removed, the feature/column
#' means become observations of covariance.
#'
#' @param mat A matrix with observations (e.g. individuals) as rows and features
#'   as columns.
#' @param delimiter A delimiter for the pairwise feature names.
#'
#' @return A matrix with all pairwise interactions (pairwise column/feature
#'   products) from the input matrix. Column names will show what features were
#'   used for the pariwise product.
#' @export
#'
#' @examples
#' # Generate covariance matrix
#' library(mvtnorm)
#'
#' cov_mat <- rWishart(1, 10, diag(10))[, , 1]
#'
#' # Generate correlated data with means 1-10
#' dat <- rmvnorm(100, mean = 1:10, sigma = cov_mat)
#' colnames(dat) <- LETTERS[1:10]
#'
#' # Remove mean from matrix
#' centered_dat <- scale(dat, scale = FALSE)
#' colMeans(centered_dat)
#'
#' # Generate pairwise matrix
#' product_mat <- generate_all_interaction(centered_dat)
#'
#' # Check if covariance estimates are close
#' colMeans(product_mat)
#' cov(dat)[1:7, 1:7]
generate_all_interaction <- function(mat, delimiter = "_") {
  dat <- mat
  do.call(cbind, utils::combn(colnames(dat), 2, FUN = function(x) {
    list(stats::setNames(
      data.frame(dat[, x[1]] * dat[, x[2]]),
      paste(x, collapse = delimiter)
    ))
  }))
}

#' Generate squared matrix
#'
#' This function generates a matrix where the columns/features of the input
#' matrix features have been squared. The new matrix will retain the dimensions
#' of the input matrix. On a centered matrix, this function will turn the
#' observations into observations of variance.
#'
#' @param mat A matrix with observations (e.g. individuals) as rows and features
#'   as columns.
#'
#' @return A matrix with each column/feature squared.
#' @export
#'
#' @examples
#' # Generate data
#' obs <- sapply(1:10, function(x) {
#'   rnorm(10000, mean = x, sd = x^2)
#' })
#'
#' # Center matrix
#' centered_obs <- scale(obs, scale = FALSE)
#'
#' # Generate squared matrix
#' sq_obs <- generate_squared_terms(centered_obs)
#'
#' # Compare variance estimates
#' colMeans(sq_obs)
#' apply(obs, 2, var)
generate_squared_terms <- function(mat) {
  return(apply(mat, 2, function(x) {
    x^2
  }))
}

#' Run Correlation by Individual Level Product (CILP)
#'
#' Runs an efficient implementation of CILP from Lea et al. 2019. CILP detects
#' differential correlation of pairwise features between two groups.
#'
#' @param prod_mat A matrix with all class-specific mean-removed pairwise
#'   interactions to be tested.
#' @param groupings A numeric vector denoting class/group membership of the
#'   observations.
#'
#' @return A data frame of effect size estimates and p-values for each feature
#'   pair.
#' @export
#'
#' @examples
#' # Generate covariance matrix
#' library(mvtnorm)
#'
#' cov_mats <- rWishart(2, 10, diag(10))
#'
#' # Generate correlated data with means 1-10 and 10-1
#' dat <- rbind(rmvnorm(100, mean = 1:10, sigma = cov_mats[, , 1]), rmvnorm(100, mean = 10:1, sigma = cov_mats[, , 2]))
#' colnames(dat) <- LETTERS[1:10]
#'
#' # Create group vector
#' groups <- c(rep(0, 100), rep(1, 100))
#'
#' # Remove mean from matrix
#' centered_dat <- scale(dat, scale = FALSE)
#' colMeans(centered_dat)
#'
#' # Generate pairwise matrix
#' product_mat <- generate_all_interaction(centered_dat)
#'
#' # Run CILP
#' res <- CILP(product_mat, groups)
#'
CILP <- function(prod_mat, groupings) {
  pvals <- apply(prod_mat, 2, function(x) {
    ctest <- stats::cor.test(scale(x), groupings)
    return(c(ctest$estimate, ctest$p.value))
  })
  pvals <- data.frame(t(pvals))
  colnames(pvals) <- c("Effect Size", "p-value")
  return(pvals)
}

#' Calculate eigengene
#'
#' @param mat A matrix for which the eigengene is to be calculated.
#'
#' @return The eigengene of the input matrix.
#' @export
#'
#' @examples
#' # Generate data
#' obs <- rbind(
#'   sapply(1:10, function(x) {
#'     rnorm(100, mean = x)
#'   }),
#'   sapply(10:1, function(x) {
#'     rnorm(100, mean = x)
#'   })
#' )
#'
#' # Calculate Eigengene
#' eg <- calculateEigengene(obs)
#'
#' # Eigengene should reflect mean differences
#' plot(eg, ylab = "Eigengene Value")
calculateEigengene <- function(mat) {
  eigengene <- irlba::irlba(scale(mat), nv = 1)$u
  rownames(eigengene) <- rownames(mat)
  avg_expr <- rowMeans(mat)
  if (stats::cor(avg_expr, eigengene[, 1]) < 0) {
    eigengene[, 1] <- -eigengene[, 1]
  }
  return(eigengene)
}
