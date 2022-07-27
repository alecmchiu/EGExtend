
#' Remove class-specific mean from each class in matrix
#'
#' Remove the class-specific mean for each feature/column.
#'
#' @param matrix A matrix with observations (e.g. individuals) as rows and
#'   features as columns.
#' @param class_df A data frame with two columns. The first columns are
#'   observation identifiers and the second column denotes class/group.
#' @param scale A logical (TRUE/FALSE) determining whether variance should be
#'   scaled to one as well.
#'
#' @return A matrix with the class-specific mean removed from each observation
#'   in its appropriate group.
#' @seealso [base::scale()]
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
demean <- function(matrix, class_df, scale=FALSE) {
  if (is.null(ncol(class_df))){
    stop(sprintf("Class data frame does not have any columns?"))
  }
  if (ncol(class_df) < 2){
    stop(sprintf("Class data frame does not conform to specifications."))
  }
  mat2 <- matrix
  for (each in unique(class_df[, 2])) {
    ids <- subset(class_df, class_df[, 2] == each)[, 1]
    mat2[ids, ] <- scale(mat2[ids, ], scale = scale)
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
#' The pairwise interaction matrix allows one to estimate covariance when paired with the demean function. This relies on the fact that `Cov(X,Y) = E[XY]-E[X]E[Y]`, but by removing means, we eliminate the second term.
#'
#' @param matrix A matrix with observations (e.g. individuals) as rows and
#'   features as columns.
#' @param delimiter A delimiter for the pairwise feature names.
#' @param interactions_df A two-column data frame of interactions to be
#'   calculated. The first column denotes the first feature and the second
#'   column denotes the second feature of the interaction.
#'
#' @return A matrix with all pairwise interactions (pairwise column/feature
#'   products) from the input matrix. Column names will show what features were
#'   used for the pairwise product.
#' @export
#' @seealso [stats::cov()]
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
#' rownames(dat) <- 1:100
#'
#' # Remove mean from matrix
#' centered_dat <- scale(dat, scale = FALSE)
#' colMeans(centered_dat)
#'
#' # Generate pairwise matrix
#' product_mat <- generate_interactions(centered_dat)
#'
#' # Check if covariance estimates are close
#' demean_est <- colSums(product_mat)/(nrow(product_mat)-1)
#'
#' samp_cov <- sapply(strsplit(colnames(product_mat),split="_"),function(x){cov(dat)[x[1],x[2]]})
#' names(samp_cov) <- colnames(product_mat)
#'
#' head(demean_est)
#' head(samp_cov)
generate_interactions <- function(matrix, delimiter = "_", interactions_df = NA) {
  if (is.null(colnames(matrix))){
    stop(sprintf("Column names must be set."))
  }
  if (is.null(rownames(matrix))){
    stop(sprintf("Row names must be set."))
  }
  if (all(is.na(interactions_df))){
    res <- do.call(cbind, utils::combn(colnames(matrix), 2, FUN = function(x) {
      list(stats::setNames(
        data.frame(matrix[, x[1]] * matrix[, x[2]]),
        paste(x, collapse = delimiter)
      ))
    }))
  }
  else{
    res <- apply(interactions_df,1,function(x){matrix[,x[1]] * matrix[,x[2]]})
    rownames(res) <- rownames(matrix)
    colnames(res) <- apply(interactions_df,1,paste0,collapse=delimiter)
  }
  return(res)
}

#' Generate squared matrix
#'
#' This function generates a matrix where the columns/features of the input
#' matrix features have been squared. The new matrix will retain the dimensions
#' of the input matrix. On a centered matrix, this function will turn the
#' observations into observations of variance.
#' 
#' #' Using the squared matrix with the demean function allows one to estimate variance. This relies on the fact that `Var(X) = E[X^2]-E[X]^2`, but by removing means, we eliminate the second term.
#'
#' @param matrix A matrix with observations (e.g. individuals) as rows and features
#'   as columns.
#'
#' @return A matrix with each column/feature squared.
#' @export
#' @seealso [stats::var()]
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
#' sq_obs <- generate_squared_matrix(centered_obs)
#'
#' # Compare variance estimates
#' colSums(sq_obs)/(nrow(sq_obs)-1)
#' apply(obs, 2, var)
generate_squared_matrix <- function(matrix) {
  return(apply(matrix, 2, function(x) {
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
#' dat <- rbind(rmvnorm(100, mean = 1:10, sigma = cov_mats[, , 1]),
#'  rmvnorm(100, mean = 10:1, sigma = cov_mats[, , 2]))
#' colnames(dat) <- LETTERS[1:10]
#' rownames(dat) <- 1:200
#'
#' # Create group data frame
#' groups <- cbind(1:100,c(rep(0, 100), rep(1, 100)))
#'
#' # Remove mean and rescale matrix
#' rescaled_dat <- demean(dat,groups,scale = TRUE)
#' colMeans(rescaled_dat)
#'
#' # Generate pairwise matrix
#' product_mat <- generate_interactions(rescaled_dat)
#'
#' # Run CILP
#' res <- CILP(product_mat, groups[,2])
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
#' Calculates the eigengene for a given matrix. We use the \emph{irlba} package to efficiently calculate the eigengene using the Lanzcos algorithm. The eigengene is aligned to match the orientation of the input matrix and will retain the row names of the input matrix.
#'
#' The eigengene is canonically known as the first principal component. However, principal components are rotation-invariant, meaning that the result can be sometimes inverted (\emph{i.e. the results times -1}). We orient the eigengene to the input matrix by aligning the correlation between the eigengene and the vector of averages (\emph{row means}) of the input matrix.
#'
#' @param matrix A matrix for which the eigengene is to be calculated.
#'
#' @return The eigengene of the input matrix.
#' @seealso [irlba::irlba()]
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
calculateEigengene <- function(matrix) {
  eigengene <- irlba::irlba(scale(matrix), nv = 1)$u
  rownames(eigengene) <- rownames(matrix)
  avg_expr <- rowMeans(matrix)
  if (stats::cor(avg_expr, eigengene[, 1]) < 0) {
    eigengene[, 1] <- -eigengene[, 1]
  }
  return(eigengene)
}
