test_that("Demean sucessfully removes mean", {
  
  # Generate 20000 observation by 10 features matrix
  # 10000 samples of class one with features means with 1-10
  # 10000 samples of class two with feature means 10-1
  obs <- rbind(
    sapply(1:10, function(x) {
      rnorm(10000, mean = x)
    }),
    sapply(10:1, function(x) {
      rnorm(10000, mean = x)
    })
  )

  # Create class data frame
  classes <- cbind(1:20000, c(rep(0, 10000), rep(1, 10000)))
  colnames(classes) <- c("ID", "Class")

  # Remove class-specific means from each group
  demean_obs <- demean(obs, classes)

  # Check if means are near 0 for each class
  c1_means <- all(abs(colMeans(demean_obs[1:10000, ])) < 1e-8)
  c2_means <- all(abs(colMeans(demean_obs[10001:20000, ])) < 1e-8)
  expect_true(all(c(c1_means,c2_means)))
})

test_that("Demean sucessfully rescales variance", {
  
  # Generate 20000 observation by 10 features matrix
  # 10000 samples of class one with features means and variance with 1-10
  # 10000 samples of class two with feature means and variance 10-1
  obs <- rbind(
    sapply(1:10, function(x) {
      rnorm(10000, mean = x, sd = x)
    }),
    sapply(10:1, function(x) {
      rnorm(10000, mean = x, sd = x)
    })
  )
  
  # Create class data frame
  classes <- cbind(1:20000, c(rep(0, 10000), rep(1, 10000)))
  colnames(classes) <- c("ID", "Class")
  
  # Remove class-specific means from each group
  demean_obs <- demean(obs, classes, scale=T)
  
  # Check if means are near 0 for each class
  c1_means <- all(abs(colMeans(demean_obs[1:10000, ])) < 1e-8)
  c2_means <- all(abs(colMeans(demean_obs[10001:20000, ])) < 1e-8)
  
  # Check if means are near 0 for each class
  c1_var <- all((apply(demean_obs[1:10000, ],2,var) - 1) < 1e-8)
  c2_var <- all((apply(demean_obs[10001:20000, ],2,var) - 1) < 1e-8)
  expect_true(all(c(c1_means,c2_means,c1_var,c2_var)))
})

test_that("Pairwise demeaned matrix estimates covariance",{
  cov_mat <- rWishart(1, 10, diag(10))[, , 1]

  # Generate correlated data with means 1-10
  dat <- mvtnorm::rmvnorm(10000, mean = 1:10, sigma = cov_mat)
  colnames(dat) <- LETTERS[1:10]
  rownames(dat) <- 1:10000

  # Remove mean from matrix
  centered_dat <- scale(dat, scale = FALSE)

  # Generate pairwise matrix
  product_mat <- generate_interactions(centered_dat)

  # Check if covariance estimates are close
  demean_est <- colSums(product_mat)/(nrow(product_mat)-1)
  samp_cov <- sapply(strsplit(colnames(product_mat),split="_"),function(x){cov(dat)[x[1],x[2]]})
  names(samp_cov) <- colnames(product_mat)
  expect_true(all(abs(demean_est - samp_cov) < 1e-8))
})

test_that("Pairwise matrix function - column names not set",{
  cov_mat <- rWishart(1, 10, diag(10))[, , 1]
  
  # Generate correlated data with means 1-10
  dat <- mvtnorm::rmvnorm(1000, mean = 1:10, sigma = cov_mat)
  
  # Remove mean from matrix
  centered_dat <- scale(dat, scale = FALSE)
  
  # Generate pairwise matrix
  expect_error(generate_interactions(centered_dat))
})

test_that("Pairwise matrix function - row names not set",{
  cov_mat <- rWishart(1, 10, diag(10))[, , 1]
  
  # Generate correlated data with means 1-10
  dat <- mvtnorm::rmvnorm(1000, mean = 1:10, sigma = cov_mat)
  
  # Remove mean from matrix
  centered_dat <- scale(dat, scale = FALSE)
  colnames(centered_dat) <- LETTERS[1:10]
  
  # Generate pairwise matrix
  expect_error(generate_interactions(centered_dat))
})

test_that("Pairwise interactions - subset works",{
  cov_mat <- rWishart(1, 10, diag(10))[, , 1]
  
  # Generate correlated data with means 1-10
  dat <- mvtnorm::rmvnorm(10000, mean = 1:10, sigma = cov_mat)
  colnames(dat) <- LETTERS[1:10]
  rownames(dat) <- 1:10000
  
  # Create subset of interactions
  i_df <- cbind(rep('A',10),LETTERS[1:10])
  
  # Remove mean from matrix
  centered_dat <- scale(dat, scale = FALSE)
  
  # Generate pairwise matrix
  product_mat <- generate_interactions(centered_dat,interactions_df = i_df)
  
  # Check if covariance estimates are close
  demean_est <- colSums(product_mat)/(nrow(product_mat)-1)
  samp_cov <- sapply(strsplit(colnames(product_mat),split="_"),function(x){cov(dat)[x[1],x[2]]})
  names(samp_cov) <- colnames(product_mat)
  expect_true(all(abs(demean_est - samp_cov) < 1e-8))
})

test_that('Squared terms variance estimation',{
  obs <- sapply(1:10, function(x) {
    rnorm(10000, mean = x, sd = x^2)
  })

  # Center matrix
  centered_obs <- scale(obs, scale = FALSE)

  # Generate squared matrix
  sq_obs <- generate_squared_matrix(centered_obs)

  # Compare variance estimates
  expect_true(all(abs(colSums(sq_obs)/(nrow(sq_obs)-1) - apply(obs, 2, var)) < 1e-8))
})

