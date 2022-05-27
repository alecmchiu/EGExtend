
demean <- function(mat,class_df){
  mat2 <- mat
  for (each in unique(class_df[,2])){
        ids <- subset(class_df,class_df[,2] == each)[,1]
        mat2[ids,] <- scale(mat2[ids,],scale=F)
  }
  return(mat2)
}

generate_all_interaction <- function(X_mat) {
  dat <- X_mat
  colnames(dat) <- 1:ncol(dat)
  do.call(cbind, combn(colnames(dat), 2, FUN= function(x)
    list(setNames(data.frame(dat[,x[1]]*dat[,x[2]]),
                  paste(x, collapse="_")) )))
}

generate_squared_terms <- function(g){
  return(apply(g,2,function(x){x^2}))
}


