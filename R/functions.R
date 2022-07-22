
demean <- function(mat,class_df){
  mat2 <- mat
  for (each in unique(class_df[,2])){
        ids <- subset(class_df,class_df[,2] == each)[,1]
        mat2[ids,] <- scale(mat2[ids,],scale=F)
  }
  return(mat2)
}

generate_all_interaction <- function(mat) {
  dat <- mat
  colnames(dat) <- 1:ncol(dat)
  do.call(cbind, combn(colnames(dat), 2, FUN= function(x)
    list(setNames(data.frame(dat[,x[1]]*dat[,x[2]]),
                  paste(x, collapse="_")) )))
}

generate_squared_terms <- function(g){
  return(apply(g,2,function(x){x^2}))
}

CILP <- function(prod_mat){
    pvals <- apply(prod_mat,2,function(x){
    ctest <- cor.test(scale(x),groupings)
    return(c(ctest$estimate,ctest$p.value))
    })
}

calculateEigengene <- function(mat){
  eigengene <- irlba::irlba(scale(mat),nv=1)$u
  rownames(eigengene) <- rownames(mat)
  avg_expr <- rowMeans(mat)
  if (cor(avg_expr,eigengene[,1]) < 0){
    eigengene[,1] <- -eigengene[,1]
  }
  return(eigengene)
}