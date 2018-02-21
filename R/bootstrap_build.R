rand.two.graphs <- function(edges1, edges2) {
  clust1 <- as.factor(edges1[upper.tri(edges1)])
  clust2 <- as.factor(edges2[upper.tri(edges2)])
  rand <- sum(clust1 == clust2) / length(clust1)
  return(rand)
}

train.oneboot <- function(expr, nkeep1, args) {
  # sample at random
  sel.ind <- sample(1:nrow(expr), nrow(expr), replace=TRUE)
  not.selected <- setdiff(1:nrow(expr), sel.ind)
  oob.data <- expr[not.selected, ]
  
  # train
  args$expr <- expr[sel.ind,]
  oneboot.res <- do.call("build.network", args)
  kept.ind <- which(oneboot.res$n.edges > nkeep1)
  if (length(kept.ind)>0) {
    sel.ind <- min(kept.ind, na.rm=TRUE)
  } else {
    sel.ind <- length(oneboot.res$network)
    warning(paste0("'nkeep1' not reached. I advice you decrease its value.\n
  Current 'nkeep1' is ",nkeep1))
  }
  kept.res <- (oneboot.res$network[[sel.ind]] != 0)
  mse <- 0
  mse_mu <- 0
  first_list <- oneboot.res$options$first.list
  second_list <- oneboot.res$options$second.list
  for (ind in 1:ncol(args$expr)) {
    mse <- mse + mean((oob.data[,ind,drop=FALSE] - 
                         as.matrix(oob.data[,-ind]) %*% 
                         oneboot.res$betas[[sel.ind]][,ind,drop=FALSE])^2)
    mse_mu <- mse_mu + mean((oob.data[,ind,drop=FALSE] - 
                               as.matrix(oob.data[,-ind]) %*% 
                               oneboot.res$betas[[sel.ind]][,ind,drop=FALSE])^2)
    if (!is.null(first_list) & 
          (colnames(expr)[ind] %in% union(first_list[,1], first_list[,2]))) {
      where_r <- first_list[,1] == colnames(expr)[ind]
      if (sum(where_r) != 0)
        mse_mu <- mse_mu + oneboot.res$options$mu^2/2 * 
        (oneboot.res$betas[[sel.ind]][match(first_list[where_r,2],colnames(expr)[-ind]),ind] - 1)^2
      where_r <- first_list[,2] == colnames(expr)[ind]
      if (sum(where_r) != 0)
        mse_mu <- mse_mu + oneboot.res$options$mu^2/2 * 
        (oneboot.res$betas[[sel.ind]][match(first_list[where_r,1],colnames(expr)[-ind]),ind] - 1)^2
    }
    
    if (!is.null(second_list) & 
          (colnames(expr)[ind] %in% union(second_list[,1], second_list[,2]))) {
      where_r <- second_list[,1] == colnames(expr)[ind]
      if (sum(where_r) != 0)
        mse_mu <- mse_mu + oneboot.res$options$mu^2 *
        (oneboot.res$betas[[sel.ind]][match(second_list[where_r,2],colnames(expr)[-ind]),ind])^2
      where_r <- second_list[,2] == colnames(expr)[ind]
      if (sum(where_r) != 0)
        mse_mu <- mse_mu + oneboot.res$options$mu^2 * 
        (oneboot.res$betas[[sel.ind]][match(second_list[where_r,1],colnames(expr)[-ind]),ind])^2
    }
  }
  
  return(list("edges" = kept.res, "mse" = mse, "mse_mu" = mse_mu))
}

bootstrap.build <- function(expr, nboot=10, nkeep1=round(0.1*nrow(expr)^2/2),
                            nkeep2=round(0.5*nboot), parallel=FALSE, ncores=10,
                            seeds = NULL, ...) {
  
  args <- list(...)
  p <- ncol(expr)
  if (is.null(args$verbose)) args$verbose <- FALSE
  
  if (parallel) {
    registerDoMC(cores=ncores)
    
    res <- foreach(rep=1:nboot) %dopar% {
      if (!is.null(seeds)) set.seed(seeds[rep])
      stability <- 0
      cat("bootstrap sample ", rep, "\n")
      kept.res <- train.oneboot(expr, nkeep1, args)
      return(kept.res)
    }
    cat("calculate stability...\n")
    stability <- 0
    for (rep in 2:nboot) {
      for (rep2 in 1:(rep-1)) {
        stability <- stability + 
          rand.two.graphs(res[[rep]]$edges, res[[rep2]]$edges)
      }
    }
  } else {
    res <- list()
    stability <- 0
    for (rep in 1:nboot) {
      cat("bootstrap sample ", rep, "\n")
      res[[rep]] <- train.oneboot(expr, nkeep1, args)
      if (rep > 1) {
        for (rep2 in 1:(rep-1)) {
          stability <- stability + 
            rand.two.graphs(res[[rep]]$edges, res[[rep2]]$edges)
        }
      }
    }
    
  }
  stability <- stability / (nboot * (nboot -1) / 2)
  total.mse <- mean(unlist(lapply(res, function(alist) alist$mse)))
  edge.res <- lapply(res, function(alist) alist$edges)
  total.res <- Reduce("+", edge.res)
  instability <- total.res[upper.tri(total.res)] / nboot
  instability <- 2 * mean(instability) / (p * (p-1))
  final.res <- total.res > nkeep2
  
  
  return(list("edges" = final.res, "mse" = total.mse, 'adjacency'=total.res,
              "stars" = instability, "stability" = stability))
}

# data(cancer)
# # no constraint
# res <- bootstrap.build(expr, nboot=3, nkeep2=1)
# # constraints
# res <- bootstrap.build(expr, nboot=3, nkeep2=1,
#                        first.list=matrix(c("AMFR","BB_S4","BECNI","BTG3"),
#                                          ncol=2, byrow=TRUE),
#                        second.list=matrix(c("AMFR","E2F3"), ncol=2),
#                        mu=1, r.mu=2)
# print(res)