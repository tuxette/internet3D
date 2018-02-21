##### This file is the main file of therese
# it contains the function that allows you to infer networks
# it is widely inspired from the code of Julien Chiquet's R package 'simone'

build.network <- function(expr, ...) {
  # Initialize options
  args <- list(...)
  options <- do.call("set.options", args)
  
  # Centering (always) and normalizing if required (TRUE by default)
  X <- scale.default(expr, TRUE, options$scale)
  
  # Setting unprovided parameters
  if (is.null(options$max.it)) 
    options$max.it <- 10 * min(c(nrow(expr), ncol(expr)))
  r.max <- max(apply(X,2,var))
  if (is.null(options$mu)) options$mu <- 0.5 * r.max
  if (is.null(options$penalty.max)) options$penalty.max <- r.max
  options$penalty.max <- min(options$penalty.max, r.max)
  if (is.null(options$penalties)) {
    options$penalties <- seq(options$penalty.max, options$penalty.min,
                             length=options$n.penalties)
    keep.all <- FALSE
  } else keep.all <- TRUE
  options$penalties[options$penalties < options$penalty.min] <-
    options$penalty.min
  options$penalties[options$penalties > options$penalty.max] <-
    options$penalty.max
  
  # use 'first.list' and 'second.list' to set up the constraint matrix
  constraint.mat <- matrix(0, nrow=ncol(expr), ncol=ncol(expr))
  if (!is.null(options$first.list)) {
    where.first <- match(options$first.list[,1], colnames(expr))
    where.second <- match(options$first.list[,2], colnames(expr))
    all.signs <- apply(cbind(where.first, where.second), 1,
                       function(arow) sign(cor(expr[ ,arow[1]],expr[,arow[2]])))
    constraint.mat[cbind(where.first,where.second)] <- all.signs * options$mu^2
    constraint.mat[cbind(where.second,where.first)] <- all.signs * options$mu^2
  }
  if (!is.null(options$second.list)) {
    where.first <- match(options$second.list[,1], colnames(expr))
    where.second <- match(options$second.list[,2], colnames(expr))
    all.signs <- apply(cbind(where.first, where.second), 1,
                       function(arow) sign(cor(expr[ ,arow[1]],expr[,arow[2]])))
    constraint.mat[cbind(where.first,where.second)] <- Inf
    constraint.mat[cbind(where.second,where.first)] <- Inf
  }
  if (sum(constraint.mat!=0)==0) constraint.mat <- NULL
  
	if (options$verbose) {
    cat("\nNetwork Inference for Internet3D... \n")
    print(options)
    cat("\n\n")
	}
  
	networks <- list()
  betas    <- list()
  pcors    <- list()
	BIC      <- c()
	AIC      <- c()
	lambdas  <- c()
	n.edges  <- c()
	loglik   <- c()
	loglik.pen <- c()
	last.edges <- -1
	last.crit  <- -Inf
	if (options$verbose) cat(format(c("|  penalty","|    edges","| BIC"),
                                  width=10, justify="right"),"\n\n")

  # Loop on penalties...
	for (lambda in options$penalties) {
    ## The main inference function is called here...
		out <- infer.edges(X, constraint.mat, lambda, options)
		if (sum(is.na(out$n.edges))) break 
		if (sum(out$n.edges  - last.edges) < 0) break
		if (!keep.all) {
			if (all(out$n.edges == last.edges)) {
				lambdas[length(lambdas)] <- lambda
				next
			}
		}

		# Gather the results together
		options$initial.guess <- out$Beta
		BIC  <- c(BIC,out$BIC)
		AIC  <- c(AIC,out$AIC)
		loglik  <- c(loglik,out$loglik)
		loglik.pen  <- c(loglik.pen,out$loglik.pen)
		lambdas <- c(lambdas,lambda)
		last.edges <- out$n.edges
		n.edges    <- rbind(n.edges,out$n.edges)
		last.crit  <- out$loglik.pen
		networks[[length(networks)+1]] <- out$Theta
    pcors[[length(networks)]] <- out$partial.cor
    betas[[length(networks)]] <- out$Beta
    
		if (options$verbose) {
			cat(format(list(lambda,paste(last.edges,collapse=","),out$BIC), width=10,
                 digits=4, justify="right"),"\n")
		}
		if (max(last.edges) > options$edges.max) break
	}

	# Return results
	res <- structure(list(data=expr, networks=networks, betas=betas, pcors=pcors, 
                        used.lambdas=as.vector(lambdas), 
                        n.edges=as.matrix(n.edges), BIC=as.vector(BIC),
	                      AIC=as.vector(AIC), loglik=as.vector(loglik),
	                      loglik.pen=as.vector(loglik.pen), options=options),
	                 class="internet3D")
  if (options$verbose) {cat("\n\n"); summary(res)}
  return(res)
}

##### S3 method for 'internet3D'
print.internet3D <- function(x,...) {
  cat("Internet3D object inferred for", ncol(x$data), "variables observed for",
      nrow(x$data), "individuals.\n\n")
}

summary.internet3D <- function(object,...) {
  cat("Internet3D object inferred for", ncol(object$data),
      "variables observed for", nrow(object$data), "individuals.\n\n")
  p <- ncol(object$data)
  
  print(object$options)

  cat("\n***** Results obtained by Internet3D\n\n")
  cat("    Best BIC for network number", which.min(object$BIC), "with",
      object$n.edges[which.min(object$BIC),],"edges (densities:",
      format(2*object$n.edges[which.min(object$BIC),]/p/(p-1),digits=2), ").\n")
  cat("    Best AIC for network number", which.min(object$AIC), "with",
      object$n.edges[which.min(object$AIC),],"edges (densities:",
      format(2*object$n.edges[which.min(object$AIC),]/p/(p-1),digits=2), ").\n")
  cat("    Best penalized log-likelihood for network number",
      which.max(object$loglik.pen), "with",
      object$n.edges[which.max(object$loglik.pen),], "edges. (densities:",
      format(2*object$n.edges[which.max(object$loglik.pen),]/p/(p-1), digits=2),
      ").\n")
}

## Tests
# data(cancer)
# # no constraint
# res <- build.network(expr)
# # constraints
# res <- build.network(expr, 
#                      first.list=matrix(c("AMFR","BTG3","BECNI","BTG3"),
#                                        ncol=2, byrow=TRUE),
#                      second.list=matrix(c("AMFR","E2F3"), ncol=2),
#                      mu=1, r.mu=2)
# print(res)