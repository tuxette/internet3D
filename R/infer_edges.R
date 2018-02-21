##### This file is the main file of therese

# calculate the number of edges from a matrix (loops removed)
nb.edges <- function(x) {
  edges <- sum(abs(x)>0)    
  edges <- (edges - sum(abs(diag(x))>0))/2
  edges
}

# calculate penalized log-Likelihood and its differential
LL <- function(beta, S, s, Lambda) {
  L <- 0.5*t(beta) %*% S %*% beta + t(beta) %*% s + sum(Lambda*abs(beta))
  L
}

differentialLL <- function(beta, S, s, Lambda) {
  dL <- S %*% beta + s
  zero <- beta==0
  # terms for which beta is different from zero
  if (sum(!zero)!=0)
    dL[!zero] <- dL[!zero] + Lambda[!zero] * sign(beta[!zero])
  # which ones can be null?
  if (sum(zero)!=0) {
    # terms for which beta is zero
    can.be.null <- abs(dL[zero]) < Lambda[zero]
    
    if (sum(can.be.null)!=0)
      dL[zero][can.be.null] <- 0
    
    if (sum(!can.be.null)!=0)
      dL[zero][!can.be.null] <- dL[zero][!can.be.null] -
        Lambda[zero][!can.be.null]*sign(dL[zero][!can.be.null])  
  }
  
  dL
}

# Inferrence for one variable
cLasso <- function(S, s, rho, beta, max.size, max.it) {
  # Initialization
  p <- ncol(S)
  sigma <- abs(beta) > 0
  theta <- sign(beta)
  it <- 0
  rho <- rep(rho, ncol(S))
  abs.dL.min <- abs(differentialLL(beta, S, s, rho))
  eps <- 1e-8

  if (all(abs.dL.min < eps) | sum(sigma) >= max.size) {
    return(list(beta=beta,converged=TRUE))
  } else {
    l <- which.max(abs.dL.min)
    nabla.f <- S %*% beta + s
    sigma[l] <- TRUE
    theta[l] <- -sign(nabla.f[l])
  }
  
  # Optimization
  repeat {
    it <- it+1
    sign.feasible <- FALSE
    
    while (!sign.feasible) {
      # Optimization over the active set
      h <- rep(0,p)
      x <- try(solve(S[sigma,sigma],cbind(rho*theta+s)[sigma]), silent=TRUE)
      if (is.vector(x)) {
        h[sigma] <- - beta[sigma] - x 
      } else {
        out <- optim(beta[sigma], method="BFGS", fn=LL, gr=differentialLL,
                     S=S[sigma,sigma], s=s[sigma], Lambda=rho[sigma])
        h[sigma] <- out$par-beta[sigma]
      }
      
      # Update active set
      sign.feasible <- all(sign(beta+h)[sigma]==theta[sigma])
      if (!sign.feasible) {
        ind <- which(sign(beta+h)[sigma]!=theta[sigma])
        gamma <- -beta[sigma][ind]/h[sigma][ind]
        gamma[is.nan(gamma)] <- Inf
        k <- which(gamma == min(gamma))
        if (gamma[k[1]] == Inf) {
          return(list(beta=beta,converged=FALSE))
        }
        beta <- beta + gamma[k[1]]*h
        sigma[sigma][ind][k] <- FALSE
      } else {
        beta <- beta+h
      }
      theta[sigma] <- sign(beta)[sigma]
    }
    
    abs.dL.min <- abs(differentialLL(beta,S,s,rho))
    if (all(abs.dL.min[!sigma] < eps) | it > max.it | sum(sigma) >= max.size) {
      if (it > max.it & any(abs.dL.min[!sigma] > eps)) {
        return(list(beta=beta,converged=FALSE))
      } else {
        return(list(beta=beta,converged=TRUE))
      }
    } else {
      l <- which.max(abs.dL.min)
      nabla.f <- S %*% beta + s
      sigma[l] <- TRUE
      theta[l] <- -sign(nabla.f[l])
    }
  }
}

# Inferrence for one penalty (loop over variables inside)
infer.edges <- function(X, constraint.mat, penalty, options) {
	# Initialization
	n <- nrow(X)
	p <- ncol(X) 
  
	if (is.null(options$initial.guess)) {
	  Beta <- matrix(0,(p-1), p)
	} else Beta <- options$initial.guess
  
	S.t <- var(X,na.rm=TRUE) 
	
	corrected.constraint <- constraint.mat
	corrected.constraint[constraint.mat == Inf] <- 0
	
	# Variable loop: build one LM for each variable
	for (k in 1:p) {
	  C.11 <- S.t[-k,-k]
	  C.12 <- S.t[-k,k]
    if (!is.null(constraint.mat)) {
      constraint.vect <- rep(0, p-1)
      constraint.vect[constraint.mat[k,-k] != 0] <- 1
      C.11 <- C.11 + diag(constraint.vect)*options$mu^2
      
      C.12 <- C.12 + corrected.constraint[k,-k]
	  }
	  
	  ## Function that performs the inference for one variable
	  results <- cLasso(C.11, C.12, penalty, beta=Beta[,k], max.size=min(n,p),
	                    max.it=options$max.it)
        
	  if (results$converged) {
	    Beta[,k] <- results$beta
	  } else  {
	    cat("out of convergence... stopping here")
	    return(list(Theta=NA, Beta=NA, n.edges=NA, loglik=NA, loglik.pen=NA,
                  BIC=NA, AIC=NA))
	  }
	}

  # Gathering all outputs together
  Theta       <- list()
  Theta.and   <- list()
  Theta.or    <- list()
  loglik      <- 0
  loglik.pen  <- 0
  total.df    <- 0
  AIC         <- 0
  
  # Build parameters list
  Theta <- matrix(0,p,p)
  dimnames(Theta) <- list(colnames(X), colnames(X))
  for (k in 1:p) {
    Theta[k,k]  <- 1/(S.t[k,k])
    Theta[-k,k] <- Beta[,k] * Theta[k,k] 
  }
  D  <- diag(Theta)
  Theta.tilde <- Theta %*% diag(D^(-1/2))
  loglik <- 0.5*(log(prod(D)) - 
                   sum(t(Theta.tilde) %*% var(X)%*%Theta.tilde))
	if (is.null(constraint.mat)) {
	  loglik.pen <- loglik - sum(abs(penalty * Theta))
	} else {
    cons <- matrix(constraint.mat[lower.tri(constraint.mat)|
                                        upper.tri(constraint.mat)],
                   nrow=nrow(constraint.mat)-1)
    loglik.pen <- loglik - sum(abs(penalty * Theta)) - options$mu/n*
      sum(sweep((Beta-cons)^2,2,diag(S.t),"/"))
	}  
  df <- (sum(abs(Theta)>0)- sum(abs(D)>0))/2
  AIC <- -2*loglik + 2*df
  loglik <- loglik + loglik
  loglik.pen <- loglik.pen + loglik.pen
  total.df <- total.df + df
    
  # Post-symetrization with the AND/OR rules
  Theta.and <- sign(Theta) * pmin(abs(Theta),t(abs(Theta)))
  Theta.or <- pmax(Theta,t(Theta)) - pmax(-Theta,-t(Theta))
  diag(Theta.or) <- diag(Theta.or)/2
  BIC <- -2*loglik + total.df * log(n)
  
  if (options$symmetrization=="AND") {
    Theta <- Theta.and
  } else Theta <- Theta.or
  
  # Get the number of edges inferred
  n.edges <- nb.edges(Theta)
  partial.cor <- - sweep(sweep(Theta, 1, sqrt(diag(Theta)), "/"), 2,
                         sqrt(diag(Theta)), "/")
  diag(partial.cor) <- 1
  
	return(list(Theta=Theta, partial.cor=partial.cor, Beta=Beta, n.edges=n.edges,
              loglik=loglik, loglik.pen=loglik.pen, BIC=BIC, AIC=AIC))
}

## Tests
# X <- matrix(rnorm(100),ncol=5)
# options <- set.options()
# X <- scale.default(X, TRUE, options$scale)
# penalty <- 0.1
# options$mu <- 0.5
# options$max.it <- 10
# out <- infer.edges(X, NULL, penalty, options)