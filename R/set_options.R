##### Handling options

set.options <- function(first.list=NULL, second.list=NULL, scale=TRUE,
                        verbose=TRUE, mu=NULL, penalties=NULL, penalty.min=1e-2,
                        penalty.max=NULL, n.penalties=100, edges.max=Inf, 
                        symmetrization=c("AND","OR"), initial.guess=NULL, 
                        max.it=NULL) {
  
  if (!is.null(penalties)) penalties <- sort(penalties, decreasing=TRUE)
  penalty.min <- max(penalty.min,.Machine$double.eps)
  
	all.options <- list(first.list=first.list, second.list=second.list, 
                      scale=scale, verbose=verbose, mu=mu, penalties=penalties,
	                    penalty.min=penalty.min, penalty.max=penalty.max,
	                    n.penalties=n.penalties, edges.max=edges.max, 
                      symmetrization=match.arg(symmetrization),
                      initial.guess=initial.guess, max.it=max.it)
  
  class(all.options) <- "internet3DOptions"
  return(all.options)
}

## Test
# set.options()


##### S3 method for internet3DOptions

print.internet3DOptions <- function(x,...) {
  cat("\n***** Parameters for internet3D\n\n")
  cat("    scaling                        : ", x$scale, "\n")
  if (!is.null(x$mu))
    cat("    mu                             : ", x$mu, "\n")
    cat("    number of penalties            : ", x$n.penalties, "\n")
    if (!is.null(x$penalty.min)&&!is.null(x$penalty.max))
      cat("    with min/max                   : ", x$penalty.min, "-",
          x$penalty.max, "\n")
  if (!(x$edges.max==Inf))
    cat("    maximum number of edges        : ", x$edges.max, "\n")
  if (x$symmetrization=="AND")
    cat("    symmetrization rule            :  AND\n")
  else cat("    symmetrization rule            :  OR\n")
}

## Test
# print.internet3DOptions(set.options())