#' Set and print options for network inference.
#'
#' @param first.list list of edges supposed to be present (array with two 
#' columns and a number of rows corresponding to the number of edges with a 
#' prior)
#' @param second.list list of edges supposed to be absent (array with two 
#' columns and a number of rows corresponding to the number of edges with a 
#' prior)
#' @param scale logical. Should the data be scaled to unit variance prior the 
#' analysis
#' @param verbose logical. Should messages be printed during the learning
#' process
#' @param mu positive number. Regularization parameter for the L2 penalty on
#' edges with a prior assumption
#' @param penalties a sequence of decreasing positive numbers to control the L1
#' regularization
#' @param penalty.min if \code{penalties=NULL}, minimum value for the L1
#' regularization parameter
#' @param penalty.max if \code{penalties=NULL}, maximum value for the L1
#' regularization parameter
#' @param n.penalties if \code{penalties=NULL}, number of values to include in
#' the sequence of positive numbers that control the L1 regularization
#' @param edges.max maximum number of edges to stop the learning process
#' @param symmetrization symmetrization rule, to be chosen between \code{"AND"}
#' and \code{"OR"}
#' @param initial.guess initial guess (if existing) for the solution of the
#' optimization problem
#' @param max.it maximum number of iterations allowed to solve the optimization 
#' problem
#' @return object of type \code{internet3DOptions} that can be used within the
#' functions \code{\link{build.network}} or \code{\link{bootstrap.build}}
#' @seealso \code{\link{build.network}}, \code{\link{bootstrap.build}}
#' @export
#' @examples
#' set.options()
#' set.options(mu = 1)

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

##### S3 method for internet3DOptions

#' @rdname set.options
#' @param x a \code{internet3DOptions} object
#' @param ... not used
#' @export

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
