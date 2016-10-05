#' @title Geometric mean of a vector
#'
#' @authors Jacques van Helden and Mustafa AbuElQumsan
#'
#' @description \code{GeomMean} computes the geometric mean of a vector of numerical
#' values, using a robust implementation: (1) NA values are optionally removed;
#' (2) the computation is done on the log-transformed values to avoid limitations of floating
#'  point computation, such as exceeding MAXFLOAT (product of many values > 1),
#'  or rounding to 0 (product of many values < 1).
#'
#' @param x A vector of numerical values from which the geometric mean will be computed.
#' @param na.rm=FALSE if TRUE, NA values will be removed from the vector before computing
#' the geometric mean.
#' @return a numerical value with the geometric mean of the input vector.
#' If the input data contains NA values and the na.rm option is not activated,
#' the function returns NA.
#' @examples
#' # Use cases for the GeomMean() function: generate some representative vectors of
#' # values and compute their geometric mean
#'
#' # Test 1: a small vector of strictly positive real values
#' x <- rchisq(n=10, df=3)
#' print(x)
#' print(prod(x)^(1/length(x))) # Straight computation of the geometric mean
#' print(GeomMean(x)) # Computation with our function should give the same result
#'
#' # Test 2: include NA values in the vector
#' x <- c(rchisq(n=10, df=3), NA, NA)
#' print(x)
#' print(prod(x)^(1/length(x))) # Straight computation of the geometric mean
#' print(GeomMean(x)) # Computation with our function should give the same result
#'
#' # Test 3: include NA values in the vector, but remove them with the option na.rm
#' x <- c(rchisq(n=10, df=3), NA, NA)
#' print(x)
#' print(prod(x, na.rm = TRUE)^(1/length(na.omit(x)))) #' Straight computation of the geometric mean
#' print(GeomMean(x, na.rm = TRUE)) #' Computation with our function should give the same result
#'
#' # Test 4: show that GeomMean() circumvents a problem with product of many values >1
#' #' exceeding max floating point number
#' x <- c(rchisq(n=1000, df=12))
#' print(x)
#' print(prod(x)^(1/length(x))) #' Straight computation of the geometric mean
#' print(GeomMean(x)) #' Computation with our function should give the same result
#'
#' # Test 5: show that GeomMean() circumvents a problem with product of many
#' # values < 1 leading to 0 (another limitation of floating point).
#' x <- c(rchisq(n=1000, df=1))
#' print(x)
#' print(prod(x)^(1/length(x))) # Straight computation of the geometric mean
#' print(GeomMean(x)) # Computation with our function should give the same result
#' @export
GeomMean <- function(x, na.rm=FALSE) {

  #  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) ## Attention, this formula is from internet and it is false, because Zero values are discarded !!!

  # Filter NA values if specified
  if (na.rm) {
    x <- na.omit(x)
  }

  # If there is at least one zero value, the geometric mean is zero
  if (sum(x==0, na.rm=TRUE) > 0) {
    return(0)
  }

  # Negative values are not compatible with the geometric mean, because one would need to compute the root of a negative number
  if (sum(x < 0, na.rm=TRUE) > 0) {
    stop("Cannot compute geometric mean for a vector with negative values. ")
  }

  # Warning: the "naive" way to compute the geometric mean leads rapidly to Inf,
  # because the product of all values can exceeed the maximal floating point
  # number. For this reason, we do the computation on log-transformed values.
  # We first compute the arithmetic mean of log-transformed values, and then
  # take the exponential of this mean.
  return(exp(mean(log(x))))

  # Just for validation, the following code should give the same result when
  # the product of the values does not exceed MAXFLOAT.
  # (prod(x))^(1/length(x)) ## Not working if the product exceeds max float
}




