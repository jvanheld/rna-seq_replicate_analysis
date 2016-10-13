x<- c(100,200,50,60,NA)

test_that("Test GeomMean(x)",{

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


})

source("/home/el-qunsam/rna-seq_replicate_analysis/MultiVariateNGS/R/run_tests.R")
test_results <- test_dir("/home/el-qunsam/rna-seq_replicate_analysis/MultiVariateNGS/R", reporter="summary")
