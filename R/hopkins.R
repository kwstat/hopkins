# hopkins.R

## ---------------------------------------------------------------------------

#' Hopkins statistics for clustering tendency
#' 
#' @name hopkins
#' @aliases hopkins package-hopkins
#' @author Kevin Wright
#' @docType package
NULL

## ---------------------------------------------------------------------------

#' @title Hopkins statistic for clustering tendency
#'
#' @description Calculate Hopkins statistic for given data.
#'
#' Calculated values 0-0.3 indicate regularly-spaced data.
#' Values around 0.5 indicate random data.
#' Values 0.7-1 indicate clustered data.
#' 
#' CAUTION: This function does NOT center and scale the columns of X.
#' You may need to do this manually before using this function.
#' 
#' You should NOT set The parameter 'd'. It is included here to allow for
#' comparisons of \code{hopkins::hopkins()} and \code{clustertend::hopkins()}.
#'
#' The data \code{U} is also not normally set by the user.  It is included
#' here to allow for unit testing and also for customization of the
#' uniformly-sampled points (e.g. enlarged by 5 percent as suggested by
#' some authors).
#' 
#' Some authors suggest sampling less than 10 percent of points.
#' Others suggest m>10 points to avoid small-sample problems.
#' The distribution of Hopkins statistic requires that nearest neighbors
#' to the selected points be mutually independent, so that only a few of
#' the points can be marked.  The distribution of Hopkins statistic is
#' Beta(m,m), independent of the dimensionality of the data d.
#'
#' Cross & Jain say "The m sampling points are few enough in number,
#' relative to n (the number of events), that their presence does not
#' materially affect  the overall density.  Ratios of at least 10 to 1
#' and preferably 20 to 1 are used in the literature.
#' On the other hand, it seems that m should be at least 10 in order
#' to avoid any small sample problems with the distributions of the
#' statistics.  This effectively  limits the methods to problems with
#' at least 100 events.  In high dimensions, very little can be said
#' about data sets that are sparser than that."
#' 
#' Note:
#' 
#' Comparison of \code{hopkins::hopkins()} and \code{clustertend::hopkins()}.
#' 
#' The `hopkins::hopkins()` function uses distances^d (where "distance"
#' is the Euclidean distance between points and "d" is the number of
#' columns in the data).
#' The value returned is: Hopkins statistic.
#' 
#' The `clustertend::hopkins()` function uses distances^1.
#' The value returned is: 1 - Hopkins statistic.
#' 
#' @param X Data (matrix or data.frame) to check clusterability.
#' 
#' @param m Number of rows to sample from X. Default is 1/10th the number of rows of X.
#' 
#' @param d Dimension of the data (number of columns of X).
#' 
#' @param U Data containing \code{m} uniformly-sampled points.
#' 
#' @return The value of Hopkins statistic.
#' 
#' @author Kevin Wright
#' 
#' @examples
#' set.seed(1)
#' hopkins(iris[, -5], m=15) # .9952293
#' 
#' @references
#' Hopkins, B. and Skellam, J.G., 1954.
#' A new method for determining the type of distribution of plant individuals.
#' Annals of Botany, 18(2), pp.213-227.
#' 
#' Cross, G. R., and A. K. Jain. (1982).
#' Measurement of clustering tendency.
#' Theory and Application of Digital Control. Pergamon, 1982. 315-320.
#'
#' @importFrom pdist pdist
#' @importFrom stats runif
#' @export 
hopkins <- function (X, m=ceiling(nrow(X)/10), d=ncol(X), U=NULL) {
  
  if (!(is.matrix(X)) & !(is.data.frame(X))) 
    stop("X must be data.frame or matrix")

  if (m >= nrow(X)) 
    stop("m must be no larger than num of samples")

  if(missing(U)) {
    # U is a matrix of column-wise uniform values sampled from the space of X
    colmin <- apply(X, 2, min)
    colmax <- apply(X, 2, max)    
    U <- matrix(0, ncol = ncol(X), nrow = m)
    for (i in 1:ncol(X)) {
      U[, i] <- runif(m, min = colmin[i], max = colmax[i])
    }
  } else {
    # The user has provided the uniform values.
  }

  # Random sample of m rows in X (without replacement)
  k <- sample(1:nrow(X), m)
  W <- X[k, , drop=FALSE]   # Need 'drop' in case X is single-column
  
  # distance between each row of W and each row of X
  dwx <- as.matrix(pdist(W,X))
  # Caution: W[i,] is the same point as X[k[i],] and the distance between them is 0,
  # but we do not want to consider that when calculating the minimum distance
  # between W[i,] and X, so change the distance from 0 to Inf
  for(i in 1:m) dwx[i,k[i]] <- Inf
  # distance from each row of W to the NEAREST row of X
  dwx <- apply(dwx, 1, min)
  
  # distance between each row of U and each row of X
  dux <- as.matrix(pdist(U,X)) # rows of dux refer to U, cols refer to X
  # distance from each row of U to the NEAREST row of X
  dux <- apply(dux, 1, min)

  # You would think this would be faster, but it is not for our test cases:
  # stat = 1 / (1 + sum(dwx^d) / sum( dux^d ) )
  
  return( sum(dux^d) / sum( dux^d + dwx^d ) )
}

# ----------------------------------------------------------------------------

#' @title Calculate the p-value for Hopkins statistic
#'
#' Calculate the p-value for Hopkins statistic
#'
#' Under null hypothesis of spatial randomness, Hopkins statistic has a
#' Beta(m,m) distribution, where 'm' is the number of events/points sampled.
#' This function calculates the p-value for the statistic.
#'
#' @param x Observed value of Hopkins statistic
#' @param n Number of events/points sampled.
#' @return A p-value between 0 and 1.
#' @author Kevin Wright
#' @examples
#' hopkins.pval(0.21, 10) # .00466205
#' @references 
#' Michael T. Gastner (2005).
#' Spatial distributions: Density-equalizing map projections, facility location, and two-dimensional networks.
#' Ph.D. dissertation, Univ. Michigan (Ann Arbor, 2005).
#' http://hdl.handle.net/2027.42/125368
#'
#' @importFrom stats pbeta
#' @export 
hopkins.pval <- function(x,n) {
  if(x > 0.5)
    1 - (pbeta(x, n, n) - pbeta(1-x, n, n) )
  else
    1 - (pbeta(1-x, n, n) - pbeta(x, n, n) )
}

