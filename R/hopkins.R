# hopkins.R

#' Hopkins statistics for clustering tendency
#' 
#' @name hopkins
#' @aliases hopkins package-hopkins
#' @author Kevin Wright
NULL

## ---------------------------------------------------------------------------

#' Hopkins statistic for clustering tendency
#'
#' Calculate Hopkins statistic for given data.
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
#' @param k kth nearest neighbor to find.
#' 
#' @param U Data containing \code{m} uniformly-sampled points.
#'
#' @param method Either "simple" or "torus".
#' 
#' @return The value of Hopkins statistic.
#' 
#' @author Kevin Wright
#' 
#' @examples
#' set.seed(1)
#' hopkins(iris[, 1:4], m=15) # .9952293
#'
#' hop <- rep(NA, 100)
#' for(i in 1:100){
#'   hop[i] <- hopkins(iris[,1:4], m=8)
#' }
#' mean(hop)
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
#' @importFrom donut nnt
#' @importFrom RANN nn2
#' @importFrom pdist pdist
#' @importFrom stats runif
#' @export 
hopkins <- function (X, m=ceiling(nrow(X)/10), d=ncol(X), k=1, U=NULL, method="simple") {
  
  if (!(is.matrix(X)) & !(is.data.frame(X))) 
    stop("X must be data.frame or matrix")

  if (m >= nrow(X)) 
    stop("m must be MUCH SMALLER than the number of samples")

  if(missing(U)) {
    # U is a matrix of column-wise uniform values sampled from the space of X
    colmin <- apply(X, 2, min)
    colmax <- apply(X, 2, max)    
    U <- matrix(0, ncol = ncol(X), nrow = m)
    for (i in 1:ncol(X)) {
      U[, i] <- runif(m, min = colmin[i], max = colmax[i])
    }
  } else {
    # The user has provided the uniform values in U.
  }

  # Random sample of m rows in X (without replacement)
  j <- sample(1:nrow(X), m)
  W <- X[j, , drop=FALSE]   # Need 'drop' in case X is single-column

  if(method=="simple") {
    # distance between each row of W and each row of X
    dwx <- as.matrix(pdist(W,X))
    # Caution: W[i,] is the same point as X[j[i],] and the distance between them is 0,
    # but we do not want to consider that when calculating the minimum distance
    # between W[i,] and X, so change the distance from 0 to Inf
    for(i in 1:m) dwx[i,j[i]] <- Inf
    # distance from each row of W to the NEAREST row of X
    dwx <- apply(dwx, 1, min)
    
    # distance between each row of U and each row of X
    dux <- as.matrix(pdist(U,X)) # rows of dux refer to U, cols refer to X
    # distance from each row of U to the NEAREST row of X
    dux <- apply(dux, 1, min)
    
    hop <- sum(dux^d) / sum( dux^d + dwx^d )
  }
  if(method=="torus") {    
    rng <- t(apply(X,2,range))

    # Note: Since W is a sample from X, the 1st nearest point in X will
    # always be the same point with distance 0, so add 1 to k.
    nearw <- donut::nnt(X, W, k=k+1, torus=1:ncol(W), ranges=rng )
    dwx <- nearw$nn.dists[,k+1]

    # For U, find the 1st nearest point in X, k=1.
    nearu <- donut::nnt(X, U, k=k, torus=1:ncol(W), ranges=rng )
    dux <- nearu$nn.dists[,k]
    
    hop <- sum(dux^d) / sum( dux^d + dwx^d )
  }
  
  if(method=="boundedsphere" | method=="boundedcube"){
    # distance between each row of W and each row of X
    dwx <- as.matrix(pdist(W,X))
    for(i in 1:m) dwx[i,j[i]] <- Inf
    # distance from each row of W to the NEAREST row of X
    dwx <- apply(dwx, 1, min)
    
    # distance between each row of U and each row of X
    dux <- as.matrix(pdist(U,X)) # rows of dux refer to U, cols refer to X
    dux <- apply(dux, 1, min)

    if(method=="boundedsphere") const <- 1
    if(method=="boundedcube") const <- 2
    ukd <- (const*dux)^d
    wkd <- (const*dwx)^d
    N <- nrow(X)

    # I *think* this is the formula of Rotondi, page 560-561.
    # However, what happens if ukd=1? Then division by zero.
    # Also, this is supposed to have Beta(kM,kM) distribution, so it should
    # be between [0,1], but my example has values outside [0,1].
    # Should X be scaled to unit hypercube/hypersphere???
    hop <- (N-k+1) * sum( ukd/(1-ukd) ) /
      sum( (N-k+1)*ukd/(1-ukd) + (N-k)*(wkd/(1-wkd)) )
  }
  
  # You would think this would be faster, but it is not for our test cases:
  # stat = 1 / (1 + sum(dwx^d) / sum( dux^d ) )
  
  return( hop )
}

# ----------------------------------------------------------------------------

#' Calculate the p-value for Hopkins statistic
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

if(0){
  D <- 8 # dimension of data, columns(X)
  N <- 5000 # number of events, rows(X)
  M <- 8 # number of events sampled
  B <- 1000 # number of simulations

  #scale01 <- function(x){ (x-min(x))/(max(x)-min(x)) }
  set.seed(12345)
  hop1 <- hop2 <- NULL
  for(ii in 1:B){
    X <- matrix(runif(N*D, min=-1, max=1), ncol=D)
    # calculate radial distance from origin for each point
    rad <- apply(X, 1, function(x) sqrt(sum((x-0)^2)))
    X <- X[rad < 1.0,]

    # We need to scale the data into unit hypercube
    # X <- apply(X, 2, scale01)
    # Since this is a simulation study, we can use the first M rows
    # as random origins in the unit hypersphere
    hop1[ii] <- hopkins::hopkins(X[-c(1:M),], U=X[1:M,], m=M, d=D)
    hop2[ii] <- hopkins::hopkins(X[-c(1:M),], U=X[1:M,], m=M, d=D, method="boundedsphere")
  }

  # Now the plots
  plot(density(hop1), col="red", xlim=c(0,1), main="", xlab="", ylim=c(0,4))
  lines(density(hop2), col="blue")
  xv <- seq(0,1,length=100)
  lines(xv, dbeta( xv, M, M) , col="black", lwd=2)
  legend("topleft",
         c("Hopkins", "Modified Hopkins", "Beta(M,M)"),
         text.col=c("red","blue","black")
         )

}
