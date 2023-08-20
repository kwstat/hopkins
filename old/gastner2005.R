# Gastner2005.R

# The examples here come from:
# Michael T. Gastner (2005).
# Spatial distributions: Density-equalizing map projections, facility location, and two-dimensional networks
# Ph.D. dissertation, Univ. Michigan (Ann Arbor, 2005).
# http://hdl.handle.net/2027.42/125368

Three examples of spatial point patterns, reproduced from [13], are shown in Fig. A.1.

Panel (a) shows the locations of 65 Japanese black pine saplings in a square of side 5.7 m [179].

Panel (b) the locations of 62 redwood seedlings in a square of side 23 m [180, 181],

Panel (c) the pattern formed by the centers of mass of 42 cells from insect tissue [182, 181]

libs(spatstat.data)

# Systematic data, evenly spaced
data(cells)
cells = data.frame(x=cells$x, y=cells$y)
plot(cells, main="cells")

# Random data
data(japanesepines)
japanesepines <- data.frame(x=japanesepines$x, y=japanesepines$y)
plot(japanesepines, main="japanesepines")

# Clustered data
data(redwood)
redwood = data.frame(x=redwood$x, y=redwood$y)
plot(redwood, main="redwood")

hopkins(cells, 10) # .16-.28
hopkins(japanesepines, 10) # .38-.62
hopkins(redwood, 10) # .55-.90




# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

show_hopkins <- function (X, m, d=ncol(X), U=NULL) {
  
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

  # Random sample of rows in X
  k <- sample(1:nrow(X), m)
  W <- X[k, , drop=FALSE]   # Need 'drop' in case X is single-column

  # distance between each row of W and each row of X
  dwx <- as.matrix(pdist(W,X))
  # Caution: W[i,] is the same point as X[k[i],] and the distance between them is 0,
  # but we do not want to consider that when calculating the minimum distance
  # between W[i,] and X, so change the distance from 0 to Inf
  for(i in 1:m) dwx[i,k[i]] <- Inf

  plot(X, main="", col="gray20")
  # Show nearest events in X
  points(W, col="blue", pch="W", cex=.5)
  # From W, nearest points in X
  # Run the line setting some distances to Inf. Then:
  x_near_w <- apply(dwx, 1, which.min)
  x0 <- W[,1]; y0 <- W[,2]
  x1 <- X[x_near_w, 1] ; y1 <- X[x_near_w, 2]
  arrows(x0, y0, x1, y1, col="skyblue", length=0.1)
  w_distances <- sqrt((x1-x0)^2 + (y1-y0)^2)
  w_mids <- cbind((x0+x1)/2, (y0+y1)/2)
  text(w_mids+.02, labels=round(w_distances,3), cex=1, col="blue")
  
  
  # distance from each row of W to the NEAREST row of X
  dwx <- apply(dwx, 1, min)
  
  # distance between each row of U and each row of X
  dux <- as.matrix(pdist(U,X)) # rows of dux refer to U, cols refer to X

  # From U, nearest points in X
  points(U, col="red", pch="U", cex=.5)
  x_near_u <- apply(dux, 1, which.min)
  x0 <- U[,1]; y0 <- U[,2]
  x1 <- X[x_near_u, 1] ; y1 <- X[x_near_u, 2]
  arrows(x0, y0, x1, y1, col="pink", length=0.1)
  u_distances <- sqrt((x1-x0)^2 + (y1-y0)^2)
  u_mids <- cbind((x0+x1)/2, (y0+y1)/2)
  text(u_mids+.02, labels=round(u_distances,3), cex=1, col="red")

  # distance from each row of U to the NEAREST row of X
  dux <- apply(dux, 1, min)

  title( main=paste0(substitute(X),
                     ": Hopkins: ",
                     round( sum(dux^d) / sum( dux^d + dwx^d ) , 2) ))

  return(  sum(dux^d) / sum( dux^d + dwx^d ) )
}

i=0
par(mfrow=c(2,2))

i=i+1
#i=15
set.seed(i)
show_hopkins(cells, m=3)
show_hopkins(japanesepines, m=3)
show_hopkins(redwood, m=3)
frame()

mean_sd <- function(x) {round( c(mean(x), sd(x) ), 2)}


# 0.21 .05 systematic 
unlist(lapply(1:100, function(x) hopkins(cells,10))) %>% mean_sd
show_hopkins(cells, m=3)
#1 - (pbeta(.79, 10, 10) - pbeta(.21, 10, 10) ) # .00466
hopkins.pval(.79, 10)

# 0.50 .1 random
unlist(lapply(1:100, function(x) hopkins(japanesepines,10))) %>% mean_sd()
show_hopkins(japanesepines, m=3)
#1 - (pbeta(.51, 10, 10) - pbeta(.49, 10, 10) ) # .9296
hopkins.pval(.51, 10)

# .84 .05 clustered
unlist(lapply(1:100, function(x) hopkins(redwood,10))) %>% mean_sd       
show_hopkins(redwood, m=3)
#1 - (pbeta(.84, 10, 10) - pbeta(.16, 10, 10) ) # .000498
hopkins.pval(.84, 10)


