
libs(FNN, pdist, microbenchmark, mvtnorm)

X <- matrix(rnorm(1000*5), ncol=5)
XU <- matrix(rnorm(10*5), ncol=5)

microbenchmark(out1 <- knnx.index(X, XU, k=1))
microbenchmark(out2 <- get.knnx(X, XU, k=1))
microbenchmark(out3 <- pdist(XU, X))
head(out3@dist)

lib
    sigma<- function(v, r, p)
    {
      	V<- matrix(r^2, ncol=p, nrow=p)
    	  diag(V)<- 1
        V*v
    }

    X<- rmvnorm(1000, mean=rep(0, 20), sigma(1, .5, 20))
    print(system.time(knn.dist(X)) )
    print(system.time(knn.dist(X, algorithm = "kd_tree")))
