\name{clustertend-package}
\alias{clustertend-package}
\alias{clustertend}
\docType{package}
\title{
Check the Clustering Tendency
}
\description{
Calculate some statistics aiming to help analyzing the clustering tendency of the given data. In the first version, Hopkins' statistic is implemented.
}
\details{
\tabular{ll}{
Package: \tab clustertend\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2015-05-17\cr
License: \tab GPL (>= 2)\cr
}
Preprocess your data into a dataframe or matrix form. Then several statistics about clustering tendency can be calculated. In the fisrt version, we only provided calculating function of Hopkins' statistic.
}
\author{
Luo YiLan, Zeng RuTong

Maintainer: Zeng RuTong <670273197@qq.com>
}
\references{
Lawson, R.G. and Jurs, P.C.(1990) New index for clustering tendency and its application to chemical problems. Journal of Chemical Information and Computer Sciences. (Journal of Chemical Information and Computer Sciences, 1990, 30(1):36-41)
}
\keyword{ package }
\examples{
x<-matrix(runif(200,1,100),50,4);
hopkins(x,n=10)
}
