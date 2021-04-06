\name{hopkins}
\alias{hopkins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Hopkins' statistic
}
\description{
Calculate the Hopkins' statistic of given data. 'n' can be set to see whether this statistic converges.
}
\usage{
hopkins(data, n, byrow = F, header = F)
}
\arguments{
  \item{data}{
a data frame or a matrix of the sample
}
  \item{n}{
an integer, the number of points selected from sample space which is also the number of points selected from the given sample(data) 
}
  \item{byrow}{
logical. If FALSE(the default)the variables is taken by columns, otherwise the variables is taken by rows.
}
  \item{header}{
logical. If FALSE(the default) the first column(or row) will be deleted in the calculation 
}
}
\details{
Sample data must be preprocessed into dataframe or matrix form before given as the value of parameter "data".
}
\value{
the number of Hopkins' statistic will be shown in the CW.
}
\references{
Lawson, R.G. and Jurs, P.C.(1990) New index for clustering tendency and its application to chemical problems. Journal of Chemical Information and Computer Sciences. (Journal of Chemical Information and Computer Sciences, 1990, 30(1):36-41)
}
\author{
Luo YiLan, Zeng RuTong 670273197@qq.com
}
\examples{
x<-matrix(runif(200,1,100),50,4);
hopkins(x,n=10)
}
\keyword{ Statistics }
\keyword{ cluster }
