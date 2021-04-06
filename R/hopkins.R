
#' @title Calculate the Hopkins' statistic
#'
#' Calculate the Hopkins' statistic of given data. 'n' can be set to see whether this statistic converges.
#'
#' Sample data must be preprocessed into dataframe or matrix form before given as the value of parameter "data".
#'
#' @param data a data frame or a matrix of the sample
#' @param n an integer, the number of points selected from sample space which is also the number of points selected from the given sample(data) 
#' @param byrow logical. If FALSE(the default)the variables is taken by columns, otherwise the variables is taken by rows.
#' @param header logical. If FALSE(the default) the first column(or row) will be deleted in the calculation 
#' @return the number of Hopkins' statistic will be shown in the CW.
#' @author Luo YiLan, Zeng RuTong 670273197
#' 
#' @examples 
#' x<-matrix(runif(200,1,100),50,4);
#' hopkins(x,n=10)
#' 
#' @references 
#' Lawson, R.G. and Jurs, P.C.(1990) New index for clustering tendency and its application to chemical problems. Journal of Chemical Information and Computer Sciences. (Journal of Chemical Information and Computer Sciences, 1990, 30(1):36-41)
#' @importFrom stats dist runif
#' @export
#' 
hopkins <- function(data,n,byrow=F,header=F) 
{
  if(is.data.frame(data))
    data<-as.matrix(data)
  if (!(is.matrix(data)))
    stop("data must be data.frame or matrix") 
  if(n>=nrow(data))
    stop("n must be no larger than num of samples")
  if(byrow==T) data<-t(data)
  if(header==T) data<-data[-1,]
  c<-apply(data,2,min)#minimum value per colume
  d<-apply(data,2,max)
  p<-matrix(0,ncol=ncol(data),nrow=n)#n vectors of space
  for(i in 1:ncol(data))
  {
    p[,i]<-runif(n,min=c[i],max=d[i])
  }
  
  #k<-round(runif(n,1,nrow(data)))
  k <- sample(1:nrow(data), n)
  
  q<-as.matrix(data[k,])
  distp=rep(0,nrow(data))
  #distq=rep(0,nrow(data)-1)
  distq=0;
  minp=rep(0,n)
  minq=rep(0,n)
  for(i in 1:n)
  {
    distp[1]<-dist(rbind(p[i,],data[1,]))
    minqi<-dist(rbind(q[i,],data[1,]))
    for(j in 2:nrow(data))
    {
      distp[j]<-dist(rbind(p[i,],data[j,]))
      error<-q[i,]-data[j,]
      if(sum(abs(error))!=0)
      {
        #distq[j]<-dist(rbind(q[i,],data[j,]))
        distq<-dist(rbind(q[i,],data[j,]))
        if(distq<minqi)
          minqi<-distq;
      }
    }
    minp[i]<-min(distp)
   # minq[i]<-apply(distq,1,min)
   minq[i]<-minqi;
  }
  list(H=(sum(minq)/(sum(minp)+sum(minq))))
  
}
