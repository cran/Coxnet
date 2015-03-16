\name{Coxnet-package}
\alias{Coxnet-package}
\docType{package}
\title{
Regularized Cox Model
}
\description{
This package fits Cox model regularized with net (L1 and Laplacian), elastic-net (L1 and L2) or lasso (L1) penalty. In addition, it efficiently solves an approximate L0 variable selection based on truncated likelihood function. Moreover, it can also handle the adaptive version of these regularization forms, such as adaptive lasso and net adjusting for signs of linked coefficients.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients.
}
\details{
\tabular{ll}{
Package: \tab Coxnet\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-03-16\cr
License: \tab GPL (>= 2)\cr
}
Two functions:
\code{\link{Coxnet}}, \code{\link{print.Coxnet}}, \code{\link{coxsplit}}\cr
}
\author{
Xiang Li, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@columbia.edu>
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent, Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}\cr
  Sun, H., Lin, W., Feng, R., and Li, H. (2014)
  \emph{Network-regularized high-dimensional cox regression for analysis of genomic data, Statistica Sinica.}\cr
  \url{http://www3.stat.sinica.edu.tw/statistica/j24n3/j24n319/j24n319.html}
}

\keyword{Cox model}
\keyword{Regularization}
\keyword{package}
\examples{
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)
fiti=Coxnet(x,y,penalty="Lasso")  # Lasso
}

