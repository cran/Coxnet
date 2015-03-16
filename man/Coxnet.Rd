\name{Coxnet}
\alias{Coxnet}
\title{
Fit a Cox Model with Various Regularization Forms
}
\description{
Fit a Cox model regularized with net, elastic-net or lasso penalty.
In addition, it efficiently solves an approximate L0 variable selection based on truncated likelihood function (L0-Trunc).
Moreover, it can also handle the adaptive version of these regularization forms, such as adaptive lasso and net adjusting for signs of linked coefficients.
The regularization path is automatically computed.
}
\usage{
Coxnet(x, y, Omega = NULL, penalty = c("Lasso","Enet", "Net"),
  alpha = 1, lambda = NULL, nlambda = 50, rlambda = NULL,
  nfolds = 1, foldid = NULL, itrunc = TRUE, adaptive = c(FALSE,TRUE), aini = NULL,
  alambda = NULL, nalambda = 10, isd = TRUE, ifast = TRUE, keep.beta = FALSE,
  thresh = 1e-7, maxit = 1e+5)
}
\arguments{
  \item{x}{input matrix, each row is an observation vector.
}
  \item{y}{response variable. \code{y} should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored.
}
  \item{Omega}{adjancy matrix with zero diagonal, used for \code{penalty = "Net"} to calculate Laplacian matrix.
}
  \item{penalty}{penalty type. Can choose \code{"Net"}, \code{"Enet"} and \code{"Lasso"}. For \code{"Net"}, need to specify \code{Omega}; otherwises, \code{"Enet"} is performed.
}
  \item{alpha}{ratio between L1 and Laplacian for \code{"Net"}, or between L1 and L2 for \code{"Enet"}. Can be zero and one. For \code{penalty = "Net"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*(\beta^{T}L\beta)},}
where \eqn{L} is a Laplacian matrix calculated from \code{Omega}. For \code{adaptive = c(FALSE,TRUE)}, its calculation is also based on an initial estimate of \eqn{\beta} from regularized Cox model with \code{penalty = "Enet",alpha = 1} to adjust for signs of coefficents. For \code{penalty = "Enet"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*||\beta||_2}}.
}
  \item{lambda}{a user supplied decreasing sequence. If \code{lambda = NULL}, a sequency of \code{lambda} is generated based on \code{nlambda} and \code{rlambda}.
}  
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{lambda.max} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{nfolds}{number of folds. Default is 10. Smallest value allowable is \code{nfolds = 3}.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
  \item{itrunc}{logical flag for approximate L0 based on truncated likelihood function. Default is \code{itrunc = TRUE}.
}
  \item{adaptive}{logical flags for adaptive version. Default is \code{adaptive = c(FALSE, TRUE)}. The first element is for adaptive on \eqn{\beta} in L1 and the second for adjusting for signs of linked coefficients in Laplacian matrix.
}
  \item{aini}{a user supplied initial estimate of \eqn{\beta}. It is a list including \code{wbeta} for adaptive L1 and \code{sgn} for adjusting Laplacian matrix. If \code{aini = NULL}, \code{aini} is generated from regularized Cox model with \code{penalty = "Enet",alpha = 1}.
}
  \item{alambda}{a user supplied decreasing sequence used in initial estimate. If \code{alambda = NULL}, the sequency of \code{alambda} is generated based on \code{nalambda} and \code{rlambda}.
}
  \item{nalambda}{number of \code{alambda} values. Default is 10.
}

  \item{isd}{logical flag for \code{x} variable standardization, prior to fitting the model sequence. The coefficients \eqn{\beta} are always returned on the original scale. Default is \code{isd = TRUE}.
}
  \item{ifast}{logical flag for efficient calculation of risk set updates. Default is \code{ifast = TRUE}.
}
  \item{keep.beta}{logical flag for returning all the estimates of \eqn{\beta} for each \code{lambda} value. For \code{keep.beta = FALSE}, only return the estimate with the largest cross-validation partial likelihood and another based on one standard error rule.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-7}.
}
  \item{maxit}{Maximum number of iterations for coordinate descent. Default is \code{10^5}.
}
}
\details{
  One-step coordinate descent algorithm is applied for each \code{lambda}. \code{ifast = TRUE} adopts an efficient way to update risk set and sometimes the algorithm ends before all \code{nlambda} values of \code{lambda} have been evaluated. To evaluate small values of \code{lambda}, \code{ifast = FALSE} can be used. The two methods only affect the efficiency of algorithm, not the estimates.
  
  Cross-validation partial likelihood is used for tuning parameters. For \code{itrunc =  TRUE}, an approximate L0 variable selection based on truncated likelihood function (L0-Trunc) is also performed.
}
\value{
An object with S3 class \code{"Coxnet"}.
  \item{Beta}{a sparse Matrix of coefficients, stored in class "dgCMatrix".}
  \item{Beta0}{coefficients based on L0-Trunc, for \code{itrunc = TRUE}.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. For cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse} and one \code{index} with \code{**} and \code{*} indicating the largest \code{cvm} and one standard error respectively. \code{index = ***} implies that the largest \code{cvm} is same as the one based on one standard error.}
  \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and \code{nzero} based on L0-Trunc, for \code{itrunc = TRUE}.}
  
  \item{lambda.max}{value of \code{lambda} that gives maximum \code{cvm}.}
  \item{lambda.1se}{value of \code{lambda} based on one standard error.}
  \item{lambda.opt}{value of \code{lambda} based on L0-Trunc.}
  
  \item{cv.nzero}{\code{cvm} with length of number of non-zero components of \code{Beta0}. The kth value of \code{cv.nzero} corresponds to retaining k largest non-zero coefficients (absolute values) in \code{Beta0}. The optimal number of non-zero is selected by the maximum value of \code{cv.nzero} at \code{lambda = lambda.opt}.}
  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging).}
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
  \url{http://www3.stat.sinica.edu.tw/statistica/j24n3/j24n319/j24n319.html}\cr
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.Coxnet}}, \code{\link{coxsplit}}
}
\examples{
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=Coxnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)
}
\keyword{Cox model}
\keyword{Regularization}
