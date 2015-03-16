

#################################################
#####  Initial values for adaptive methods  #####
#################################################

coxini=function(x, y, alambda=NULL, nalambda=10, rlambda=NULL, isd=TRUE, ifast=TRUE){
  
  N0=nrow(x);p=ncol(x)
  ifast=as.integer(ifast)
  
  xscale=rep(1, p)
  if (isd) {
    tem=scaleC(x)
    xscale=tem$sd;x=tem$x;
    rm(tem)
  }
  
  if (is.null(alambda)) {
    prep0=coxprep(x, y);wbeta=rep(1, p)
    ### Lambda path
    lambda_max=max_lambdaC(prep0$x, prep0$tevent, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, 1, wbeta, N0)
    lambda_min=ifelse(is.null(rlambda), ifelse(N0>=p, lambda_max*0.0001, lambda_max*0.01), lambda_max*rlambda)
    alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  }
  
  repeat {
    outi=coxEnet(x, y, alpha=0.0, lambda=alambda, keep.beta=T)
    if(!is.null(outi))break
    alambda=alambda*2.0
  }
  
  indexi=ncol(outi$Beta)
  beta0=outi$Beta[, indexi]
  wbeta=1/abs(beta0);sgn=sign(beta0[1:p])
  return(list(wbeta=wbeta, sgn=sgn, lambda=alambda[indexi]))
}



###############
###  local  ###
###############

# locoxini=function(x, y, w, w0, h, alambda, nalambda, rlambda=NULL, alinear=FALSE, isd=TRUE, thresh3=1e-15){
#   # h=0.1;rlambda=NULL;nfolds=1;p=ncol(x);alambda=NULL;thresh3=1e-15;isd=T
#   nw0=length(w0);p=ncol(x)
#   
#   ### Lambda path
#   if(is.null(alambda)){
#     alambda=locoxini.lambda(x, y, w, w0, h, nalambda, rlambda, alinear, isd)
#   }else{nalambda=length(alambda)}
#   
#   outi=list();thresh3i=thresh3;ithresh3=0
#   repeat{
#     nalambdai=nalambda;alambdai=alambda;ithresh3=ithresh3+1
#     for(iw0 in 1:nw0){
#       if(alinear){u=cbind(x, (w-w0[iw0])*x)}else{u=x}
#       outi[[iw0]]=locoxt.enet(u, y, w, w0[iw0], h, alpha=0.0, lambda=alambdai, isd=isd, thresh2=thresh3)
#       nalambdai=length(outi[[iw0]]$lambda);alambdai=alambda[1:nalambdai]
#       if(nalambdai==0)break
#     }
#     thresh3i=thresh3i/10
#     if(all(sapply(outi, function(x){length(x$lambda)>0})))break
#     if(ithresh3==10)stop("Need larger lambda!")
#   }
#   
#   ### re-fit using thresh2=0
#   outi=list()
#   for(iw0 in 1:nw0){
#     if(alinear){u=cbind(x, (w-w0[iw0])*x)}else{u=x}
#     outi[[iw0]]=locoxt.enet(u, y, w, w0[iw0], h, alpha=0.0, lambda=alambdai, keep.beta=T, isd=isd, thresh2=0)
#   }
#   
#   beta0=sapply(outi, function(x){x$Beta[, nalambdai]})[1:p, ]
#   wbeta=matrix(1/abs(beta0), ncol=nw0);sgn=matrix(sign(beta0), ncol=nw0)
#   
#   return(list(wbeta=wbeta, sgn=sgn, lambda=alambdai[nalambdai]))
# }
# 
# 
# locoxini.lambda=function(x, y, w, w0, h, nalambda=10, rlambda=NULL, alinear=FALSE, isd=TRUE){
#   
#   ### Lambda path
#   nw0=length(w0);N0=nrow(x);p=ncol(x)
#   tem_lambda_max=numeric(nw0)
#   if(alinear){p=p*2}
#   
#   for(iw0 in 1:nw0){
#     
#     if(alinear){u=cbind(x, (w-w0[iw0])*x)}else{u=x}
#     ### scaleC and standardized
#     xscale=rep(1, p)
#     if(isd){
#       tem=scaleC(u)
#       xscale=tem$sd;u=tem$x
#     }
#     
#     ###  Full data  ###
#     prep0=locoxprept(u, y, w, w0[iw0], h)
#     tem_lambda_max[iw0]=max_loclambdaC(prep0$x, prep0$tevent, prep0$Kh, prep0$Kh1, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, 1.0, rep(1.0, p), N0)
#   }
#   
#   lambda_max=max(tem_lambda_max);temm=min(tem_lambda_max)
#   lambda_min=ifelse(is.null(rlambda), ifelse(N0>=p, temm*0.0001, temm*0.05), temm*rlambda)
#   lambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
#   return(lambda)
# }


