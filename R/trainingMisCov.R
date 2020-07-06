#' trainingMisCov
#'
#' This is train algorithm, we consider three different cases
#' 1. CC; 2.WCC; 3.DR
#' The algorithm uses the hinge loss function to do minimization procedure
#' If CC case, return the estimated alpha based on complete data estimated parameters
#' If WCC case, return the estimated alpha based on complete data and
#' for propensity score which is also used in testing procedure;
#' If DR case, return (1) the contrained quadratic program solution alpha
#' (2) estimated parameters for propensity score (3) estimated parameters
#' of conditional distribution of x2 given X1 and Y. (4) Covariates of the
#' dataset after imputation- the large dataset
#'
#' @param misCovDat data set with missing covariates.
#'                  (X,V,R,Y1,Y), X is the fully observed covariates; V is the popential missing covariates; R is the missing indicator; Y1=Y the binary response
#' @param px dimsension of covariates X (totally observed)
#' @param kerType type of kerenl function, "RBF"(if choose) or "linear" (else choose, automatically)
#' @param kerMethod method to use "CC", "WCC", and "DR"
#' @param lambda tuning parameter
#' @param sigma parameter for RBF kernel, also act as tuning parameter
#' @param PSFunPath Path of function used to estimated propensity score
#' @param IMPFunPath Path of imputation function used to generte imputation data
#' @param B The imputation time for DR (B>=1)
#' @param testPurpose test purpose, for crossivalidation, the loss function is chosen as the phi loss;
#'                                  for test (default value) , the loss function is chosen as the classification loss.
#'
#' @return trainRes
#'
#' @import e1071
#' @import modopt.matlab
#' @importFrom kernlab rbfdot vanilladot kernelMatrix
#' @importFrom dplyr tibble
#' @import stats
#' @export
#'
trainingMisCov=function(misCovDat,px,kerType,kerMethod,lambda,sigma,PSFunPath,IMPFunPath,B,testPurpose){
  pTotal=dim(misCovDat)[2] # dimension of (X,V,R,Y1,Y)
  misCovDat=as.matrix(misCovDat)
  #colnames(misCovDat)=NULL
  X=misCovDat[,1:px] # the fully observed covariates
  V=misCovDat[,(px+1):(pTotal-3)] # the popential missing covariates
  Z=cbind(X,V) # the covariates
  p=ncol(Z) # dimension of covariates
  R=misCovDat[,(pTotal-2)]   #  missing indicator
  Y=misCovDat[,pTotal] # the response (-1 and 1)
  n=nrow(misCovDat) # sample size of the traing data

  ## Get kernel matirx
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel


  YCC=Y[R==1] # Y corresponding complete cases
  ZCC=Z[R==1,] # covariates corresponding complete cases
  XCC=Z[R==1,1:px] # Y corresponding complete cases
  VCC=Z[R==1,(px+1):p] # covariates corresponding complete cases
  YCCMatrix=diag(YCC)
  nCC=length(YCC)  # sample size of the traing data

  C=1/(2*nCC*lambda) # the parameter in the minimization process.

  # Get kernel matrix for WCC and CC
  if(kerMethod!="DR"){
    K=kernelMatrix(kerFun,ZCC)
    # K=matrix(0,nCC,nCC)
    # for(i in 1:nCC){
    #   for(j in 1:i){K[i,j]=kerFun(ZCC[i,],ZCC[j,]);K[j,i]=K[i,j]}} # end of kernel matrix computation
  }# end of if

  ## CC method
  if(kerMethod=="CC"){
    ## Constrainted optmiazed program for CC
    # L=(1/2)*alpha'*D*K*D*alpha-alpha'1, subgect to 0<alpha<C*1,
    # where D=diag(YCC), 1 is a vector, 0 is a vector.
    res=quadprog(YCCMatrix%*%(K+diag(0.0001,nCC))%*%YCCMatrix,rep(-1,nCC),A=NULL,b=NULL,Aeq=rep(0,nCC),beq=0,
                 lb=rep(0,nCC),ub=rep(C,nCC))
    # Get the estimator of alpha
    hatAlphaCC=res$x
    trainRes=list(hatAlphaCC=hatAlphaCC,ZCC=Z[R==1,],YCC=YCC,trainSampleSize=n,
                  sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                  px=px,testPurpose=testPurpose)
  } # end of CC


  ## Estimate the propensity score model for WCC and DR
  if(kerMethod!="CC"){
    PSFun= readRDS(PSFunPath)
    PSmodel=PSFun(R,X,Y)
    piHat=PSmodel(X,Y)}

  ## WCC method
  if(kerMethod=="WCC"){
    ## Constrainted optmiazed program for WCC
    # L=alpha'1-(1/2)*alpha'*D*K*D*alpha, subgect to 0<alpha<C*W*1,
    # where D=diag(YCC), W=diago(piHat), 1 is a vector, 0 is a vector.
    W=diag(1/piHat[R==1]) # weight matrix diag(1/pi_1,...1/pi_nCC)
    res=quadprog(YCCMatrix%*%(K+diag(0.0001,nCC))%*%YCCMatrix,rep(-1,nCC),A=NULL,b=NULL,Aeq=rep(0,nCC),beq=0,
                 lb=rep(0,nCC),ub=W%*%rep(C,nCC))
    # Get the estimator of alpha
    hatAlphaWCC=res$x
    trainRes=list(hatAlphaWCC=hatAlphaWCC,ZCC=ZCC,YCC=YCC,R=R,X=X,Y=Y,trainSampleSize=n,
                  sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                  PSFunPath=PSFunPath,px=px,testPurpose=testPurpose)
  } # end of WCC

  ## DR method
  if(kerMethod=="DR"){

    ## Imputation process
     # imputation time is B
    IMPFun=readRDS(IMPFunPath)
    modelIMP=IMPFun(XCC,YCC,VCC)

    pv=p-px
    df2=matrix(0,n*B,pv+1)
    for(i in 1:n){
      if(is.vector(X)){matIMP=modelIMP(B,X[i],Y[i])}else{matIMP=modelIMP(B,X[i,],Y[i])}
      df2[((i-1)*B+1):(i*B),]=as.matrix(matIMP)}
    df2=as.data.frame(df2)
    names(df2) <-c(paste0("v",1:pv), "y")
   
    repFun=function(x){return(rep(x,each=B))}
    if(is.vector(X)){ImpData=as.matrix(cbind(repFun(X),df2[,1:pv],repFun(R),df2$y,repFun(piHat)))}else{
    ImpData=as.matrix(cbind(apply(X,2,repFun),df2[,1:pv],repFun(R),df2$y,repFun(piHat)))}


    ## Combine the original missing data and imputating data.
    #  The new data set will be dimension of [(n+1)*B]*4
    N=n+n*B
    fulData=matrix(0,N,pTotal)
    fulData[1:n,]=cbind(Z,R,Y,piHat)
    fulData[(n+1):N,]=ImpData

    Zful=fulData[,1:(pTotal-3)] # covariate in the new data
    Rful=fulData[,(pTotal-2)] # missing indicator in the new data
    Yful=fulData[,(pTotal-1)] # response in the new data
    piHatful=fulData[,pTotal] # estimated propensity score in the new data



    ## Get vecMu and vecNu
    vecMu=rep(0,N)
    vecNu=rep(0,N)
    # mu=(B+1)*R_i/PiHat_i*I(Y_i==1) 1<=i<=n;
    #   =N/(n*B)*[R_i*(1-PiHat_i)/PiHat_i*I(Y_i==-1)+(1-M_i)*I(Y_i==1)], (n+1)<=i<=n*(B+1);
    vecMu[1:n]=(B+1)*R/piHat*as.numeric(Y==1)
    vecMu[(n+1):N]=N/(n*B)*(Rful[(n+1):N]*((1-piHatful[(n+1):N])/piHatful[(n+1):N])*as.numeric(Yful[(n+1):N]==-1)+(1-Rful[(n+1):N]**as.numeric(Yful[(n+1):N]==1)))

    # nu=(B+1)*M_i/PiHat_i*I(Y_i==0), 1<=i<=n;
    #   =N/(n*B)*[M_i*(1-PiHat_i)/PiHat_i*I(Y_i==1)+(1-M_i)*I(Y_i==-1)], (n+1)<=i<=n*(B+1).
    vecNu[1:n]=(B+1)*R/piHat*as.numeric(Y==-1)
    vecNu[(n+1):N]=N/(n*B)*(Rful[(n+1):N]*((1-piHatful[(n+1):N])/piHatful[(n+1):N])*as.numeric(Yful[(n+1):N]==1)+(1-Rful[(n+1):N]**as.numeric(Yful[(n+1):N]==-1)))


    # Get kernel matrix
    # K=matrix(0,N,N)
    # for(i in 1:N){
    #   for(j in 1:i){K[i,j]=kerFun(Zful[i,],Zful[j,]);K[j,i]=K[i,j]}}
    K=kernelMatrix(kerFun,Zful)

    # Q=(K,-K) Q is the matrix of quadratic optimization problem
    #   (-K,K)
    Qmatrix=matrix(0,2*N,2*N)
    Qmatrix[1:N,1:N]=K
    Qmatrix[1:N,(N+1):(2*N)]=-K
    Qmatrix[(N+1):(2*N),1:N]=-K
    Qmatrix[(N+1):(2*N),(N+1):(2*N)]=K

    ## Constrainted optmiazed program
    res=quadprog(Qmatrix+diag(0.0001,2*N), rep(-1,2*N), NULL, NULL,Aeq=rep(0,2*N),beq=0,
                 lb=rep(0,2*N),ub=1/(2*N*lambda)*c(vecMu,vecNu)+rep(0.001,2*N))

    # Get the estimator of alpha
    hatAlphaDR=res$x[1:N]-res$x[(N+1):(2*N)]

    trainRes=list(hatAlphaDR=hatAlphaDR,Zful=Zful,ZCC=ZCC,YCC=YCC,R=R,X=X,Y=Y,
                  trainSampleSize=n,sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                  impTime=B,PSFunPath=PSFunPath,IMPFunPath=IMPFunPath,px=px,testPurpose=testPurpose)
  }# end of DR

  return(trainRes)
}#end of function

