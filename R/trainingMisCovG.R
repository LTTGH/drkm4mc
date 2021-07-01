#' trainingMisCovG
#'
#' training function for kernel machies with missing covariates for hepatitis dataset
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

trainingMisCovG=function(misCovDat,px,kerType,kerMethod,
                         lambda,sigma,PSFunPath,IMPFunPath,B,testPurpose){
  pTotal=dim(misCovDat)[2] # dimension of (X,V,R,Y1,Y)
  misCovDat=as.matrix(misCovDat)
  #colnames(misCovDat)=NULL
  X=misCovDat[,1:px] # the fully observed covariates
  V=misCovDat[,(px+1):(pTotal-4)] # the popential missing covariates
  Z=cbind(X,V) # the covariates
  p=ncol(Z) # dimension of covariates
  R1=misCovDat[,(pTotal-3)] 
  R2=misCovDat[,(pTotal-2)]   #  missing indicator
  Y=misCovDat[,pTotal] # the response (-1 and 1)
  n=nrow(misCovDat) # sample size of the traing data
  
  ## Get kernel matirx
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel
  
  
  YCC=Y[((R1==1)*(R2==1))==1] # Y corresponding complete cases
  ZCC=Z[((R1==1)*(R2==1))==1,] # covariates corresponding complete cases
  XCC=Z[((R1==1)*(R2==1))==1,1:px] # Y corresponding complete cases
  VCC=Z[((R1==1)*(R2==1))==1,(px+1):p] # covariates corresponding complete cases
  YCCMatrix=diag(YCC)
  nCC=length(YCC)  # sample size of the traing data
  
  C=1/(2*nCC*lambda) # the parameter in the minimization process.
  
  # Get kernel matrix for WCC and CC
  if(kerMethod!="DR"){
    K=kernelMatrix(kerFun,ZCC)}# end of if
  
  ## CC method
  if(kerMethod=="CC"){
    ## Constrainted optmiazed program for CC
    # L=(1/2)*alpha'*D*K*D*alpha-alpha'1, subgect to 0<alpha<C*1,
    # where D=diag(YCC), 1 is a vector, 0 is a vector.
    res=quadprog(YCCMatrix%*%(K+diag(0.0001,nCC))%*%YCCMatrix,rep(-1,nCC),A=NULL,b=NULL,Aeq=rep(0,nCC),beq=0,
                 lb=rep(0,nCC),ub=rep(C,nCC))
    # Get the estimator of alpha
    hatAlphaCC=res$x
    trainRes=list(hatAlphaCC=hatAlphaCC,ZCC=ZCC,YCC=YCC,R1=R1,R2=R2,trainSampleSize=n,
                  sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                  px=px,testPurpose=testPurpose)
  } # end of CC
  
  
  ## Estimate the propensity score model for WCC and DR
  if(kerMethod!="CC"){
    PSFun=readRDS(PSFunPath)
    ## WCC method
    if(kerMethod=="WCC"){
      ## Constrainted optmiazed program for WCC
      # L=alpha'1-(1/2)*alpha'*D*K*D*alpha, subgect to 0<alpha<C*W*1,
      # where D=diag(YCC), W=diago(piHat), 1 is a vector, 0 is a vector.
      PSmodel=PSFun(R2,X,Y) 
      piHat=PSmodel(X,Y)
      W=diag(1/piHat[R2==1]) # weight matrix diag(1/pi_1,...1/pi_nCC)
      res=quadprog(YCCMatrix%*%(K)%*%YCCMatrix,rep(-1,nCC),A=NULL,b=NULL,Aeq=rep(0,nCC),beq=0,
                   lb=rep(0,nCC),ub=W%*%rep(C,nCC))
      # Get the estimator of alpha
      hatAlphaWCC=res$x
      trainRes=list(hatAlphaWCC=hatAlphaWCC,ZCC=ZCC,YCC=YCC,R1=R1,R2=R2,X=X,Y=Y,trainSampleSize=n,
                    sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                    PSFunPath=PSFunPath,px=px,testPurpose=testPurpose)
    } # end of WCC
    
    ## DR method
    if(kerMethod=="DR"){
      
      ## R1=0 and R2=0
      piHat=rep(0,n)
      MissingCase1=ifelse((R1==0)*(R2==0)==1,0,1)
      PSmodel1=PSFun(MissingCase1,X,Y) 
      piHat[MissingCase1==0]=PSmodel1(X[MissingCase1==0,],Y[MissingCase1==0])
      
      ## R1=1 and R2=0
      XV_particalComplete=cbind(X[MissingCase1==1,],V[MissingCase1==1,1])
      Y_particalComplete=Y[MissingCase1==1]
      PSmodel2=PSFun(R2[MissingCase1==1],XV_particalComplete,Y_particalComplete) 
      piHat[MissingCase1==1]=PSmodel2(XV_particalComplete,Y_particalComplete) 
      
      ## Imputation process
      IMPFun=readRDS(IMPFunPath)
      modelIMP_CC=IMPFun(XCC,YCC,VCC) ## For R1=0 and R2=0
      modelIMP_PC=IMPFun(cbind(XCC,VCC[,1]),YCC,VCC[,2])
      
      pv=p-px
      df2=matrix(0,n*B,pv+1)
      for(i in 1:n){
        if(MissingCase1[i]==0){
          if(is.vector(X)){matIMP=modelIMP_CC(B,X[i],Y[i])}else{
            matIMP=modelIMP_CC(B,X[i,],Y[i])}
          df2[((i-1)*B+1):(i*B),]=as.matrix(matIMP)}else{
            if(is.vector(X)){matIMP=modelIMP_PC(B,c(X[i],V[i,1]),Y[i])}else{
              matIMP=modelIMP_PC(B,c(X[i,],V[i,1]),Y[i])}
            df2[((i-1)*B+1):(i*B),]=cbind(rep(V[i,1],B),as.matrix(matIMP))}# end of if else
      }# end of for 
      df2=as.data.frame(df2)
      names(df2) <-c(paste0("v",1:pv), "y")
      df2[,1][df2[,1]<0]=0.5
      df2[,2][df2[,2]<0]=0.5
      
      repFun=function(x){return(rep(x,each=B))}
      if(is.vector(X)){ImpData=as.matrix(cbind(repFun(X),df2[,1:pv],repFun(R1),repFun(R2),df2$y,repFun(piHat)))}else{
        ImpData=as.matrix(cbind(apply(X,2,repFun),df2[,1:pv],repFun(R1),repFun(R2),df2$y,repFun(piHat)))}
      
      
      ## Combine the original missing data and imputating data.
      #  The new data set will be dimension of [(n+1)*B]*4
      N=n+n*B
      fulData=matrix(0,N,pTotal)
      fulData[1:n,]=cbind(Z,R1,R2,Y,piHat)
      fulData[(n+1):N,]=ImpData
      
      Zful=fulData[,1:(pTotal-4)] # covariate in the new data
      R1ful=fulData[,(pTotal-3)] # missing indicator in the new data
      R2ful=fulData[,(pTotal-2)] # missing indicator in the new data
      Yful=fulData[,(pTotal-1)] # response in the new data
      piHatful=fulData[,pTotal] # estimated propensity score in the new data
      
      
      
      ## Get vecMu and vecNu
      vecMu=rep(0,N)
      vecNu=rep(0,N)
      # mu=(B+1)*R_i/PiHat_i*I(Y_i==1) 1<=i<=n;
      #   =N/(n*B)*[R_i*(1-PiHat_i)/PiHat_i*I(Y_i==-1)+(1-M_i)*I(Y_i==1)], (n+1)<=i<=n*(B+1);
      vecMu[1:n]=(B+1)*(R1==1)*(R2==1)/piHat*as.numeric(Y==1)
      vecMu[(n+1):N]=N/(n*B)*(((R1ful[(n+1):N]==1)*(R2ful[(n+1):N]==1))*((1-piHatful[(n+1):N])/piHatful[(n+1):N])*as.numeric(Yful[(n+1):N]==-1)
                              +(1-((R1ful[(n+1):N]==1)*(R2ful[(n+1):N]==1)))*as.numeric(Yful[(n+1):N]==1))
      
      # nu=(B+1)*M_i/PiHat_i*I(Y_i==0), 1<=i<=n;
      #   =N/(n*B)*[M_i*(1-PiHat_i)/PiHat_i*I(Y_i==1)+(1-M_i)*I(Y_i==-1)], (n+1)<=i<=n*(B+1).
      vecNu[1:n]=(B+1)*(R1==1)*(R2==1)/piHat*as.numeric(Y==-1)
      vecNu[(n+1):N]=N/(n*B)*(((R1ful[(n+1):N]==1)*(R2ful[(n+1):N]==1))*((1-piHatful[(n+1):N])/piHatful[(n+1):N])*as.numeric(Yful[(n+1):N]==1)
                              +(1-((R1ful[(n+1):N]==1)*(R2ful[(n+1):N]==1)))*as.numeric(Yful[(n+1):N]==-1))
      
      
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
      
      trainRes=list(hatAlphaDR=hatAlphaDR,Zful=Zful,ZCC=ZCC,YCC=YCC,R1=R1,R2=R2,X=X,Y=Y,V=V,
                    MissingCase=MissingCase1,XV_particalComplete=XV_particalComplete,
                    Y_particalComplete=Y_particalComplete,
                    trainSampleSize=n,sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                    impTime=B,PSFunPath=PSFunPath,IMPFunPath=IMPFunPath,px=px,testPurpose=testPurpose)
    }# end of DR
  } # end of !CC  
    return(trainRes)
}#end of function

