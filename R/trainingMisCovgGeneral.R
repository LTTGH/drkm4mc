#' trainingMisCovGeneral 
#' 
#' training function for kernel machies with non-monotone type of missing covariates
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

trainingMisCovGeneral=function(misCovDat,px,kerType,kerMethod,
                         lambda,sigma,PSFunPath,IMPFunPath,B,testPurpose){
  pTotal=dim(misCovDat)[2] # dimension of (X,V,R,Y1,Y)
  misCovDat=as.matrix(misCovDat)
  #colnames(misCovDat)=NULL
  X=misCovDat[,1:px] # the fully observed covariates
  V=misCovDat[,(px+1):(pTotal-5)] # the popential missing covariates
  Z=cbind(X,V) # the covariates
  p=ncol(Z) # dimension of covariates
  R=misCovDat[,(pTotal-4)] 
  R1=misCovDat[,(pTotal-3)] 
  R2=misCovDat[,(pTotal-2)]   #  missing indicator
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
      PSmodel=PSFun(R,X,Y) 
      piHat=PSmodel(X,Y)
      W=diag(1/piHat[R==1]) # weight matrix diag(1/pi_1,...1/pi_nCC)
      res=quadprog(YCCMatrix%*%(K)%*%YCCMatrix,rep(-1,nCC),A=NULL,b=NULL,Aeq=rep(0,nCC),beq=0,
                   lb=rep(0,nCC),ub=W%*%rep(C,nCC))
      # Get the estimator of alpha
      hatAlphaWCC=res$x
      trainRes=list(hatAlphaWCC=hatAlphaWCC,ZCC=ZCC,YCC=YCC,R=R,R1=R1,R2=R2,X=X,Y=Y,trainSampleSize=n,
                    sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                    PSFunPath=PSFunPath,px=px,testPurpose=testPurpose)
    } # end of WCC
    
    ## DR method
    if(kerMethod=="DR"){
      
      ## X2 is missing
      piHat=rep(0,n)
      PSmodel0=PSFun(R,X,Y) 
      index0=(R==0)-(R==0)*(R1==0)-(R==0)*(R2==0)##index for X2 is missng
      piHat[index0==1]=PSmodel0(X[index0==1,],Y[index0==1])
      
      ## X21 is missing (x16,x17,x19,x20)
      PSmodel1=PSFun(R1[R1>=0],cbind(X[R1>=0,],V[R1>=0,c(1,4)]),Y[R1>=0]) 
      piHat[R1==0]=PSmodel1(cbind(X[R1==0,],V[R1==0,c(1,4)]),Y[R1==0])
      pi1=PSmodel1(cbind(X[R1==1,],V[R1==1,c(1,4)]),Y[R1==1])
      
      ## X22 is missing (x15,x18)
      PSmodel2=PSFun(R2[R2>=0],cbind(X[R2>=0,],V[R2>=0,c(2,3,5,6)]),Y[R2>=0]) 
      piHat[R2==0]=PSmodel2(cbind(X[R2==0,],V[R2==0,c(2,3,5,6)]),Y[R2==0])
      pi2=PSmodel2(cbind(X[R2==1,],V[R2==1,c(2,3,5,6)]),Y[R2==1])
      
      
      ## Imputation process
      IMPFun=readRDS(IMPFunPath)
      modelIMP_CC=IMPFun(XCC,YCC,VCC) ## For R1=0 and R2=0
      modelIMP_PC1=IMPFun(cbind(XCC,VCC[,c(1,4)]),YCC,VCC[,c(2,3,5,6)])
      modelIMP_PC2=IMPFun(cbind(XCC,VCC[,c(2,3,5,6)]),YCC,VCC[,c(1,4)])
      
      pv=p-px
      df2=matrix(0,n*B,pv+1)
      df3=matrix(0,n*B,pv+1)
      for(i in 1:n){
        if(index0[i]==1){matIMP=modelIMP_CC(B,X[i,],Y[i])
                          df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=as.matrix(matIMP)}
        if(R1[i]==0){matIMP=modelIMP_PC1(B,c(X[i,],V[i,c(1,4)]),Y[i])
                     matIMP=as.matrix(matIMP)
                         df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=cbind(rep(V[i,1],B),matIMP[,1:2],
                                                   rep(V[i,4],B),matIMP[,3:4],matIMP[,5])}
        if(R1[i]==1){matIMP=modelIMP_PC1(B,c(X[i,],V[i,c(1,4)]),Y[i])
                     matIMP=as.matrix(matIMP)
                     df2[((i-1)*B+1):(i*B),]=cbind(rep(V[i,1],B),matIMP[,1:2],
                                      rep(V[i,4],B),matIMP[,3:4],matIMP[,5])}
        
        if(R2[i]==1){matIMP=modelIMP_PC2(B,c(X[i,],V[i,c(2,3,5,6)]),Y[i])
                     matIMP=as.matrix(matIMP)
                     df3[((i-1)*B+1):(i*B),]=cbind(matIMP[,1],rep(V[i,2],B),rep(V[i,3],B),matIMP[,2],
                                      rep(V[i,5],B),rep(V[i,6],B),matIMP[,3])}
        
        if(R2[i]==0){matIMP=modelIMP_PC2(B,c(X[i,],V[i,c(2,3,5,6)]),Y[i])
                     matIMP=as.matrix(matIMP)
                     df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=c(matIMP[,1],rep(V[i,2],B),rep(V[i,3],B),matIMP[,2],
                                                   rep(V[i,5],B),rep(V[i,6],B),matIMP[,3])}
      } # end of for
      df2=as.data.frame(df2)
      names(df2) <-c(paste0("v",1:pv), "y")
      df3=as.data.frame(df3)
      names(df3) <-c(paste0("v",1:pv), "y")
      piHat1=piHat
      piHat1[R==1]=pi1
      piHat2=piHat
      piHat2[R==1]=pi2
      
      
     repFun=function(x){return(rep(x,each=B))}
     RImpute=repFun(R)
     X_N=rbind(apply(X,2,repFun),apply(X[R==1,],2,repFun))
     df_N=rbind(df2[,1:pv],df3[RImpute==1,1:pv])
     R_N=c(repFun(R),repFun(R)[RImpute==1])
     R1_N=c(repFun(R1),repFun(R1)[RImpute==1])
     R2_N=c(repFun(R2),repFun(R2)[RImpute==1])
     Y_N=c(df2$y,df3$y[RImpute==1])
     piHatN=c(repFun(piHat1),repFun(piHat2)[RImpute==1])
     ImpData=as.matrix(cbind(X_N,df_N,R_N,R1_N,R2_N,Y_N,piHatN))
      
      ## Combine the original missing data and imputating data.
      #  The new data set will be dimension of [(n+1)*B]*4
      N=n+dim(ImpData)[1]
      fulData=matrix(0,N,pTotal);
      fulData[1:n,]=cbind(Z,R,R1,R2,Y,piHat1)
      fulData[(n+1):N,]=ImpData
      
      Zful=fulData[,1:(pTotal-5)] # covariate in the new data
      Rful=fulData[,(pTotal-4)]
      R1ful=fulData[,(pTotal-3)] # missing indicator in the new data
      R2ful=fulData[,(pTotal-2)] # missing indicator in the new data
      Yful=fulData[,(pTotal-1)] # response in the new data
      piHatful=fulData[,pTotal] # estimated propensity score in the new data
      
      ## Get vecMu and vecNu
      vecMu=rep(0,N)
      vecNu=rep(0,N)
      # mu=(B+1)*R_i/PiHat_i*I(Y_i==1) 1<=i<=n;
      #   =N/(n*B)*[R_i*(1-PiHat_i)/PiHat_i*I(Y_i==-1)+(1-M_i)*I(Y_i==1)], (n+1)<=i<=n*(B+1);
      vecMu[1:n]=(B+1)*(R==1)/(2*piHat1)*as.numeric(Y==1)+(B+1)*(R==1)/(2*piHat2)*as.numeric(Y==1)
      vecMu[(n+1):N]=N/(n*B)*((Rful[(n+1):N]==1)*((1-piHatful[(n+1):N])/(2*piHatful[(n+1):N]))*as.numeric(Yful[(n+1):N]==-1)+
                              (1-(Rful[(n+1):N]==1))*as.numeric(Yful[(n+1):N]==1))
      
      # nu=(B+1)*M_i/PiHat_i*I(Y_i==0), 1<=i<=n;
      #   =N/(n*B)*[M_i*(1-PiHat_i)/PiHat_i*I(Y_i==1)+(1-M_i)*I(Y_i==-1)], (n+1)<=i<=n*(B+1).
      vecNu[1:n]=(B+1)*(R==1)/(2*piHat1)*as.numeric(Y==-1)+(B+1)*(R==1)/(2*piHat2)*as.numeric(Y==-1)
      vecNu[(n+1):N]=N/(n*B)*((Rful[(n+1):N]==1)*((1-piHatful[(n+1):N])/(2*piHatful[(n+1):N]))*as.numeric(Yful[(n+1):N]==1)+
                              (1-(Rful[(n+1):N]==1))*as.numeric(Yful[(n+1):N]==-1))
      
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
      
      trainRes=list(hatAlphaDR=hatAlphaDR,Zful=Zful,ZCC=ZCC,YCC=YCC,R=R,R1=R1,R2=R2,X=X,Y=Y,V=V,
                    trainSampleSize=n,sigma=sigma,kerMethod=kerMethod,kerType=kerType,
                    impTime=B,PSFunPath=PSFunPath,IMPFunPath=IMPFunPath,px=px,testPurpose=testPurpose)
    }# end of DR
  } # end of !CC  
  return(trainRes)
}#end of function

