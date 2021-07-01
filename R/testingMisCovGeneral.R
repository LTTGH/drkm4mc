#' testingMisCovGeneral
#'
#' Compute the empirical risk of cv or testing for kernel machies with non-monotone type of missing covariates
#'
#' @param trainRes  Results results from trainging results
#' @param testData testing data set
#' @return YHat the empirical risk of the test
#'
#' @import e1071
#' @importFrom kernlab rbfdot vanilladot
#' @import stats
#' @export
#'
#'
testingMisCovGeneral=function(trainRes,testData){

  ##Get the parameter from training function.
  nTrain=trainRes$trainSampleSize # sample size of original data set
  kerMethod=trainRes$kerMethod
  sigma=trainRes$sigma
  kerType=trainRes$kerType
  testPurpose=trainRes$testPurpose
  px=trainRes$px #dimension of fully observed covariate

  ZCCTrain=trainRes$ZCC # Complete covariates in training data set
  YCCTrain=trainRes$YCC # Complete covariates in training data set
  XCCTrain=ZCCTrain[,1:px] # Y corresponding complete cases
  VCCTrain=ZCCTrain[,-c(1:px)] # covariates corresponding complete cases

  ## data structure is (X,V,R,Y)
  testData=as.matrix(testData)

  nTest=nrow(testData) # Sample size of testing data
  pTotal=dim(testData)[2] # dimensions of (X,V,R,Y)
  XTest=testData[,1:px]
  VTest=testData[,(px+1):(pTotal-4)]
  ZTest=cbind(XTest,VTest) # the covariates
  p=ncol(ZTest)
  RTest=testData[,(pTotal-3)] 
  R1Test=testData[,(pTotal-2)] 
  R2Test=testData[,(pTotal-1)]   #  missing indicator
  YTest=testData[,pTotal] # the response (-1 and 1)

  ## Get kernel matirx
  if (kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel

  ## Estimate propensity for WCC and DR
  if(kerMethod!="CC"){
    PSFunTestPath=trainRes$PSFunPath
    PSFunTest=readRDS(PSFunTestPath)
    RTrain=trainRes$R
    R1Train=trainRes$R1
    R2Train=trainRes$R2
    XTrain=trainRes$X
    YTrain=trainRes$Y}# end of if kerMethod!="CC"

  if(testPurpose=="cv"){
    if(kerMethod=="DR"){
      VTrain=trainRes$V
      B=trainRes$impTime # imputation time
      ZfulTrain=trainRes$Zful # Covariates in imputation training data set, it ia used to compute kernel
      NTrain=dim(ZfulTrain)[1] # sample size of training dataset after imputation
      
      piHatTest=rep(0,nTest)
      PSmodel0Test=PSFunTest(RTrain,XTrain,YTrain) 
      index0Test=(RTest==0)-(RTest==0)*(R1Test==0)-(RTest==0)*(R2Test==0)##index for X2 is missng
      piHatTest[index0Test==1]=PSmodel0Test(XTest[index0Test==1,],YTest[index0Test==1])
      
      ## X21 is missing (x16,x17,x19,x20)
      PSmodel1Test=PSFunTest(R1Train[R1Train>=0],cbind(XTrain[R1Train>=0,],VTrain[R1Train>=0,c(1,4)]),YTrain[R1Train>=0]) 
      piHatTest[R1Test==0]=PSmodel1Test(cbind(XTest[R1Test==0,],VTest[R1Test==0,c(1,4)]),YTest[R1Test==0])
      pi1Test=PSmodel1Test(cbind(XTest[R1Test==1,],VTest[R1Test==1,c(1,4)]),YTest[R1Test==1])
      
      ## X22 is missing (x15,x18)
      PSmodel2Test=PSFunTest(R2Train[R2Train>=0],cbind(XTrain[R2Train>=0,],VTrain[R2Train>=0,c(2,3,5,6)]),YTrain[R2Train>=0]) 
      piHatTest[R2Test==0]=PSmodel2Test(cbind(XTest[R2Test==0,],VTest[R2Test==0,c(2,3,5,6)]),YTest[R2Test==0])
      pi2Test=PSmodel2Test(cbind(XTest[R2Test==1,],VTest[R2Test==1,c(2,3,5,6)]),YTest[R2Test==1])
      
      
      IMPFunTestPath=trainRes$IMPFunPath
      IMPFunTest=readRDS(IMPFunTestPath)
      modelIMPTest_CC=IMPFunTest(XCCTrain,YCCTrain,VCCTrain) ## For R1=0 and R2=0
      modelIMPTest_PC1=IMPFunTest(cbind(XCCTrain,VCCTrain[,c(1,4)]),YCCTrain,VCCTrain[,c(2,3,5,6)])
      modelIMPTest_PC2=IMPFunTest(cbind(XCCTrain,VCCTrain[,c(2,3,5,6)]),YCCTrain,VCCTrain[,c(1,4)])
      
      pv=p-px
      df2=matrix(0,nTest*B,pv+1)
      df3=matrix(0,nTest*B,pv+1)
      
      for(i in 1:nTest){
        if(index0Test[i]==1){matIMPTest=modelIMPTest_CC(B,XTest[i,],YTest[i])
        df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=as.matrix(matIMPTest)}
        if(R1Test[i]==0){matIMPTest=modelIMPTest_PC1(B,c(XTest[i,],VTest[i,c(1,4)]),YTest[i])
        matIMPTest=as.matrix(matIMPTest)
        df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=cbind(rep(VTest[i,1],B),matIMPTest[,1:2],
                                                              rep(VTest[i,4],B),matIMPTest[,3:4],matIMPTest[,5])}
        if(R1Test[i]==1){matIMPTest=modelIMPTest_PC1(B,c(XTest[i,],VTest[i,c(1,4)]),YTest[i])
        matIMPTest=as.matrix(matIMPTest)
        df2[((i-1)*B+1):(i*B),]=cbind(rep(VTest[i,1],B),matIMPTest[,1:2],
                                      rep(VTest[i,4],B),matIMPTest[,3:4],matIMPTest[,5])}
        
        if(R2Test[i]==1){matIMPTest=modelIMPTest_PC2(B,c(XTest[i,],VTest[i,c(2,3,5,6)]),YTest[i])
        matIMPTest=as.matrix(matIMPTest)
        df3[((i-1)*B+1):(i*B),]=cbind(matIMPTest[,1],rep(VTest[i,2],B),rep(VTest[i,3],B),matIMPTest[,2],
                                      rep(VTest[i,5],B),rep(VTest[i,6],B),matIMPTest[,3])}
        
        if(R2Test[i]==0){matIMPTest=modelIMPTest_PC2(B,c(XTest[i,],VTest[i,c(2,3,5,6)]),YTest[i])
        matIMPTest=as.matrix(matIMPTest)
        df3[((i-1)*B+1):(i*B),]=df2[((i-1)*B+1):(i*B),]=c(matIMPTest[,1],rep(VTest[i,2],B),rep(VTest[i,3],B),matIMPTest[,2],
                                                          rep(VTest[i,5],B),rep(VTest[i,6],B),matIMPTest[,3])}
      } # end of for
      df2=as.data.frame(df2)
      names(df2) <-c(paste0("v",1:pv), "y")
      df3=as.data.frame(df3)
      names(df3) <-c(paste0("v",1:pv), "y")
      
      piHatTest1=piHatTest
      piHatTest1[RTest==1]=pi1Test
      piHatTest2=piHatTest
      piHatTest2[RTest==1]=pi2Test
      

      repFun=function(x){return(rep(x,each=B))}
      RImputeTest=repFun(RTest)
      X_NTest=rbind(apply(XTest,2,repFun),apply(XTest[RTest==1,],2,repFun))
      df_NTest=rbind(df2[,1:pv],df3[RImputeTest==1,1:pv])
      R_NTest=c(repFun(RTest),repFun(RTest)[RImputeTest==1])
      R1_NTest=c(repFun(R1Test),repFun(R1Test)[RImputeTest==1])
      R2_NTest=c(repFun(R2Test),repFun(R2Test)[RImputeTest==1])
      Y_NTest=c(df2$y,df3$y[RImputeTest==1])
      piHatNTest=c(repFun(piHatTest1),repFun(piHatTest2)[RImputeTest==1])
      ImpDataTest=as.matrix(cbind(X_NTest,df_NTest,R_NTest,R1_NTest,R2_NTest,Y_NTest,piHatNTest))
      ## Combine the original missing data and imputating data.
      #  The new data set will be dimension of [(n+1)*B]*pTotal
      NTest=nTest+dim(ImpDataTest)[1]
      fulDataTest=matrix(0,NTest,pTotal+1)
      fulDataTest[1:nTest,]=cbind(ZTest,RTest,R1Test,R2Test,YTest,piHatTest)
      fulDataTest[(nTest+1):NTest,]=ImpDataTest

      ZfulTest=fulDataTest[,1:(pTotal-4)] # covariate in the new data
      RfulTest=fulDataTest[,(pTotal-3)]
      R1fulTest=fulDataTest[,(pTotal-2)] # missing indicator in the new data
      R2fulTest=fulDataTest[,(pTotal-1)]
      YfulTest=fulDataTest[,pTotal] # response in the new data
      piHatfulTest=fulDataTest[,pTotal+1] # estimated propensity score in the new data
      

      ## Get vecMu and vecNu
      vecMu=rep(0,NTest)
      vecNu=rep(0,NTest)
      # mu=(B+1)*R_i/PiHat_i*I(Y_i==1) 1<=i<=n;
      #   =N/(n*B)*[R_i*(1-PiHat_i)/PiHat_i*I(Y_i==-1)+(1-R_i)*I(Y_i==1)], (n+1)<=i<=n*(B+1);
      vecMu[1:nTest]=(B+1)*(RTest==1)/(2*piHatTest1)*as.numeric(YTest==1)+
                     (B+1)*(RTest==1)/(2*piHatTest2)*as.numeric(YTest==1)
      vecMu[(nTest+1):NTest]=NTest/(nTest*B)*((RfulTest[(nTest+1):NTest]==1)*((1-piHatfulTest[(nTest+1):NTest])/(2*piHatfulTest[(nTest+1):NTest]))*as.numeric(YfulTest[(nTest+1):NTest]==-1)
                                              +(1-(RfulTest[(nTest+1):NTest]==1))*as.numeric(YfulTest[(nTest+1):NTest]==1))
      # nu=(B+1)*M_i/PiHat_i*I(Y_i==0), 1<=i<=n;
      #   =N/(n*B)*[M_i*(1-PiHat_i)/PiHat_i*I(Y_i==1)+(1-M_i)*I(Y_i==0)], (n+1)<=i<=n*(B+1).
      vecNu[1:nTest]=(B+1)*(RTest==1)/(2*piHatTest1)*as.numeric(YTest==-1)+
                     (B+1)*(RTest==1)/(2*piHatTest2)*as.numeric(YTest==-1)
      vecNu[(nTest+1):NTest]=NTest/(nTest*B)*((RfulTest[(nTest+1):NTest]==1)*((1-piHatfulTest[(nTest+1):NTest])/(2*piHatfulTest[(nTest+1):NTest]))*as.numeric(YfulTest[(nTest+1):NTest]==1)
                                              +(1-(RfulTest[(nTest+1):NTest]==1))*as.numeric(YfulTest[(nTest+1):NTest]==-1))

      # Get kernel matirx for test data set
      KTest=matrix(0,NTest,NTrain)
      for(i in 1:NTest){
        for(j in 1:NTrain)
        {KTest[i,j]=kerFun(ZfulTest[i,],ZfulTrain[j,])} # end of j
      } # end of i

      # Estimate f(x) in testing data set
      fHatDR=KTest%*%trainRes$hatAlphaDR # f=sum_{j=1}^{N}alpha_j*K(x,X[j])

      # Compute the convex surrogate loss
      phi_f=apply(as.data.frame(cbind(rep(0,NTest),1-fHatDR)),1,max)
      phi_mf=apply(as.data.frame(cbind(rep(0,NTest),1+fHatDR)),1,max)

      # Compute the empirical risk
      DR_empRiskCV=sum(vecMu%*%phi_f+vecNu%*%phi_mf)/NTrain
      YHat=c(DR_empRiskCV,rep(0,nTest-1))
    } else {

      ## complete data in the testing data
      ZCCTest=ZTest[RTest==1,]# Complete covariates in testing data set
      XCCTest=ZTest[RTest==1,1:px]
      YCCTest=YTest[RTest==1]
      nCCTrain=nrow(ZCCTrain) # sample size of complete data in training data set
      nCCTest=nrow(ZCCTest) # sample size of complete data in testing data set

      # Get kernel matrix
      KTest=matrix(0,nCCTest,nCCTrain)
      for(i in 1:nCCTest){
        for(j in 1:nCCTrain)
        {KTest[i,j]=kerFun(ZCCTest[i,],ZCCTrain[j,])}}
      if(kerMethod=="WCC"){
        PSmodelTest=PSFunTest(RTrain,XTrain,YTrain)
        piHatTest=PSmodelTest(XTest,YTest)
        # Estimate f(x) in testing data set
        fHatWCC=KTest%*%diag(YCCTrain)%*%trainRes$hatAlphaWCC # f=sum_{j=1}^{N}alpha_j*Y_j*K(x,X[j])
        ## get the empirical risk
        piHatTestCC=piHatTest[RTest==1]
        hingLossWCC=apply(as.data.frame(cbind(rep(0,nCCTest),YCCTest*(1-fHatWCC))),1,max) # hingle loss
        WCC_empRiskCV=mean((1/piHatTestCC)*hingLossWCC)
        YHat=c(WCC_empRiskCV,rep(0,nTest-1))
      }else{
        ## CC method for CV
        # Estimate f(x) in testing data set
        fHatCC=KTest%*%diag(YCCTrain)%*%trainRes$hatAlphaCC # f=sum_{j=1}^{N}Y_jalpha_j*K(x,X[j])}}
        hingLossCC=apply(as.data.frame(cbind(rep(0,nCCTest),YCCTest*(1-fHatCC))),1,max) # hingle loss
        CC_empRiskCV=mean(hingLossCC)
        YHat=c(CC_empRiskCV,rep(0,nTest-1))} # end of CC CV
    } # end of WCC and CC CV
  } else {
    ## Testing purpose
    if(kerMethod=="DR"){
      B=trainRes$impTime # imputation time
      ZfulTrain=trainRes$Zful # Covariates in imputation training data set, it ia used to compute kernel
      NTrain=dim(trainRes$Zful)[1] # sample size of training dataset after imputation
  
      # Get kernel matirx for test data set
      KTest=matrix(0,nTest,NTrain)
      for(i in 1:nTest){
        for(j in 1:NTrain)
        {KTest[i,j]=kerFun(ZTest[i,],ZfulTrain[j,])}
      }

      # Estimate f(x) in testing data set
      fHatDR=KTest%*%trainRes$hatAlphaDR # f=sum_{j=1}^{N}alpha_j*K(x,X[j])
      YHat=c(sign(fHatDR))}else{
        ZCCTrain=trainRes$ZCC # Complete covariates in training data set
        ## complete data in the testing data
        ZCCTest=ZTest[RTest==1,]# Complete covariates in testing data set
        XCCTest=ZTest[RTest==1,1:px]
        YCCTest=YTest[RTest==1]
        nCCTrain=nrow(ZCCTrain) # sample size of complete data in training data set

        # Get kernel matrix
        KTest=matrix(0,nTest,nCCTrain)
        for(i in 1:nTest){
          for(j in 1:nCCTrain)
          {KTest[i,j]=kerFun(ZTest[i,],ZCCTrain[j,])}}

        if(kerMethod=="WCC"){
          # Estimate f(x) in testing data set
          fHatWCC=KTest%*%diag(YCCTrain)%*%trainRes$hatAlphaWCC # f=sum_{j=1}^{N}alpha_j*K(x,X[j])
          YHat=c(sign(fHatWCC))
        }else{if(kerMethod=="CC"){
          fHatCC=KTest%*%diag(YCCTrain)%*%trainRes$hatAlphaCC # f=sum_{j=1}^{N}alpha_j*K(x,X[j])
          YHat=c(sign(fHatCC))
        }# end of 3rd ifelse of kerMethod CC
        }# end of 2nd ifelse of kerMethod WCC
      }# end of 1st ifelse of kerMethod DR
  } # end of testing purpose
  return(YHat)} # end of function
