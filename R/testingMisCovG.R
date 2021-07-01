#' testingMisCovG
#'
#' Compute the empirical risk of cv or testing for kernel machies with missing covariates for hepatitis dataset
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
testingMisCovG=function(trainRes,testData){

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
  VTest=testData[,(px+1):(pTotal-3)]
  ZTest=cbind(XTest,VTest) # the covariates
  p=ncol(ZTest)
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
    R1Train=trainRes$R1
    R2Train=trainRes$R2
    XTrain=trainRes$X
    YTrain=trainRes$Y}# end of if kerMethod!="CC"

  if(testPurpose=="cv"){
    if(kerMethod=="DR"){
      VTrain=trainRes$V
      MissingCase_Train=trainRes$MissingCase
      XVparticalComplete_Train=trainRes$XV_particalComplete
      YparticalComplete_Train=trainRes$Y_particalComplete
      B=trainRes$impTime # imputation time
      ZfulTrain=trainRes$Zful # Covariates in imputation training data set, it ia used to compute kernel
      NTrain=nTrain+nTrain*B # sample size of training dataset after imputation
      
      piHatTest=rep(0,nTest)
      PSmodelTest1=PSFunTest(MissingCase_Train,XTrain,YTrain) 
      MissingCase_Test=ifelse((R1Test==0)*(R2Test==0)==1,0,1)
      piHatTest[MissingCase_Test==0]=PSmodelTest1(XTest[MissingCase_Test==0,],YTest[MissingCase_Test==0])
      
      
      XVparticalComplete_Test=cbind(XTest[MissingCase_Test==1,],VTest[MissingCase_Test==1,1])
      YparticalComplete_Test=YTest[MissingCase_Test==1]
      PSmodelTest2=PSFunTest(R2Train[MissingCase_Train==1],XVparticalComplete_Train,YparticalComplete_Train) 
      piHatTest[MissingCase_Test==1]=PSmodelTest2(XVparticalComplete_Test,YparticalComplete_Test) 
      
      
      IMPFunTestPath=trainRes$IMPFunPath
      IMPFunTest=readRDS(IMPFunTestPath)
      modelIMPTest_CC=IMPFunTest(XCCTrain,YCCTrain,VCCTrain) ## For R1=0 and R2=0
      modelIMPTest_PC=IMPFunTest(cbind(XCCTrain,VCCTrain[,1]),YCCTrain,VCCTrain[,2])
      pv=p-px
      df2=matrix(0,nTest*B,pv+1)
      
      for(i in 1:nTest){
        if(MissingCase_Test[i]==0){
          if(is.vector(XTest)){matIMPTest=modelIMPTest_CC(B,XTest[i],YTest[i])}else{
            matIMPTest=modelIMPTest_CC(B,XTest[i,],YTest[i])}
          df2[((i-1)*B+1):(i*B),]=as.matrix(matIMPTest)}else{
            if(is.vector(XTest)){matIMPTest=modelIMPTest_PC(B,c(XTest[i],VTest[i,1]),YTest[i])}else{
              matIMPTest=modelIMPTest_PC(B,c(XTest[i,],VTest[i,1]),YTest[i])}
            df2[((i-1)*B+1):(i*B),]=cbind(rep(VTest[i,1],B),as.matrix(matIMPTest))}# end of if else
      }# end of for 
      df2=as.data.frame(df2)
      names(df2) <-c(paste0("v",1:pv), "y")
      df2[,1][df2[,1]<0]=0.5
      df2[,2][df2[,2]<0]=0.5

      repFun=function(x){return(rep(x,each=B))}
      if(is.vector(XTest)){ImpDataTest=as.matrix(cbind(repFun(XTest),df2[,1:pv],repFun(RTest),df2$y,repFun(piHatTest)))}else{
      ImpDataTest=as.matrix(cbind(apply(XTest,2,repFun),df2[,1:pv],repFun(R1Test),repFun(R2Test),df2$y,repFun(piHatTest)))}

      ## Combine the original missing data and imputating data.
      #  The new data set will be dimension of [(n+1)*B]*pTotal
      NTest=nTest+nTest*B
      fulDataTest=matrix(0,NTest,pTotal+1)
      fulDataTest[1:nTest,]=cbind(ZTest,R1Test,R2Test,YTest,piHatTest)
      fulDataTest[(nTest+1):NTest,]=ImpDataTest

      ZfulTest=fulDataTest[,1:(pTotal-3)] # covariate in the new data
      R1fulTest=fulDataTest[,(pTotal-2)] # missing indicator in the new data
      R2fulTest=fulDataTest[,(pTotal-1)]
      YfulTest=fulDataTest[,pTotal] # response in the new data
      piHatfulTest=fulDataTest[,pTotal+1] # estimated propensity score in the new data

      ## Get vecMu and vecNu
      vecMu=rep(0,NTest)
      vecNu=rep(0,NTest)
      # mu=(B+1)*R_i/PiHat_i*I(Y_i==1) 1<=i<=n;
      #   =N/(n*B)*[R_i*(1-PiHat_i)/PiHat_i*I(Y_i==-1)+(1-R_i)*I(Y_i==1)], (n+1)<=i<=n*(B+1);
      vecMu[1:nTest]=(B+1)*(R1Test==1)*(R2Test==1)/piHatTest*as.numeric(YTest==1)
      vecMu[(nTest+1):NTest]=NTest/(nTest*B)*(((R1fulTest[(nTest+1):NTest]==1)*(R2fulTest[(nTest+1):NTest]==1))*((1-piHatfulTest[(nTest+1):NTest])/piHatfulTest[(nTest+1):NTest])*as.numeric(YfulTest[(nTest+1):NTest]==-1)
                                              +(1-((R1fulTest[(nTest+1):NTest]==1)*(R2fulTest[(nTest+1):NTest]==1)))*as.numeric(YfulTest[(nTest+1):NTest]==1))
      # nu=(B+1)*M_i/PiHat_i*I(Y_i==0), 1<=i<=n;
      #   =N/(n*B)*[M_i*(1-PiHat_i)/PiHat_i*I(Y_i==1)+(1-M_i)*I(Y_i==0)], (n+1)<=i<=n*(B+1).
      vecNu[1:nTest]=(B+1)*(R1Test==1)*(R2Test==1)/piHatTest*as.numeric(YTest==-1)
      vecNu[(nTest+1):NTest]=NTest/(nTest*B)*(((R1fulTest[(nTest+1):NTest]==1)*(R2fulTest[(nTest+1):NTest]==1))*((1-piHatfulTest[(nTest+1):NTest])/piHatfulTest[(nTest+1):NTest])*as.numeric(YfulTest[(nTest+1):NTest]==1)
                                              +(1-((R1fulTest[(nTest+1):NTest]==1)*(R2fulTest[(nTest+1):NTest]==1)))*as.numeric(YfulTest[(nTest+1):NTest]==0))

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
      Index=((R1Test==1)*(R2Test==1))==1
      ZCCTest=ZTest[Index,]# Complete covariates in testing data set
      XCCTest=ZTest[Index,1:px]
      YCCTest=YTest[Index]
      nCCTrain=nrow(ZCCTrain) # sample size of complete data in training data set
      nCCTest=nrow(ZCCTest) # sample size of complete data in testing data set

      # Get kernel matrix
      KTest=matrix(0,nCCTest,nCCTrain)
      for(i in 1:nCCTest){
        for(j in 1:nCCTrain)
        {KTest[i,j]=kerFun(ZCCTest[i,],ZCCTrain[j,])}}
      if(kerMethod=="WCC"){
        PSmodelTest=PSFunTest(R2Train,XTrain,YTrain)
        piHatTest=PSmodelTest(XTest,YTest)
        # Estimate f(x) in testing data set
        fHatWCC=KTest%*%diag(YCCTrain)%*%trainRes$hatAlphaWCC # f=sum_{j=1}^{N}alpha_j*Y_j*K(x,X[j])
        ## get the empirical risk
        piHatTestCC=piHatTest[R2Test==1]
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
      NTrain=nTrain+nTrain*B # sample size of training dataset after imputation
  
      # Get kernel matirx for test data set
      KTest=matrix(0,nTest,NTrain)
      for(i in 1:nTest){
        for(j in 1:NTrain)
        {KTest[i,j]=kerFun(ZTest[i,],ZfulTrain[j,])}
      }

      # Estimate f(x) in testing data set
      fHatDR=KTest%*%trainRes$hatAlphaDR # f=sum_{j=1}^{N}alpha_j*K(x,X[j])
      YHat=c(sign(fHatDR))}else{
        Index=((R1Test==1)*(R2Test==1))==1
        ZCCTrain=trainRes$ZCC # Complete covariates in training data set
        ## complete data in the testing data
        ZCCTest=ZTest[Index==1,]# Complete covariates in testing data set
        XCCTest=ZTest[Index==1,1:px]
        YCCTest=YTest[Index==1]
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
