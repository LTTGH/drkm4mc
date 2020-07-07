# drkm4mc
[R/drkm4mc](https://github.com/LTTGH/drkm4mc) is an [R](https:/www.r-project.org) package. 
This package is used to construct kernel machines with missing covariates for the problem of classification. A doubly robust kernel machine and a weighted-complete-case kernel machies are given.

#### Installation
You can install it from its [GitHub repository](https://github.com/LTTGH/drkm4mc). You first need to install the [devtools](https://github.com/hadley/devtools) package.

install R/drkm4mc using the `install_github` function in the [devtools](https://github.com/hadley/devtools) package.

```r
library(devtools)
install_github("LTTGH/drkm4mc")
```

#### Example use
The example is the veterbral dataset, this dataset is available from the UCI machine learning repository (www.ics.uci.edu/˜mlearn/MLRepository.html) and used in the paper
kernel machines for missing covariates. The dataset was split into two sets, 
the first dataset is used to train the doubly-robust kernel machine with samplesize 200.
The remaining dataset is used to test the  doubly-robust kernel machine with samplesize 110. 


```{r}
library(drkm4mc)
vData=VertebralData
trainSampleSize=200
trainLable=sample(c(1:310),trainSampleSize,replace = F)
trainData=vData[trainLable,] 
testLabel=setdiff(c(1:310),trainLable)
testData=vData[testLabel,]
```

First generate the missing mechanism in the training dataset. Assume that 
x4 and x5 are the covariates that are popentially missing.
```{r}
alpha=2
beta=-11
stdx1=(trainData$x1-mean(trainData$x1))/sd(trainData$x1)
stdx2=(trainData$x2-mean(trainData$x2))/sd(trainData$x2)
PScov=(stdx1+stdx2)/2
probR=exp(alpha+beta*(PScov*trainData$Y))/(1+exp(alpha+beta*(PScov*trainData$Y)))
R=rbinom(trainSampleSize,1,probR)
```

Vertebral data with missing covariates is
```{r}
misCovDat=as.data.frame(cbind(trainData$x1,trainData$x2,trainData$x3,trainData$x6,
                              trainData$x4,trainData$x5,R,trainData$Y,trainData$Y1))
```


Set the pathes of imputation model and propensity score model 
```{r}
psPath=url(("https://github.com/LTTGH/drkm4mc/raw/master/ModelPS1.RDS"))
impPath=url(("https://github.com/LTTGH/drkm4mc/raw/master/ModelImp1.RDS"))
```

Train the doubly-robust kernel machines
```{r}
options(warn=-1)
trainRes=trainingMisCov(misCovDat,px=4,kerType="RBF",kerMethod="DR",lambda=1,sigma=0.1,psPath,impPath,B=2,testPurpose="testing")

```

Using the doubly-robust kernel machines to predict Y in the testing dataset
and show the classification error
```{r}
dataTst=as.data.frame(cbind(testData$x1,testData$x2,testData$x3,
                        testData$x6,testData$x4,testData$x5,rep(1,110),testData$Y))
1-mean(testingMisCov(trainRes,dataTst)==testData$Y)
```
