### Modified-Weight-Function (R-Code)
### Step-1: Data reshaping, sorting and cleaning

### Gene Expression Omnibus Dataset (accession number GSE 65622)
### After downloading the dataset in CSV format, our first job is to make a transposition. Secondly, the irrelevant columns (like patient-id) are removed from the data and finally, we rename the file as "wdata3.csv" and save in a folder.

setwd("C:/Users/Souvik/Desktop/Clinical Research/Weight-paper") #store the data
wdata3 = read.csv("wdata3.csv") #read the data in R
wdat = data.frame(wdata3[,14:520])
dim(wdat) #dimension of the dataset

count = c() #count the number of null values within each variable
for(j in 1:ncol(wdat))
{
  count[j] = 0
  for(i in 1:nrow(wdat))
  {
    if(wdat[i,j] == "null")
    {
      count[j] = count[j] + 1
    }
  }
}

ntotal = 0 #total number of features who have less than 5% missing values
for(j in 1:ncol(wdat))
{
  if(count[j] <= ncol(wdat)*0.05)
    ntotal = ntotal + 1
}

temp = c() #select those features who have more than 5% missing values
for(j in 1:ncol(wdat))
{
  if(count[j] > ncol(wdat)*0.05)
  {
    wdat[,j] = NA
    temp[j] = j
  }
}

u1 = c(na.omit(temp))
wdt = wdat[,-u1] #dataset of those selected features who have less than 5% missing values
d0 = matrix(0, nrow(wdt), (ncol(wdt)))
d1 = matrix(0, nrow(wdt), (ncol(wdt)))
d2 = matrix(0, nrow(wdt), (ncol(wdt)))
for(j in 1:(ncol(wdt)))
{
  for(i in 1:nrow(wdt))
  {
    if(wdt[i,j] == "null"){
      d0[i,j] = 0
    }else{
      d0[i,j] = wdt[i,j]
    }
  }
  d1[,j] = as.numeric(d0[,j])
  s = d1[,j]
  g1 = c()
  for(i in 1:nrow(wdt))
  {
    if(s[i] == 0)
    {
      g1[i] = i
    }
  }
  g1 = na.omit(g1)
  d2[,j] = replace(d1[,j], list = g1, values = mean(s[s>0], na.rm = T))
}

### Datasets of the selected features where missing values are replaced by the mean of the variables
D = data.frame(wdata3$Sample_title, wdata3$visit_time, wdata3$OS, wdata3$Relapse, 
               wdata3$RFS, wdata3$Metastasis, wdata3$DMFS, d2) 
colnames(D) = c("Sample_title", "visit_time", "OS", "Relapse", "RFS", 
                "Metastasis", "DMFS", c(names(wdt)))
dim(D) #dimension of the dataset for final analysis

### Step-2: Function for selecting the features on the basis of modified weights 

### Required packages
library(base)
library(bdsmatrix)
library(Matrix)
library(dplyr)
library(survival)
library(glmnet)

#studyvar = Variable name of the study variable
#status = Status of the study variable (i.e. 1 = event occurred, 0 = no event occurred)
#data = Name of the high-dimensional dataset containing no missing values
#m1 = Initial column number from where variables of high dimensional data will be selected
#m2 = Ending column number till where variable of high dimensional data will be selected
#p.train = Proportion of traing dataset
#N = Number of iterations
#threshold = A threshold weight for which the features having weight more than the threshold will be included to the final model

weightfun = function(studyvar, status, data, m1, m2, p.train, N, threshold)
{
  StudyVar = pull(data, studyvar) #study variable
  Status = pull(data, status) #status of the study variable
  k = m2-m1+1 #total number of features who have less than 5% missing values
  
  Dt = matrix(0, ceiling(nrow(data)*p.train), N)
  for(j in 1:N)
  {
    Dt[,j] = sort(sample(nrow(data), 
                         size = ceiling(nrow(data)*p.train), 
                         replace = FALSE))
  }
  
  optimal_lambda.lasso.train = c()
  g = matrix(0, k, N)
  for(j in 1:N)
  {
    dt = Dt[,j]
    D.train = data[dt,] #training dataset
    trainstudy = StudyVar[dt]
    trainstatus = Status[dt]
    
    lassofit.train = glmnet(as.matrix(D.train[,m1:m2]), 
                            Surv(trainstudy, trainstatus), 
                            family="cox", alpha = 1) #lasso regression for cox model
    
    optimal_lambda.lasso.train[j] = cv.glmnet(as.matrix(D.train[,m1:m2]), 
                                              Surv(trainstudy, trainstatus), 
                                              family="cox", alpha = 1)$lambda.min 
    #cross-validation for glmnet to find the optimal lambda
    
    coeff.lassotrain = coef(lassofit.train, s = optimal_lambda.lasso.train[j]) 
    #coefficient vectors of the features
    
    for(i in 1:length(coeff.lassotrain))
    {
      if(coeff.lassotrain[i] != 0){
        g[i,j] = i #selection of i-th feature at j-th iteration
      }else{
        g[i,j] = NA
      }
    }
  }
  
  n1 = c() #it calculates how number of times each features are selected
  w = c() #weights
  CV = c() #coefficient of variation(CV)
  for(i in 1:k)
  {
    if(length(na.omit(g[i,])) == 0){
      n1[i] = NA
    }else{
      n1[i] = length(na.omit(g[i,]))
    }
    w[i] = (n1[i]/N)
    CV[i] = sd(data[,m1+i-1])/mean(data[,m1+i-1])
  }
  range.CV = max(CV) - min(CV) #range of the CV values
  
  w.modified = c() #modified weights
  for(i in 1:k)
  {
    w.modified[i] = exp(-0.5 * (1-w[i]) * (CV[i]/range.CV))
  }
  
  Dframe = data.frame(Name = names(data[,m1:m2]), OldWeight = w, CoefficientofVariation = CV, ModifiedWeight = w.modified)
  D.weight = Dframe[order(Dframe[,4], na.last = NA, decreasing = TRUE),]
  
  selectfeature = D.weight
  SelectFeature = selectfeature[selectfeature$ModifiedWeight > threshold,] 
  #feature selection for final analysis
  
  return(SelectFeature)
}

weightfun(studyvar = "RFS", status = "Relapse", data = D, m1 = 8, m2 = 334, p.train = 0.75, N = 150, threshold = 0.90)
