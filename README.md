# Modified-Weight-Function
# Gene Expression Omnibus dataset (accession number GSE 65622)
# After downloading the dataset in CSV format, our first job is to make a transposition. Secondly, the irrelevant columns (like patient-id) are removed from the data and finally, we rename the file as "wdata3.csv" and save in a folder.

# Required packages to run the code
library(base)
library(bdsmatrix)
library(Matrix)
library(dplyr)
library(survival)
library(glmnet)

# Step-1: Data reshaping, sorting and cleaning
setwd("C:/Users/Souvik/Desktop/Clinical Research/Weight-paper") #store the dataset
wdata3 = read.csv("wdata3.csv") #read the data in R
wdat = data.frame(wdata3[,14:520])

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
ntotal = 0 #total number of features which have less than 5% missing values
for(j in 1:ncol(wdat))
{
  if(count[j] <= ncol(wdat)*0.05)
    ntotal = ntotal + 1 
}
temp = c() #selection of features which have more than 5% missing values
for(j in 1:ncol(wdat))
{
  if(count[j] > ncol(wdat)*0.05)
  {
    wdat[,j] = NA
    temp[j] = j
  }
}
u1 = c(na.omit(temp))
wdt = wdat[,-u1] #dataset of the selected features which have less than 5% missing values
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
D = data.frame(wdata3$Sample_title, wdata3$visit_time, wdata3$OS, wdata3$Relapse, 
               wdata3$RFS, wdata3$Metastasis, wdata3$DMFS, d2) # Datasets of the selected features where missing values are replaced by the mean of the variables
colnames(D) = c("Sample_title", "visit_time", "OS", "Relapse", "RFS", 
                "Metastasis", "DMFS", c(names(wdt)))

# Step-2: Univariate Cox-PH Fitting and use filter method
genedata = data.frame(RFS = D$RFS, Relapse = D$Relapse, D[,8:334])
pval = c() #calculation of p-value from cox regression
for(i in 1:(dim(genedata)[2]-2))
{
  pval[i] = round(summary(coxph(Surv(RFS, Relapse)~genedata[,i+2], data = genedata))$waldtest[3],4)
}
D.pval = data.frame(names(D[,3:329]), pval)
Dpv = D.pval[D.pval[,2]<=0.05,] #list of significant variables and their p-values

DMain = data.frame(RFS = D$RFS, Relapse = D$Relapse, D[Dpv[,1]]) #dataset containing significant variables

# Step-3: Function for selecting the features on the basis of modified weights 
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
  
  Dframe = data.frame(Name = names(data[,m1:m2]), Pweight = w, 
                      CoefficientofVariation = CV, Cweight = w.modified)
  D.weight = Dframe[order(Dframe[,4], na.last = NA, decreasing = TRUE),]
  
  selectfeature = D.weight
  SelectFeature = selectfeature[selectfeature$ModifiedWeight > threshold,] 
  #feature selection for final analysis
  return(SelectFeature)
}
feat.selec = weightfun(studyvar = "RFS", status = "Relapse", data = DMain,
                        m1 = 3, m2 = 102, p.train = 0.75, N = 150, threshold = 0.90)
                        
# Cox Proportional Hazards Model Fitting
CoxPH.modified = coxph(Surv(RFS, Relapse)~DMain$FGF.5+DMain$GASP.2...WFIKKN+DMain$MSP.alpha.Chain
                       +DMain$TRAIL.R3...TNFRSF10C+DMain$L.Selectin..CD62L.+DMain$ICAM.1
                       +DMain$PDGF.BB+DMain$Eotaxin.3...CCL26+DMain$PDGF.R.beta+DMain$IFN.gamma
                       +DMain$Endoglin...CD105+DMain$FGF.13.1B, data = DMain) #for C-weight
CoxPH.old = coxph(Surv(RFS, Relapse)~DMain$FGF.5+DMain$GASP.2...WFIKKN+DMain$MSP.alpha.Chain, data = DMain) #for P-weight
