######This code fits Standard GP using DiceKriging library and multivariate regression models to a sample of EU_wide data 
rm(list=ls())

#set the path to the main R code directory eg ../JRSSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
library(DiceKriging)
library(xtable)
library(abind)
#library(mvlm)
  
  dir_create <- function(path0){{if(!dir.exists(path0)){
    dir.create(path0, recursive = TRUE)
  }}}
  RMSE <- function(z, zhat) {
        z <- as.matrix(z)
        zhat <- as.matrix(zhat)
        x <- c(z - zhat)
        u <- x[!is.na(x)]
        round(sqrt(sum(u^2)/length(u)), 4)
    }

  quad <- formula(~Emission + model + SESs + Population + GreenRed + RumLivesDem +
    NonRumiLivesDem + StructChange + TechFactor + TechChange +
    YieldFactor + IrrigEffFactor + GDP + CostsFactor + ImportFactor +
    BioEneCropDemand + ArableConservLandSc + CropInputsFactor +
    Coastal + Fluvial + HabitatRecreation + PA_Forestry + PA_Agricultural +
    FloodProtection)
	nam2=c("Food Per Capita","People Flooded 1/100year","Water Exploitation Index","Landuse Intensity Index","Landuse Diversity",
"Biodiversity Vulnerability Index","Timber Production", "Area of Artificial Surfaces"
,"Food Production","Potential Carbon Stock","Irrigation Usage","Intensive Arable","Intensive Grassland","Extensive Grassland",
"Very Extensive Grassland","Unmanaged Land","Managed Forest","Unmanaged Forest")

   
outdir <- file.path(path01,"data_out","multivariate")
outdir0 <- file.path(path01,"data_out","EU_wide")
dir_create(outdir)

j <- 1
load(file.path(outdir0,"Ytrain1"))
load(file.path(outdir0,"Xtrain1"))
s <- sample(1:nrow(Xtrain),4500)

epsilon <- 1e-3
system.time(krig_model <- km(formula=quad, design=Xtrain[s,], response=Ytrain[s,j],nugget=epsilon))
saveRDS(krig_model,file=file.path(outdir,"krig_model"))

##############Fit standard multivariate regression
load(file.path(outdir0,"Ytest"))
load(file.path(outdir0,"inptest"))

result1 <- lm(Ytrain~.,data.frame(Xtrain))
save(result,file=file.path(outdir,"result1"))

ytest <- Ytest
Xtest <- as.data.frame(abind(inptest,along=1))
colnames(Xtest) <- colnames(Xtrain)
Yhat <- predict(result1,newdata=Xtest)

f2 <- function(k){RMSE(ytest[,k],Yhat[,k])}
f3 <- function(k){1-(sum((ytest[,k] - Yhat[,k])^2)/sum((ytest[,k] - mean(ytest[,k]))^2))}
linear_reg0 <- print(signif(sapply(1:ncol(Ytrain),f2),2))
linear_reg <- print(signif(sapply(1:ncol(Ytrain),f3),2))
cbind(linear_reg,linear_reg0)

matern0 <- c(0.46,0.17,0.59,0.33,0.38,0.58,0.48,0.69,0.59,0.66,0.70,0.59,0.56,0.70,0.83,0.69,0.63,0.51)# RMSE
matern <- c(0.54,0.95,0.40,0.82,0.76,0.55,0.69,0.47,0.73,0.61,0.57,0.68,0.75,0.54,0.39,0.55,0.65,0.77) #Prop
rel_performance0 <- round(((linear_reg0 - matern0)/linear_reg0)*100,1)
rel_performance <- round(((matern - linear_reg)/linear_reg)*100,1)
res2 <- data.frame(nam2,rel_performance0)

xtable(res2)
round(((mean(linear_reg0) - mean(matern0))/mean(linear_reg0))*100,1)