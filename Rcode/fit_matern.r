###Fit Matern covariance using the Wendland function
rm(list=ls())

#system.time(source("Rcode/fit_matern.r",echo=TRUE))

#set the path to the main R code directory 
source("main.r",echo=TRUE)
path01 <- path0

source(file.path(path01,"Rcode","lik.r"),echo=TRUE)
source(file.path(path01,"Rcode","metropolis.r"),echo=TRUE)
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)

sparsity <- 0.9 # not useful
source(file.path(path01,"Rcode","config.r"),echo=TRUE)

outdir <- file.path(path01,"data_out","wendland")

library(fields)
#library(compiler)
library(sparseinv)
#library(LaplacesDemon)#

j <- 1
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")

#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
N <- 2000; r <- 6; ndim <- nrow(XX); p <- 24; m <- q <- 18; ns <- 600  #no of output 
Xsen0 <- array(XX,c(ndim,N,p))
YY <- data.matrix(fread(file.path(path01,"data_in","Ysen0.txt"),header=TRUE))# ndim times (N*m*r) matrix
Ysen0 <- array(YY,c(ndim,N,m,r))
rr <- data.matrix(read.table(file.path(path01,"data_in","sce_range.txt"),sep=",",header=FALSE)) #scenario ranges
rm(XX, YY)

s <- sample(1:N,ns)
	dat0 <- inp0 <- dat00 <- list()
	for(num in 1:ndim){
		dat0[[num]] <- (Ysen0[num,s,,j])##
		dat00[[num]] <- (Ysen0[num,-s,,j])#
	}
	inp1 <- list()
	inp0 <- ff(Xsen0[,s,])
	for(i in 1:ndim){
		 inp1[[i]]=inp0[i,,]
	}
	#
	inptest <- list()
	inp00 <- ff(Xsen0[,-s,])
	for(i in 1:ndim){
		inptest[[i]] <- inp00[i,,]
	}

	##########
	Xtrain <- abind(inp1,along=1)
	train <- abind(dat0,along=1)
	sc1 <- apply(train,2,sd)
	sc2 <- apply(train,2,mean)
	Ytrain <<- scale(train,center=sc2,scale=sc1)

##########
cor <- "Matern"
phi <- runif(p)
dd <- rdist(crossprod(t(Xtrain),diag(1/phi)))
phi.max <- max(dd)
tune <- .15*rep(2.38/p^2,p)
#options(spam.nearestdistnnz=c(4500000,400))
theta0 <- 0.05*phi.max ###taper ranges
options(spam.nearestdistnnz=c(4500000,400))

likelihood <- matern_lik(phi)
posterior <- likelihood$post
Z_rt <- Ytrain
X_rt <- Xtrain
rm(dd); gc()
###MCMC fitting
	result <- GP(NIter,Ytrain,Xtrain,Phi)
	setwd(outdir)
	save(result,file="result1")
	save(Xtrain,file="Xtrain1")
	save(Ytrain,file="Ytrain1")
	#crossvalidation
	
	############
pred <- function(xx0,mhat,omegahat,shat,phi){
	n <- nrow(Xtrain);n0=nrow(xx0)
	XX <- rbind(Xtrain,xx0)
	R0 <- stationary.taper.cov(XX,V=diag(phi), Covariance="Matern",smoothness=5/2,Taper.args= list(k=2,aRange=theta0,dimension=2))
	R <- as.dgCMatrix.spam(R0[1:n,1:n]) + Diagonal(n,nug)
	rr <- as.dgCMatrix.spam(R0[-c(1:n),-c(1:n)]) + Diagonal(n0,nug) #as.dgCMatrix.spam(rr) #+ Diagonal(n0,nug)
	r0 <- as.dgCMatrix.spam(R0[1:n,-c(1:n)])   #as.dgCMatrix.spam(r0) 
	F1 <- cholsolve(R,Ytrain-Matrix::crossprod(t(Xtrain),mhat),perm=TRUE)
	F2 <- cholsolve(R,Xtrain,perm=TRUE)
	mm <- as.matrix(Matrix::crossprod(t(xx0),mhat) + Matrix::crossprod(r0,F1))
	tq <- Matrix::t(xx0 - crossprod(r0,F2))
	vv <- as.matrix(rr - Matrix::crossprod(r0, cholsolve(R,r0,perm=TRUE))+
	Matrix::crossprod(tq, Matrix::crossprod(omegahat,tq)))
	V0 <- diag(vv)%*%t(diag(shat))/nus
	mu <- mm + sqrt(V0)* rt(n=1,df=nus) #matrix(rt(n=n0*q,df=nus),n0,q)
	return(list(mm=mm,V0=V0,mu=mu))
	#return(list(mm=mm,V0=V0))
}

	###########ESTIMATES 
	print(paste("rate=",length(unique(result$PHI[,1]))/N0))
	cf <- apply(result$PHI[-c(1:bb),],2,mean)
	mhat <- apply(result$THETA[-c(1:bb),,],c(2,3),mean)
	omegahat <- apply(result$OMEGA[-c(1:bb),,],c(2,3),mean)
	shat <- apply(result$SIGMANU[-c(1:bb),,],c(2,3),mean)
	VPP <- PP <- yy <- inp1 <- dat1 <- list()
	Ytest <- list()
	for(i in 1:ndim){
	Ytest[[i]] <- scale(dat00[[i]],center=sc2,scale=sc1)
    }
	for(i in 1:ndim){
		#print(i)
		Xtest <- inptest[[i]]
		pp <- pred(Xtest,mhat,omegahat,shat,phi)
		PP[[i]] <- pp[[1]]
		VPP[[i]] <- pp[[2]]
		yy[[i]] <- Ytest[[i]]
	}
	Yhat <- abind(PP,along=1); Ytest=abind(yy,along=1);VAR=abind(VPP,along=1)
	f2 <- function(k){RMSE(Ytest[,k],Yhat[,k])}
	f3 <- function(k){1-(sum((Ytest[,k] - Yhat[,k])^2)/sum((Ytest[,k] - mean(Ytest[,k]))^2))}
	res <- sapply(1:q,f2)
	res2 <- sapply(1:q,f3)
	ress <- print(signif(cbind(P=res2,rho=res),2))
	write.table(ress,file="model_summary.txt",row.names=FALSE)
	VES0 <- Yhat
	VES00 <- VAR
	save(VES0,file="VES0")
	save(VES00,file="VES00")
	save(Ytest,file="Ytest")
	setwd(path01)
	###END