#this function fit the MSGP model and perform the crossvalidation
mcmc_func <- function(sparsity){
  ###MCMC fitting
	result <- GP(NIter,Ytrain,Xtrain,Phi)
	setwd(outdir)
	save(result,file="result1")
	save(Xtrain,file="Xtrain1")
	save(Ytrain,file="Ytrain1")
	#crossvalidation
	########scale output
	
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
		pp <- pred(Xtest,cf,mhat,omegahat,shat)
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
     save(Xtest,file="Xtest")
	setwd(path01)
	return(ress)
	}