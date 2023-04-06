##sensitivity function
sensitivity_func  <- function(hh,myf,models="Gmodel",CC,fun,q,p0,...){
p0 <- p0
dat <- myf(n0)
F=FJ=FJI=VF=VFJ=VFJI=list()
if(models=="Gmodel"){
	#print(hh)
	#dat <- Xset1(n0)
	#p0 <- ncol(dat[[1]])
	phi <- PHI[hh,];mhat=Mhat1[hh,,];omegahat=Omegahat1[hh,,];shat=Shat1[hh,,]
	temp <- lapply(1,fa,x=dat)[[1]]
	F <- temp[[1]]; #VF=temp[[2]]
	m2 <- lapply(1:p0,fa,x=dat[[2]])
	FJ <- abind(sapply(m2,`[[`, "mu", simplify = FALSE, USE.NAMES = FALSE),along=0)
	m3 <- lapply(1:p0,fa,x=dat[[3]])
	FJI <- abind(sapply(m3,`[[`, "mu", simplify = FALSE, USE.NAMES = FALSE),along=0)
}else{
if(identical(fun,fa)){
	#print(fun)
	temp <- lapply(1,fun,x=dat)[[1]]
	F <- matrix(temp,nrow=n0,ncol=1)
	m2 <- lapply(1:p0,fun,x=dat[[2]])
	FJ <- array(abind(m2,along=0),c(p0,n0,1))
	m3 <- lapply(1:p0,fun,x=dat[[3]])
	FJI <- array(abind(m3,along=0),c(p0,n0,1))
}else{
#print(fun)
	temp <- lapply(1,fun,x=dat)[[1]]
	F <- temp
	m2 <- lapply(1:p0,fun,x=dat[[2]])
	FJ <- abind(m2,along=0)
	m3 <- lapply(1:p0,fun,x=dat[[3]])
	FJI <- abind(m3,along=0)
}}
#####
	EY = colMeans(F); EYY = colMeans(F^2)
	EYJ = abind(lapply(1:p0,function(j)colMeans(FJ[j,,]*F)),along=0)
	EY_J = abind(lapply(1:p0,function(j)colMeans(FJI[j,,]*F)),along=0)
	COV = cov(F)
	###########Without output correlation or standard sobol
	RES1 = sapply(1:p0, function(j)(EYJ - matrix(EY*EY,p0,q,byrow=TRUE))[j,]/diag(COV)) 
	RES2 = 1 - sapply(1:p0, function(j)(EY_J - matrix(EY*EY,p0,q,byrow=TRUE))[j,]/diag(COV))   
	RES3 = sapply(1:p0,function(j)tr(cov(FJ[j,,],F)))/tr(COV)#overall_main
	RES4 = 1 - sapply(1:p0,function(j)tr(cov(FJI[j,,],F)))/tr(COV)#overall_total
	COV2a = sapply(1:p0,function(j)diag(cov(FJ[j,,],F)))
	COV2b = sapply(1:p0,function(j)diag(cov(FJI[j,,],F)))
	#vector projection
	if(q==1){
	P1 = (COV2a)%*%CC%*%diag(COV)/ (c(diag(COV)%*%CC%*%diag(COV)))
	P2 = 1 - (COV2b)%*%CC%*%diag(COV)/ (c(diag(COV)%*%CC%*%diag(COV)))
	}else{
	P1 = t(COV2a)%*%CC%*%diag(COV)/ (c(diag(COV)%*%CC%*%diag(COV)))
	P2 = 1 - t(COV2b)%*%CC%*%diag(COV)/ (c(diag(COV)%*%CC%*%diag(COV)))
	}
	RES1 = ifelse(RES1<0,0,RES1); RES2 = ifelse(RES2<0,0,RES2);RES3 = ifelse(RES3<0,0,RES3);RES4 = ifelse(RES4<0,0,RES4)
	P1 = ifelse(P1<0,0,P1);P2 = ifelse(P2<0,0,P2)
	list(RES1=t(RES1),RES2=t(RES2),RES3=RES3,RES4=RES4,P1=P1,P2=P2)
}
##
parallel.func <- function(hh,dat,models="Gmodel",CC){
p0 <- ncol(dat[[1]])
phi=PHI[hh,];mhat=Mhat1[hh,,];omegahat=Omegahat1[hh,,];shat=Shat1[hh,,]
cf <- phi
print(hh)
F=FJ=FJI=VF=VFJ=VFJI=list()
if(models=="Gmodel"){
	temp <- lapply(1,fa,x=dat)[[1]]
	F <- temp[[1]]; #VF=temp[[2]]
	m2 <- lapply(1:p0,fa,x=dat[[2]]) #foreach(j=1:p0,.packages=pack,.verbose=TRUE) %dopar% {fa(j,x=dat[[2]])}
	FJ <- abind(sapply(m2,`[[`, "mu", simplify = FALSE, USE.NAMES = FALSE),along=0)
	m3 <- lapply(1:p0,fa,x=dat[[3]]) #foreach(j=1:p0,.packages=pack,.verbose=TRUE) %dopar% {fa(j,x=dat[[3]])}
	FJI <- abind(sapply(m3,`[[`, "mu", simplify = FALSE, USE.NAMES = FALSE),along=0)
	}else{
	p0 <- ncol(dat[[1]])
	temp <- lapply(1,fb,x=dat)[[1]]
	F <- temp
	m2 <- lapply(1:p0,fb,x=dat[[2]])
	FJ <- abind(m2,along=0)
	#
	m3 <- lapply(1:p0,fb,x=dat[[3]])
	FJI <- abind(m3,along=0)
	q <- dim(FJI)[3]
}
#####
	EY = colMeans(F); EYY = colMeans(F^2)
	EYJ = abind(lapply(1:p0,function(j)colMeans(FJ[j,,]*F)),along=0)
	EY_J = abind(lapply(1:p0,function(j)colMeans(FJI[j,,]*F)),along=0)
	COV = cov(F)
	###########Without output correlation
	RES1 = sapply(1:p0, function(j)(EYJ - matrix(EY*EY,p0,q,byrow=TRUE))[j,]/diag(COV)) ##== diag(cov(FJ[j,,],F)/COV) ##m
	RES2 = 1 - sapply(1:p0, function(j)(EY_J - matrix(EY*EY,p0,q,byrow=TRUE))[j,]/diag(COV))   #or#== 1-diag(cov(FJI[j,,],F)/COV)##t
	RES3 = sapply(1:p0,function(j)tr(cov(FJ[j,,],F)))/tr(COV)#overall_main
	RES4 = 1 - sapply(1:p0,function(j)tr(cov(FJI[j,,],F)))/tr(COV)#overall_total
	COV2a = sapply(1:p0,function(j)diag(cov(FJ[j,,],F)))
	COV2b = sapply(1:p0,function(j)diag(cov(FJI[j,,],F)))
	P1 = t(COV2a)%*%CC%*%diag(COV)/ c(diag(COV)%*%CC%*%diag(COV))
	P2 = 1 - t(COV2b)%*%CC%*%diag(COV)/ c(diag(COV)%*%CC%*%diag(COV))
	RES1 = ifelse(RES1<0,0,RES1); RES2 = ifelse(RES2<0,0,RES2);RES3 = ifelse(RES3<0,0,RES3);RES4 = ifelse(RES4<0,0,RES4)
	P1 = ifelse(P1<0,0,P1);P2 = ifelse(P2<0,0,P2)
	list(RES1=t(RES1),RES2=t(RES2),RES3=RES3,RES4=RES4,P1=P1,P2=P2)
}

Xset1 <- function(n0){
XXJ <- XX_J  <- list()
	X0 <- matrix(runif(n0*p0,-1,1), nrow = n0)#
	#X_j <- matrix(runif(n0*p0,-1,1), nrow = n0)
	########Xj
	for(j in 1:p0){
		XXJ[[j]] <- matrix(runif(n0*p0,-1,1), nrow = n0)
		XXJ[[j]][,j] <- X0[,j]
	}
	#
	for(j in 1:p0){
		XX_J[[j]] = X0
		XX_J[[j]][,j] = runif(n0,-1,1)
	}
	XX=list(X0,XXJ,XX_J)
}
#

fa <- function(j,x){
	colnames(x[[j]])=nam
	XX <- model.matrix(quad,as.data.frame(x[[j]]))[,-1]
	pred(XX,cf,mhat,omegahat,shat)
	}