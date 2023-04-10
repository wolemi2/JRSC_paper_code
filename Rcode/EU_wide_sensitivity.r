#######This R code performs sensitivity analysis for the EU_wide MSGP emulators

#set the path to the main R code directory eg ../JRSSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","pred.r"),echo=TRUE)
source(file.path(path01,"Rcode","sensitivity_func.r"),echo=TRUE)

listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")
sp_levels <- c("Sparsity_99","Sparsity_95","Sparsity_90","Sparsity_80")

outdir <- file.path(path01,"data_out","EU_wide")

#load fitted results
load(file.path(outdir,"result1"))
load(file.path(outdir,"Ytrain1"))
load(file.path(outdir,"Xtrain1"))

#MSGP estimated parameters
PHI <- result$PHI
Omegahat1 <- result$OMEGA
Mhat1 <- result$THETA
Shat1 <- result$SIGMANU
n <- nrow(Xtrain)
N0 <- nrow(PHI)
m <- q <- 18 
bb <- 0.75*N0
delta0 <- -q+1
nus <- df <- delta0+n
cor <- "truncpow"; nug <- .02; nu <- 2; alpha <- 3/2
cf <- apply(result$PHI[-c(1:bb),],2,mean)
mhat <- apply(result$THETA[-c(1:bb),,],c(2,3),mean)
omegahat <- apply(result$OMEGA[-c(1:bb),,],c(2,3),mean)
shat <- apply(result$SIGMANU[-c(1:bb),,],c(2,3),mean)

p0 <- p <- ncol(PHI)
ll <- sample(c(bb+1):N0,50)
n0 <- 1000 
#

clusterExport(cl, extra_func)
system.time(mod <- foreach(hh =ll,.packages=pack,.verbose=TRUE) %dopar% {sensitivity_func(hh,myf=Xset1,models="Gmodel",CC=cor(Ytrain),fun=fa,q=q,p0=p0)})
#system.time(mod <- lapply(ll,sensitivity_func,myf=Xset1,models="Gmodel",CC=cor(Ytrain),fun=fa,q=q,p0=p0))

save(mod,file=file.path(outdir,"mod"))
################################################Main effect computation
sim <- seq(-1,1,.1)
xx <- matrix(sim,nrow=length(sim),ncol=p)
effect <- function(hh){
	phi <- PHI[hh,];mhat=Mhat1[hh,,];omegahat=Omegahat1[hh,,];shat=Shat1[hh,,]
	Ej <- Vj <- list()
	xx = matrix(seq(-1,1,.1),nrow=21,ncol=p)
	p1 = pred(xx,cf,mhat,omegahat,shat)
	for(j in 1:p){
		xj = xx; xj[,-j] <- 0
		p2 = pred(xj,cf,mhat,omegahat,shat)
		EY = p1$mm; EY_uj = p2$mm
		ej = EY_uj - EY 
		vj = var(EY) + var(EY_uj) -2*cov(EY,EY_uj)
		Ej[[j]] = ej
		Vj[[j]] = vj
	}
	EJ <- abind(Ej,along=0)
	VJ <- abind(Vj, along=0)
	res <- list(EJ=EJ,VJ=VJ)
}

Effect <- lapply(ll,effect)
save(Effect,file=file.path(outdir,"Effect"))
##############END
