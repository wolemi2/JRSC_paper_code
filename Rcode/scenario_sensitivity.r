###Scenario-based sensitivity indices
#system.time(source(file.path("Rcode","scenario_sensitivity.r"),echo=TRUE))
rm(list=ls())

#set the path to the main R code directory eg ../JRSSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","pred.r"),echo=TRUE)
source(file.path(path01,"Rcode","sensitivity_func.r"),echo=TRUE)

outdir0 <- file.path(path01,"data_out","scenario_based_sensitivity")
dir_create(outdir0)

listof_scenario <- c("SSP1","SSP3","SSP4","SSP5")
sce_mod <- list()
p <- 24; m <- q <- 18; n <- 9000
delta0 <- -q+1
	nus <- df <- delta0+n
	cor <- "truncpow"; nug <- .02; nu <- 2; alpha <- 3/2
	p0 <- p; n0 <- 1000; n00 <- 25

for(t0 in 1:length(listof_scenario)){
	print(t0)
	region <- listof_scenario[t0]
	outdir <- file.path(path01,"data_out",region)
	load(file.path(outdir,"result1"))
	load(file.path(outdir,"Ytrain1"))
	load(file.path(outdir,"Xtrain1"))

	N0 <- nrow(result$PHI)
	bb <- 0.75*N0
	ll <- sample(c(bb+1):N0,n00)
	
	cf <- apply(result$PHI[-c(1:bb),],2,mean)
	mhat <- apply(result$THETA[-c(1:bb),,],c(2,3),mean)
	omegahat <- apply(result$OMEGA[-c(1:bb),,],c(2,3),mean)
	shat <- apply(result$SIGMANU[-c(1:bb),,],c(2,3),mean)
	#
	PHI <- result$PHI
	Omegahat1 <- result$OMEGA
	Mhat1 <- result$THETA
	Shat1 <- result$SIGMANU
    clusterExport(cl, extra_func)
	#system.time(mod <- foreach(hh =ll,.packages=pack,.verbose=TRUE) %dopar% {sensitivity_func(hh,myf=Xset1,models="Gmodel",CC=cor(Ytrain),fun=fa,q=q,p0=p0)})
	system.time(mod <- lapply(ll,sensitivity_func,myf=Xset1,models="Gmodel",CC=cor(Ytrain),fun=fa,q=q,p0=p0))
	sce_mod[[t0]] <- mod
	save(mod,file=file.path(outdir,"mod"))
#}
save(sce_mod,file=file.path(outdir0,"sce_mod"))