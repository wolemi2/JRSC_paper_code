######### The main code for fitting the EU_wide MSGP emulator
#system.time(source("../JRSC_paper_code/Rcode/EU_wide_MSGP.r",echo=TRUE))
rm(list=ls())
set.seed(15)

#set the path to the main R code directory eg ../JRSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","lik.r"),echo=TRUE)
source(file.path(path01,"Rcode","mcmc_func.r"),echo=TRUE)
source(file.path(path01,"Rcode","metropolis.r"),echo=TRUE)
source(file.path(path01,"Rcode","pred.r"),echo=TRUE)

#source MCMC initiliaztion parameters
sparsity <- 0.9
source(file.path(path01,"Rcode","config.r"),echo=TRUE)
#dir_create(outdir)

j <- 1
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")
sp_levels <- c(.99,.95,.90,.80)

#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
N <- 2000; r <- 6; ndim <- nrow(XX)
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
	likelihood <- lik(cf)
posterior <- likelihood$post

############Fitting the models
#Fit EU_wide models, fix sparsity to 0.90
out <- list()
t0 <- 1
	print(t0)
	region <- listof_region[t0]
	#create output dir
	outdir <- file.path(path01,"data_out",gsub("t0",t0,listof_region[t0]))
	 dir_create(outdir)
	out[[t0]] <- mcmc_func(sparsity=0.90)

###produce diagnostic plots
source(file.path(path01,"Rcode","diagnostic.r"),echo=TRUE)

#perform EU_wide sensitivity
#system.time(source(file.path(path01,"Rcode","EU_wide_sensitivity.r"),echo=TRUE))
