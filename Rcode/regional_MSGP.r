########## The main code is for fitting the regional and EU_wide MSGP emulators
#system.time(source("../JRSC_paper_code/Rcode/regional_MSGP.r",echo=TRUE))
rm(list=ls())

#set the path to the main R code directory eg ../JRSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","lik.r"),echo=TRUE)
source(file.path(path01,"Rcode","mcmc_func.r"),echo=TRUE)
source(file.path(path01,"Rcode","metropolis.r"),echo=TRUE)
source(file.path(path01,"Rcode","pred.r"),echo=TRUE)

sparsity <- 0.9 #fix
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")
sp_levels <- c(.99,.95,.90,.80)

#source MCMC initilization parameters
source(file.path(path01,"Rcode","config.r"),echo=TRUE)

#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
N <- 2000; r <- 6; ndim <- nrow(XX); m <- q <- 18; p <- 24; ns <- 600
Xsen0 <- array(XX,c(ndim,N,p))
YY <- data.matrix(fread(file.path(path01,"data_in","Ysen0.txt"),header=TRUE))# ndim times (N*m*r) matrix
Ysen0 <- array(YY,c(ndim,N,m,r))
rr <- data.matrix(read.table(file.path(path01,"data_in","sce_range.txt"),sep=",",header=FALSE)) #scenario ranges
rm(XX, YY)

#Fit 5 regional models, fix sparsity to 0.90
out <- out2 <- out3 <- list()

for(t0 in 2:length(listof_region)){
s <- sample(1:N,ns)
	dat0 <- inp0 <- dat00 <- list()
	for(num in 1:ndim){
		dat0[[num]] <- (Ysen0[num,s,,t0])##
		dat00[[num]] <- (Ysen0[num,-s,,t0])#
	}
	inp1 <- list()
	inp0 <- ff(Xsen0[,s,])
	for(i in 1:ndim){
		 inp1[[i]] <- inp0[i,,]
	}
	#
	inptest <- list()
	inp00 <- ff(Xsen0[,-s,])
	for(i in 1:ndim){
		inptest[[i]] <- inp00[i,,]
	}

	###
	Xtrain <- abind(inp1,along=1)
	train <- abind(dat0,along=1)
	sc1 <- apply(train,2,sd)
	sc2 <- apply(train,2,mean)
	Ytrain <- scale(train,center=sc2,scale=sc1)
	likelihood <- lik(cf)
posterior <- likelihood$post
#Fitting the models
	print(t0)
	region <- listof_region[t0]
	#create output dir
	outdir <- file.path(path01,"data_out",gsub("t0",t0,listof_region[t0]))
	 dir_create(outdir)
	out[[t0]] <- mcmc_func(sparsity)
}

#perform regional sensitivity
system.time(source(file.path(path01,"Rcode","regional_sensitivity.r"),echo=TRUE))
