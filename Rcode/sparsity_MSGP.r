########## This code is for fitting the different levels of sparsity to the EU_wide MSGP emulator
#system.time(source("../JRSC_paper_code/Rcode/sparsity_MSGP.r",echo=TRUE))
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
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")
sp_levels <- c(.99,.95,.90,.80)
region <- "EU_wide"

#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
N <- 2000; r <- 6; ndim <- nrow(XX); p <- 24; m <- q <- 18; ns <- 600
Xsen0 <- array(XX,c(ndim,N,p))
YY <- data.matrix(fread(file.path(path01,"data_in","Ysen0.txt"),header=TRUE))# ndim times (N*m*r) matrix
Ysen0 <- array(YY,c(ndim,N,m,r))
rr <- data.matrix(read.table(file.path(path01,"data_in","sce_range.txt"),sep=",",header=FALSE)) #scenario ranges
rm(XX, YY)
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")

out <- out2 <- out3 <- list()
#fit different sparsity levels, fix the region to "EU_wide"

t0 <- match(region,listof_region)
s <- sample(1:N,ns)
	dat0 <- inp0 <- dat00 <- list()
	for(num in 1:ndim){
		dat0[[num]] <- (Ysen0[num,s,,t0])##
		dat00[[num]] <- (Ysen0[num,-s,,t0])#
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

	###
	Xtrain <- abind(inp1,along=1)
	train <- abind(dat0,along=1)
	sc1 <- apply(train,2,sd)
	sc2 <- apply(train,2,mean)
	Ytrain <<- scale(train,center=sc2,scale=sc1)
	
for(w0 in 1:length(sp_levels)){
sparsity <- sp_levels[w0]
source(file.path(path01,"Rcode","config.r"),echo=TRUE)
likelihood <- lik(cf)
posterior <- likelihood$post
#Fitting the models
print(w0)
	   #create output dir
   outdir <- file.path(path01,"data_out",paste("Sparsity",sparsity*100,sep="_"))
   dir_create(outdir)
out2[[w0]] <- mcmc_func(sparsity)
}

#####END