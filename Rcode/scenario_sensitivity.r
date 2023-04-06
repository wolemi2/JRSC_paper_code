###Scenario-based sensitivity indices
#system.time(source(file.path("Rcode","scenario_sensitivity.r"),echo=TRUE))
rm(list=ls())

#set the path to the main R code directory eg ../JRSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","pred.r"),echo=TRUE)
source(file.path(path01,"Rcode","sensitivity_func.r"),echo=TRUE)

outdir0 <- file.path(path01,"data_out","EU_wide")
outdir <- file.path(path01,"data_out","scenario_based_sensitivity")
dir_create(outdir)

N0 <- 1000#nrow(PHI)
p <- 24; m <- q <- 18; n <- 9000
bb <- 0.75*N0
delta0 <- -q+1
nus <- df <- delta0+n
cor <- "truncpow"; nug <- .02; nu <- 2; alpha <- 3/2
#
p0 <- p 
ll <- sample(c(bb+1):N0,50)
n0 <- 1000
sc <- 4
####NEW_SCENARIO_DATA
#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
N <- 2000; r <- 6; ndim <- nrow(XX)
NT <- N*ndim #30000
Xsen0 <- array(XX,c(ndim,N,p))
YY <- data.matrix(fread(file.path(path01,"data_in","Ysen0.txt"),header=TRUE))# ndim times (N*m*r) matrix
Ysen0 <- array(YY,c(ndim,N,m,r))
rr <- data.matrix(read.table(file.path(path01,"data_in","sce_range.txt"),sep=",",header=FALSE)) #scenario ranges
sce_data <- data.matrix(fread(file.path(path01,"data_in","scenario.txt"),header=TRUE))# ndim times (N*m*sc) matrix
sce_data <- array(sce_data,c(ndim,N,m,sc))
rm(XX, YY)
gc()

#######
y0 <- sce_dat <- list()
for(i in 1:ndim){
	y0[[i]] <- sce_data[i,,,]
	}
Ys <- abind(y0,along=1)

for(i in 1:sc){
	sce_dat[[i]] <- scale(Ys[,,i],center=sc2,scale=sc1)
}

#load EU-wide fitted parameters
load(file.path(outdir0,"result1"))
	load(file.path(outdir0,"Ytrain1"))
	load(file.path(outdir0,"Xtrain1"))

	cf <- apply(result$PHI[-c(1:bb),],2,mean)
	mhat <- apply(result$THETA[-c(1:bb),,],c(2,3),mean)
	omegahat <- apply(result$OMEGA[-c(1:bb),,],c(2,3),mean)
	shat <- apply(result$SIGMANU[-c(1:bb),,],c(2,3),mean)
	#
	PHI <- result$PHI
	Omegahat1 <- result$OMEGA
	Mhat1 <- result$THETA
	Shat1 <- result$SIGMANU
#####
sce_mod <- list()
for(j in 1:length(sce_dat)){
print(j)
	clusterExport(cl, extra_func)
	system.time(mod <- foreach(hh =ll,.packages=pack,.verbose=TRUE) %dopar% {sensitivity_func(hh,myf=Xset1,models="Gmodel",CC=cor(sce_dat[[j]]),fun=fa,q=q,p0=p0)})
	#system.time(mod <- lapply(ll,sensitivity_func,myf=Xset1,models="Gmodel",CC=cor(sce_dat[[j]]),fun=fa,q=q,p0=p0))
	sce_mod[[j]] <- mod
	save(sce_mod,file=file.path(outdir,"sce_mod"))
}
save(sce_mod,file=file.path(outdir,"sce_mod"))