####This routine is for fitting the test cases and plotting Figures 2,3 and 4.
rm(list=ls())
set.seed(11)
library(sensitivity)
#set the path to the main R code directory eg ../JRSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

#source relevant R functions
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
source(file.path(path01,"Rcode","sensitivity_func.r"),echo=TRUE)

outdir <- file.path(path01,"data_out","test_cases")
dir_create(outdir)

#Test case example 1.
##########################################
# Test case: the non-monotonic Sobol g-function
# Example with a call to a numerical model
# First compute first-order indices with ranking
n0 <- 1000
X <- data.frame(matrix(runif(8 * n0), nrow = n0))
x1 <- sobolshap_knn(model = sobol.fun, X = X, U = 1, method = "rank")
print(x1)

# We can use the output sample generated for this estimation to compute total indices
# without additional calls to the model
x2 <- sobolshap_knn(model = sobol.fun, X = X, U = 0, method = "rank")
#x22 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
#tell(x2,x1$y)

saveRDS(x1,file=file.path(outdir,"x1"))
saveRDS(x2,file=file.path(outdir,"x2"))

#Sobol Analytical value
sobolf <- function(i){
    a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
	y1 <- 1
	vi <- (1/3)/(1+a[i])^2
    for (i in 1:length(a)) {
	#i <- (1/3)/(1+a[i])^2
	 y1 <- (y1*prod(1+(1/3)/(1+a[i])^2)) 
	     }
		 r=vi/(y1-1)
list(r,vi,(y1-1))	 
 }
 
ana_value1 <- sobolf(1:8)
#total sensitivity
sobolf2 <- function(i,vi){
k <- i
for(i in 1:k){
v2 <- vi[[2]][i]*prod(1+vi[[2]][-i])
}
v2/vi[[3]]
}
ana_value2 <- sapply(1:8,function(i)sobolf2(i,ana_value1))
saveRDS(ana_value1,file=file.path(outdir,"ana_value1"))
saveRDS(ana_value2,file=file.path(outdir,"ana_value2"))

Xset1 = function(n0){
XXJ = XX_J =list()
X0= matrix(runif(n0*p0,0,1), nrow = n0)#
#X_j = matrix(runif(n0*p0,0,1), nrow = n0)
########Xj
for(j in 1:p0){
XXJ[[j]] = matrix(runif(n0*p0,0,1), nrow = n0)
XXJ[[j]][,j] = X0[,j]
}
#
for(j in 1:p0){
XX_J[[j]] = X0
XX_J[[j]][,j] = runif(n0,0,1)
}
XX=list(X0,XXJ,XX_J)
}

##Example 2 / atantemp.fun
X1 <- data.frame(matrix(runif(2*n0,-7,7), nrow = n0))
X2 <- data.frame(matrix(runif(2*n0,-7,7), nrow = n0))
x22 <- sobolMultOut(model=atantemp.fun, q=100, X1, X2,MCmethod="soboljansen",ubiquitous=TRUE)
png(file.path(outdir,"f3a.png"), height =10, width=12,units="in",res=400)
plotMultOut(x22)
dev.off()
saveRDS(x22,file=file.path(outdir,"x22"))

Xset2 = function(n0){
XXJ = XX_J =list()
X0= matrix(runif(n0*p0,-7,7), nrow = n0)#
########Xj
for(j in 1:p0){
XXJ[[j]] = matrix(runif(n0*p0,-7,7), nrow = n0)
XXJ[[j]][,j] = X0[,j]
}
#
for(j in 1:p0){
XX_J[[j]] = X0
XX_J[[j]][,j] = runif(n0,-7,7)
}
XX=list(X0,XXJ,XX_J)
}

###########################
#Test case 3: #Example 3.2 in Xu et al 2019 paper
ll <- 1000 #sample(1:N0,100)
n0 <- 1000
#p0 <- 8

Xset3 <- function(n0){
	XXJ <- XX_J <- list()
	mm <- c(20000,12,0.04,9.82*10^(-4),1.34*10^7,3.35*10^8, 2*10^10,1*10^11)
	ss <- c(1400,0.12,0.0048,5.89*10^(-5),2.412*10^6,4.02*10^7,1.2*10^9,6*10^9)
	CC <- matrix(c(1,0.3122,0.3694,0.3112,1,0.0467,0.3694,0.0467,1),nrow=3)
	X0 <- sapply(1:length(mm),function(i)rnorm(n0,mm[i],ss[i]))
	p0 <- ncol(X0)
	########Xj
	for(j in 1:p0){
	XXJ[[j]] <- sapply(1:length(mm),function(i)rnorm(n0,mm[i],ss[i]))
	XXJ[[j]][,j] <- X0[,j]
	}
	#
	for(j in 1:p0){
	XX_J[[j]] = X0
	XX_J[[j]][,j] <- rnorm(n0,mm[j],ss[j]) # runif(n0,-1,1)
	}
XX <- list(X0,XXJ,XX_J)
}
ex_fun <- function(X){
	x_1 <- X[,1];x_2 <- X[,2];x_3 <- X[,3];x_4 <- X[,4]
	x_5 <- X[,5];x_6 <- X[,6];x_7 <- X[,7];x_8 <- X[,8]
	y1 <- 0.03 - (1.905*x_1*x_2^2/(x_3*x_7)) - (0.565*x_1*x_2^2/(x_4*x_8))
	y2 <- x_5*x_3 - 1.185*x_1*x_2  
	y3 <- x_6*x_4 - 0.75*x_1*x_2
	Y <- scale(cbind(y1,y2,y3))
}


modlina <- function(x){sobol.fun(x)}#
modlinb <- function(x) {atantemp.fun(x,q=q)}##multivariate output
modlinc <- function(x) {ex_fun(x)}##multivariate output

fa <- function(j,x){
	XX <- (x[[j]])
	sobol.fun(XX)
}

fb=function(j,x){
	XX <- x[[j]]
	atantemp.fun(XX,q=100)
}

fc <- function(j,x){
	XX <- x[[j]] #
	modlinc(XX)
}
p0 <- 8

mod1 <- lapply(1:ll,sensitivity_func,myf=Xset1,models="other",CC=diag(1),fun=fa,q=1,p0=8)
q <- 100
mod2 <- lapply(1:ll,sensitivity_func,myf=Xset2,models="other",CC=diag(1,q),fun=fb,q=100,p0=2)
C3 <- matrix(c(1,0.3122,0.3694,0.3112,1,0.0467,0.3694,0.0467,1),nrow=3)
mod3 <- lapply(1:ll,sensitivity_func,myf=Xset3,models="other",CC=C3,fun=fc,q=3,p0=8)

saveRDS(mod1,file=file.path(outdir,"mod1"))
saveRDS(mod2,file=file.path(outdir,"mod2"))
saveRDS(mod3,file=file.path(outdir,"mod3"))
#####
udat1 <- sen_fun(mod1)
udat2 <- sen_fun(mod2)
udat3 <- sen_fun(mod3)

tit <- c("Main indices","Total indices")
k0 <- 1
	M1 <- udat1[[k0]];p0 <- 8
	png(file.path(outdir,"f2b.png"), height =10, width=12,units="in",res=400)
	par(mfrow=c(2,1))
	#par(mai=c(0.25,0.3,0.5,0.1)+0.0)
	xvals=barplot2(M1[[1]],width=1,col="grey",main=tit[1],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.ci=TRUE,ci.u=M1[[2]],ci.l=M1[[3]],plot.grid=TRUE,ylim=range(M1))
	#lines(xvals,unlist(ana_value[[1]]),col="red")
	points(xvals,unlist(ana_value2),col="red",cex=1.5)
	mtext(text=lab0[1:p0],side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
	#
	M1 <- udat1[[k0+1]]
	xvals=barplot2(M1[[1]],width=1,col="grey",main=tit[2],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.ci=TRUE,ci.u=M1[[2]],ci.l=M1[[3]],plot.grid=TRUE,ylim=range(M1))
	#lines(xvals,ana_value2,col="red")
	points(xvals,ana_value2,col="red",cex=1.5)
	mtext(text=lab0[1:p0],side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
	dev.off()


png(file.path(outdir,"f2a0.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(2,1))
xvals=barplot2(x1$S,width=1,col="grey",main=tit[1],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.grid=TRUE)
points(xvals,unlist(ana_value2),col="red",cex=1.5)
mtext(text=lab0[1:p0],side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
xvals=barplot2(x2$S,width=1,col="grey",main=tit[2],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.grid=TRUE)
points(xvals,ana_value2,col="red",cex=1.5)
mtext(text=lab0[1:p0],side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
dev.off()

####test2
png(file.path(outdir,"f3b.png"), height =10, width=12,units="in",res=400)
	par(mfrow=c(2,1))
	for(k0 in 1:2){
	j <- 1
	M1 <- udat2[[k0]];p0 <- 2
	xvals=plot(M1[[j]][j,],col="black",main=tit[k0],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1),type="l",xlab="Outputs",ylab="")
	lines(M1[[1]][j+1,],col="red",main=tit[1],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1))
	lines(M1[[2]][j,],col="green",main=tit[1],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1),type="l")
	lines(M1[[2]][j+1,],col="green",main=tit[1],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1))
	lines(M1[[3]][j,],col="green",main=tit[1],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1),type="l")
	lines(M1[[3]][j+1,],col="green",main=tit[1],cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.5,ylim=range(M1))
	legend("topright",legend=c("x1","x2"),fill=c("black","red"))
	}
	dev.off()

#######test3
namm <- expression("S"["Y"^(1)],"S"["Y"]^(2),"S"["Y"^(3)],"S"["i"]^M,"S"["Ti"]^M,"P"["i"],"P"["Ti"])
dd <- udat3[c(1,3:6)]
d1 <- lapply(1:5,function(i)dd[[i]][[1]]);d2 <- lapply(1:5,function(i)dd[[i]][[2]]);d3 <- lapply(1:5,function(i)dd[[i]][[3]])
d11 <- t(abind(d1));d22 <- t(abind(d2));d33 <- t(abind(d3))

png(file.path(outdir,"f5.png"), height =10, width=12,units="in",res=400)
xvals=barplot2(d11,width=1,col=coll[1:7],cex=1.75,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.ci=TRUE,ci.u=d22,ci.l=d33,plot.grid=TRUE,ylim=range(cbind(d11,d22,d33)),beside=TRUE)
mtext(text=lab0[1:p0],side=1,par("usr")[3],at=xvals[1,],line=0,srt=45,padj=.5,adj=1,cex=1.75,col="blue",cex.lab=1.5,cex.axis=1.1,las=2)
legend("topright",legend=namm,fill=coll[1:7],cex=1.75)
dev.off()








