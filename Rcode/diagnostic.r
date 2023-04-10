###Routine for plotting Diagnostics in Figures 6, 7, 16 and 17.
set.seed(10)
#set the path to the main R code directory eg ../JRSSC_paper_code
source("main.r",echo=TRUE)
path01 <- path0

path1 <- file.path(path01,"data_out","EU_wide")#
outdir <- file.path(path1,"plots_1")
outdir2 <- file.path(path1,"plots_2")
dir_create(outdir)
dir_create(outdir2)

load(file.path(path1,"VES00"))
load(file.path(path1,"VES0"))
load(file.path(path1,"Ytest"))
n0 <- nrow(Ytest)/ndim #no of test data point
#
EMU <- array(VES0,c(n0,ndim,m))
VAR <- array(VES00,c(n0,ndim,m))
test <- array(Ytest,c(n0,ndim,m))

g <- c(14,8) #WAW/SSP1 and ROS
nn <- c(4,4)
ff <- c("f6.png","f7.png")
for(i in 1:2){
num <- nn[i]
print(num)
ss <- sample(1:n0,500)
emu1 <- EMU[ss,g[i],] #WAW 2050
var1 <- VAR[ss,g[i],] 
iap <- test[ss,g[i],]

png(file.path(outdir,ff[i]), height =10, width =12, units ="in", res =400)
par(mfrow=c(3,3))
for(j in 1:9){
#id <- order(var1[,j])
id <- order(emu1[,j])
out <- match(boxplot(var1[id,j],plot=FALSE)$out,var1[id,j])
upper <- (emu1[id,j] + 2*sqrt(var1[id,j]))[-out]
lower <- (emu1[id,j] - 2*sqrt(var1[id,j]))[-out]
iapp <- iap[id,j][-out]
emu <- emu1[id,j][-out]
r0 <- range(iapp,emu)
res <- cbind(iapp,emu,lower,upper)
rr <- c(r0[1]-2,r0[2]+2)
par(mar=c(5.5,4.5,1.9,1.9)+0.0)
plot(iapp,main=paste(nam2[j]),cex.lab=1.5,cex=.4,cex.main=1.5,cex.axis=1.5,ylab="Prediction",xlab="No of observation",ylim=c(-5,5),type="p")
lines(upper,col="red",cex=1.5,lwd=2.5)
lines(lower,col="red",cex=1.5,lwd=2.5)
lines(emu,col="green",cex=1.5,lwd=2.5)
legend("topleft",legend=c("IAP2 simulation","Prediction", "95% C.I"),fill=c("black","green","red"),cex=1.5)
}
dev.off()
}

#########trace plots_1
load(file.path(path1,"result1"))
PHI <- result$PHI
tvec <- c(1:24)

png(file.path(outdir,"f16.png"), height =10, width=12,units="in",res=400)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(6,4))
for(j in 1:p){
	hist(PHI[(bb+1):N0,j],col="green",xlab="",prob=TRUE,main=substitute(paste('Density of ',tau[a]),list(a=tvec[j])),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)		
}
dev.off()

png(file.path(outdir,"f17.png"), height =10, width=12,units="in",res=400)
par(mar = c(4, 4, 3, 1) + 0.1, cex = 1.2)
par(mfrow=c(6,4))
for(j in 1:p){
	plot.ts(PHI[(bb+1):N0,j],main=substitute(paste('',tau[a]),list(a=tvec[j])),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,ylab="")
}
dev.off()

#table 7
library(xtable)
s <- sample(1:N0,N0/2)
xtable(gelman.diag(list(as.mcmc(PHI[s,]),as.mcmc(PHI[-s,])),transform=TRUE,autoburnin=FALSE)$psrf)
