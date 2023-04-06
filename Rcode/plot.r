########R code for plotting Figures 1, 8, 9, 10, 11, 12, 13, 14, 15

####set the path to the main R code directory eg ../JRSC_paper_code
source("../JRSC_paper_code/main.r",echo=TRUE)
path01 <- path0
source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)

p <- 24 #number of input variables
m <- 18 #no of output 
N <- 2000 # no of original samples
T <- 3 ## time slices
Sc <- 5 # no of scenarios
ndim <- 15 #3*5
NT <- N*ndim #30000
r <- 6 #number of regions
ns <- 600 #number of samples
#Read in data
XX <- data.matrix(fread(file.path(path01,"data_in","Xsen0.txt"),header=TRUE))# ndim times (N*p) matrix
#convert to array
Xsen0 <- array(XX,c(ndim,N,p))
YY <- data.matrix(fread(file.path(path01,"data_in","Ysen0.txt"),header=TRUE))# ndim times (N*m*r) matrix
Ysen0 <- array(YY,c(ndim,N,m,r))
#source(file.path(path01,"Rcode","functions_libraries.r"),echo=TRUE)
rr <- data.matrix(read.table(file.path(path01,"data_in","sce_range.txt"),sep=",",header=FALSE)) #scenario ranges

u1 <- rr[,9:10]#
u2 <- rr[,15:16]#

id <- c(1,7,10) # 
RR <- x <- list()
for(i in 1:ndim){
	x[[i]] <- Xsen0[i,,]
}
X0 <- abind(x,along=1)
X <- Y0[,4:20]
colnames(X) <- NULL

#Extract each scenario from combined data
f1 <- function(i,dat,cond){all(cond[,1]<= dat[i,] & dat[i,] <=cond[,2])}
w1 <- sapply(1:NT,f1,dat=X[,id],u1[id,])#
w2 <- sapply(1:NT,f1,dat=X[,id],u2[id,])#
RR[[1]] <- which(w1==TRUE)
RR[[2]] <- which(w2==TRUE)

########
########unscale output
k <- 1##EU-wide; 2 - 6 denotes the regional averages
Ytest <- list()
for(i in 1:ndim){
	Ytest[[i]] <- Ysen0[i,,,k]##unscaled
}
a1 <- abind(Ytest,along=1)

####Histograms Fig 1a
png(file.path(path01,"plots","f1a.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 1:9){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)
	hist(a1[RR[[1]],i],col="green",main=paste(nam2[i],unit[i]),cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.3)
}
dev.off()

#Fig 1b
png(file.path(path01,"plots","f1b.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 10:18){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)
	hist(a1[RR[[1]],i],col="green",main=paste(nam2[i],unit[i]),cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.3)
}
dev.off()

#Fig14a
png(file.path(path01,"plots","f14a.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 1:9){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)
	hist(a1[RR[[2]],i],col="green",main=paste(nam2[i],unit[i]),cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.3)
}
dev.off()

#Fig14b
png(file.path(path01,"plots","f14b.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 10:18){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)
	hist(a1[RR[[2]],i],col="green",main=paste(nam2[i],unit[i]),cex=1.5,cex.main=1.5,cex.lab=1.5,cex.axis=1.3)
}
dev.off()

###########Fig 9
load(file.path(path01,"data_out","EU_wide","mod"))#EU-wide 
S <- sen_fun(mod)
M1 <- S[[1]]; M2 <- S[[2]];M3 <- S[[3]]; M4 <- S[[4]];M5 <- S[[5]]; M6 <- S[[6]]
#Fig 9 # first-order indices from MSGP
png(file.path(path01,"data_out","EU_wide0","f9.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 1:9){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)	
	xvals <- barplot2(M2[[1]][,i],width=1,col=coll[i],main=nam2[i],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.ci=TRUE,ci.u=M2[[2]][,i],ci.l=M2[[3]][,i],plot.grid=TRUE,ylim=range(M2))
	mtext(text=lab0,side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
}
dev.off()

#Fig 10 # total indices from MSGP
png(file.path(path01,"data_out","EU_wide0","f10.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,3))
for(i in 1:9){
	par(mai=c(0.25,0.3,0.5,0.1)+0.0)	
	xvals <- barplot2(M1[[1]][,i],width=1,col=coll[i],main=nam2[i],cex=1.5,cex.main=1.5,offset=0,cex.lab=1.5,cex.names=1.5,cex.axis=1.5,axisnames=FALSE,plot.ci=TRUE,ci.u=M1[[2]][,i],ci.l=M1[[3]][,i],plot.grid=TRUE,ylim=range(M1))
	mtext(text=lab0,side=1,par("usr")[3],at=xvals,line=0,srt=45,padj=.5,adj=1,cex=1.1,col="blue",cex.lab=1.1,cex.axis=1.1,las=2)
}
dev.off()

#Fig 11 #Generalized indices vs vector projection approach
ylim <- range(M3,M4,M5,M6)
G <- list(M3,M4,M5,M6)
png(file.path(path01,"data_out","EU_wide0","f11.png"), height =10, width=12,units="in",res=400)
par(mfcol=c(2,2))
par(mai=c(0.25,1.5,0.5,0.2)+0.0)
for(k in 1:length(G)){
	xvals <- barplot2((G[[k]][[1]]),las=2,xlab=" ",ylab="Sensitivity indices",cex=1.35,main="Main_sensitivity",cex.lab=1.5,col="green",cex.main=1.5,cex.axis=1.5,cex.names=1.5,cex=1.5,xaxt="n",plot.ci=TRUE,ci.u=G[[k]][[2]],ci.l=G[[k]][[3]],plot.grid=TRUE,ylim=ylim)
	text(xvals,par("usr")[3]-.005,srt=300,adj=1,labels=lab0,xpd=TRUE,cex=1.35,col="blue",srt=45)
} 
dev.off()

#Fig 12
load(file.path(path01,"data_out","scenario_based_sensitivity","sce_mod"))#scenario-based 
G <- list()
for(j in 1:length(sce_mod)){
	S <- sen_fun(sce_mod[[j]])
	M1 <- S[[1]]; M2 <- S[[2]];M3 <- S[[3]]; M4 <- S[[4]];M5 <- S[[5]]; M6 <- S[[6]]
	G[[j]] <-  list(rowMeans(M2[[1]]),rowMeans(M2[[2]]),rowMeans(M2[[3]]))  #M5
}
#
png(file.path(path01,"data_out","scenario_based_sensitivity","ff2.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(2,2))
par(mai=c(0.25,1.5,0.5,0.2)+0.0)
for(k in 1:length(G)){
	xvals=barplot2((G[[k]][[1]]),las=2,xlab=" ",ylab="Sensitivity indices",cex=1.35,main=nam00[1],cex.lab=1.5,col="green",cex.main=1.5,cex.axis=1.5,cex.names=1.5,cex=1.5,xaxt="n",plot.ci=TRUE,ci.u=G[[k]][[2]],ci.l=G[[k]][[3]],plot.grid=TRUE,ylim=range(G))
	text(xvals,par("usr")[3]-.005,srt=300,adj=1,labels=lab0,xpd=TRUE,cex=1.45,col="blue",srt=45)
}
dev.off()

##Fig 8 and 13/main effect
rm(mod)
load(file.path(path01,"data_out","EU_wide","Effect"))#EU-wide 
sim <- seq(-1,1,.1)
xx <- matrix(sim,nrow=length(sim),ncol=p)
ME <- (apply(abind(sapply(Effect,`[[`, "EJ", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2,3,4),mean))
#VE <- (apply(abind(sapply(Effect,`[[`, "VJ", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2,3,4),mean))
##############
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = FALSE)]
cc <- col_vector[-4]
typ <- rep(1:6,each=4) 
coll <- cc[1:p] 
pc <- 1:p

K <- list(c(1:9),c(10:18)); nam_list <- c("f8.png","f13.png")
for(k0 in 1:2){
	J <- K[[k0]]
	png(file.path(path01,"data_out","EU_wide",nam_list[k0]), height =10, width=12,units="in",res=400)
	rr <- range(ME)
	par(oma = c(6.5,1,1,1), mfrow = c(3, 3), mar = c(2, 4, 2, 1))
	for(i in J){
		j <- 1
		plot(xx[,j],ME[j,,i],xlab="Input variable",ylab="Main effects",cex=1.3,cex.main=1.3,col=coll[j],main=nam2[i],cex.lab=1.1,ylim=rr,pch=pc[j],lty=typ[j],lwd=.95)
				for(j in 1:p){
				lines(xx[,j],ME[j,,i],xlab="",ylab="Main effects",cex=1.3,cex.main=1.3,col=coll[j],main=nam2[i],ylim=rr,cex.lab=1.1,pch=pc[j],lty=typ[j],type="o",lwd=.95)
				}
	}
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
	legend("bottom",legend=nam,col=coll,lty=typ,pch=pc, cex =1.2,ncol = 6,lwd=4)
	dev.off()
}

#Fig 15 regional plots
rm(mod)
load(file.path(path01,"data_out","regional_based_sensitivity","reg_mod"))#scenario-based 
G <- list()
t <- 4
for(j in 1:length(reg_mod)){
	S <- sen_fun(reg_mod[[j]])
	G[[j]] <- S[[t]]
}

#f15
png(file.path(path01,"data_out","regional_based_sensitivity","f15.png"), height =10, width=12,units="in",res=400)
par(mfrow=c(3,2))
par(mai=c(0.25,1.5,0.5,0.2)+0.0)
for(k in 1:length(G)){
	xvals=barplot2((G[[k]][[1]]),las=2,xlab=" ",ylab="Sensitivity indices",cex=1.35,main=nam0[k],cex.lab=1.5,col="green",cex.main=1.5,cex.axis=1.5,cex.names=1.5,cex=1.5,xaxt="n",plot.ci=TRUE,ci.u=G[[k]][[2]],ci.l=G[[k]][[3]],plot.grid=TRUE,ylim=range(G))
	text(xvals,par("usr")[3]-.005,srt=300,adj=1,labels=lab0,xpd=TRUE,cex=1.45,col="blue",srt=45)
	text(xvals,par("usr")[3]-.005,srt=300,adj=1,labels=lab0,xpd=TRUE,cex=1.45,col="blue",srt=45)
}
dev.off()
