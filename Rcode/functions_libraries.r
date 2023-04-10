############Useful functions and libraries
#0 = RCP2.6; 1 = RCP4.5 ; 2 = RCP8.5
#For RCP2.6 {0 = EC-EARTH_RCA4; 1 = MPI-ESM-LR_REMO; 2 =  NorESM1-M_RCA4} ; For RCP4.5 {0 = HadGEM2-ES_RCA4; 1 = MPI-ESM-LR_CCLM4; 2 =  GFDL-ESM2M_RCA4} ; 
#For RCP8.5 {0 = HadGEM2-ES_RCA4; 1 = CanESM2_CanRCM4; 2 =  IPSL-CM5A-MR_WRF; 3 = GFDL-ESM2M_RCA4} ;
#0 = "SSP1 / We are the world"; 1 = "SSP3 / Icarus"; 2 =  "SSP5 / Should I Stay Or Should I Go" 3 = ;"SSP4 / Riders on the Storm"; 4 = "Baseline"

nam0=c("Alpine","Northern","Atlantic","Continental","Southern")
nam00=c("SSP1","SSP3","SSP4","SSP5")
nam=c("Emission","model","SESs","Population",	
"GreenRed","RumLivesDem","NonRumiLivesDem",	
"StructChange","TechFactor","TechChange",	
"YieldFactor","IrrigEffFactor","GDP",	
"CostsFactor","ImportFactor","BioEneCropDemand",	
"ArableConservLandSc","CropInputsFactor",	
"Coastal","Fluvial","HabitatRecreation","PA_Forestry",	
"PA_Agricultural","FloodProtection")

nam2=c("Food Per Capita","People Flooded 1/100year","Water Exploitation Index","Landuse Intensity Index","Landuse Diversity",
"Biodiversity Vulnerability Index","Timber Production", "Area of Artificial Surfaces"
,"Food Production","Potential Carbon Stock","Irrigation Usage","Intensive Arable","Intensive Grassland","Extensive Grassland",
"Very Extensive Grassland","Unmanaged Land","Managed Forest","Unmanaged Forest")
unit <- c("(kcal/day)","(x100persons)","(_)","(_)","(_)","(_)","(Mt)","(%)","(TJ)","(t/ha)","(mill. mA^3)","(%)","(%)","(%)","(%)","(%)","(%)","(%)")
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

lab <- expression("x"["1"],"x"["2"],"x"["3"],"x"["4"],"x"["5"],"x"["6"],"x"["7"],"x"["8"],"x"["9"],"x"["10"],
"x"["11"],"x"["12"],"x"["13"],"x"["14"],"x"["15"],"x"["16"],"x"["17"],"x"["18"],"x"["19"],"x"["20"],"x"["21"],"x"["22"],"x"["23"],"x"["24"])
lab0 <- lab
coll <- c(rainbow(7),"gray",cm.colors(3),"blueviolet","darkgreen","lightcoral",rev(heat.colors(3)))
#install.packages("myCor_1.0.tar.gz",repos=NULL)
#sourceCpp("corr3.cpp",showOutput=FALSE,verbose=FALSE)

unlink(".RData")
ncore <- nthread <- 6 
library(myCor)
library(gplots)
library(Rcpp)
#library(sensitivity)
library(data.table)
library(parallel)
library(doParallel)
library(RcppArmadillo)
library(adaptMCMC)
library(abind)
library(RColorBrewer)
library(ggplot2)
library(TSP)
library(coda)
#library(rstan)
library(Matrix)
cl <- makeCluster(ncore)
doParallel::registerDoParallel(cl)

tr <- function (a)
{
    i <- 1:nrow(a)
    return(sum(a[cbind(i, i)]))
}


 RMSE <- function(z, zhat) {
        z <- as.matrix(z)
        zhat <- as.matrix(zhat)
        x <- c(z - zhat)
        u <- x[!is.na(x)]
        round(sqrt(sum(u^2)/length(u)), 4)
    }

sen_fun <- function(modd){
par <- c("mean","sd")
R1=R2=R3=R4=R5=R6=list()
M1=M2=M3=M4=M5=M6=list()
for(i in 1:length(par)){
R1[[i]]=(apply(abind(sapply(modd,`[[`, "RES1", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2,3),par[i]))
R2[[i]]= (apply(abind(sapply(modd,`[[`, "RES2", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2,3),par[i]))
R3[[i]]=apply(abind(sapply(modd,`[[`, "RES3", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2),par[i])
R4[[i]]=apply(abind(sapply(modd,`[[`, "RES4", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2),par[i])
R5[[i]]=apply(abind(sapply(modd,`[[`, "P1", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2),par[i])
R6[[i]]=apply(abind(sapply(modd,`[[`, "P2", simplify = FALSE, USE.NAMES = FALSE),along=0),c(2),par[i])
}
i=1
M1[[i]] = R1[[i]]; M1[[i+1]] = R1[[i]] + R1[[i+1]]; M1[[i+2]] = R1[[i]] - R1[[i+1]]
M2[[i]] = R2[[i]]; M2[[i+1]] = R2[[i]] + R2[[i+1]]; M2[[i+2]] = R2[[i]] - R2[[i+1]]
M3[[i]] = R3[[i]]; M3[[i+1]] = R3[[i]] + R3[[i+1]]; M3[[i+2]] = R3[[i]] - R3[[i+1]]
M4[[i]] = R4[[i]]; M4[[i+1]] = R4[[i]] + R4[[i+1]]; M4[[i+2]] = R4[[i]] - R4[[i+1]]
M5[[i]] = R5[[i]]; M5[[i+1]] = R5[[i]] + R5[[i+1]]; M5[[i+2]] = R5[[i]] - R5[[i+1]]
M6[[i]] = R6[[i]]; M6[[i+1]] = R6[[i]] + R6[[i+1]]; M6[[i+2]] = R6[[i]] - R6[[i+1]]
list(M1,M2,M3,M4,M5,M6)
}


##functions
truncpow <- function(d, cutoff, alpha = 1){
  stopifnot(alpha %in% c(1, 3/2, 5/3))
  nu <- switch(match(alpha, c(1, 3/2, 5/3)), 1, 2, 3)
  r <- d/cutoff
  return((1-r^alpha)^nu * (r <= 1))
}
add <- function(x) Reduce("+", x)

find.tau <-
function (den, dim) 
{
    f.obj <- function(x, dim, den) {
        estden <- x^dim * (2 - x)^dim
        return((estden - den)^2)
    }
    opt <- optimize(f.obj, lower = 0, upper = 1, dim = dim, den = den)
    return(opt$minimum)
}

ff2 <- function(x){
	xx=x
	xmin=apply(x,2,min)
	xmax=apply(x,2,max)
	for(i in 1:ncol(x)){
		xx[,i]=2*((x[,i]-xmin[i])/(xmax[i]-xmin[i])) -1
		xx[is.nan(xx)] <- 0
	}
	return(xx)
}

ff <- function(x){
	xx=x
	#xmin=apply(x,2,min)
	#xmax=apply(x,2,max)
	#grand mean across scenario
	xmin=apply(x,3,min)
	xmax=apply(x,3,max)
	for(i in 1:p){
		#xx[,,i]=(x[,,i]-xmin[i])/(xmax[i]-xmin[i])
		xx[,,i]=2*((x[,,i]-xmin[i])/(xmax[i]-xmin[i])) -1
		xx[is.nan(xx)] <- 0
	}
	return(xx)
}

dir_create <- function(path0){{if(!dir.exists(path0)){
    dir.create(path0, recursive = TRUE)
  }}}
   
  
polySet <- function(p, d, w, m)
{
    tmp <- c(p, d, w, m)
    if (any(tmp < 0))
        return(NULL)
    if (any(tmp == 0))
        return(matrix(0, 1, p))
    mdm <- 0:min(d, m)
    robj <- vector("list", length(mdm))
    for (i in seq(along = mdm)) {
        j <- mdm[i]
        robj[[i]] <- cbind(as.vector(j), Recall(p - 1, d - j,
            w - (j > 0), m))
    }
    do.call("rbind", robj)
}

makeScope <- function(pset, xnms)
{
    stopifnot(is.matrix(pset), all(pset >= 0), all(pset == round(pset)))
    stopifnot(is.vector(xnms), length(xnms) == ncol(pset), is.character(xnms))
    "makeTerm" <- function(expt, nm) {
        tt <- character(0)
        if (expt > 0)
            for (i in 1:expt) if (i == 1) {
                tt <- nm
                ones <- "*1"
            }
            else {
                rest <- sprintf("I(%s%s)", nm, ones)
                tt <- sprintf("%s:%s", tt, rest)
                ones <- sprintf("%s%s", ones, "*1")
            }
        tt
    }
    terms <- tapply(pset, row(pset), function(pp) {
        mons <- lapply(seq(along = xnms), function(i) makeTerm(pp[i],
            xnms[i]))
        mons <- mons[sapply(mons, function(m) length(m) > 0)]
        paste(unlist(mons), collapse = ":")
    })
    terms[terms == ""] <- "1"
    paste(terms, collapse = " + ")
}

#######extra packages require for fitting sensitivity analysis
pack <- c("sensitivity","Rcpp","Matrix","myCor","abind","foreach")

#extra functions to import 
extra_func <- c("add","fa","PHI","Mhat1","Omegahat1","Shat1","nug","Xset1","p0","nam","quad","Xtrain","Ytrain",
"nthread","p","cf","nug","mhat","shat","omegahat","alpha","nu","tr")

quad <- formula(~Emission + model + SESs + Population + GreenRed + RumLivesDem +
    NonRumiLivesDem + StructChange + TechFactor + TechChange +
    YieldFactor + IrrigEffFactor + GDP + CostsFactor + ImportFactor +
    BioEneCropDemand + ArableConservLandSc + CropInputsFactor +
    Coastal + Fluvial + HabitatRecreation + PA_Forestry + PA_Agricultural +
    FloodProtection)