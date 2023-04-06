#########configuration parameters
############ Fixed parameters. Do not alter
p <- 24
m <- q <- 18 #no of output 
ndim <- 15
N <- 2000
ns <- 600 #number of samples
nr <- r <- 6 #number of region
n <- ns*ndim #nrow(Xtrain)
cor <- "truncpow"
#MCMC initialization can be altered
thin <- 25
NIter <- 50000
burn <- 0.5*NIter
		bb <- nburn <- 0.50*NIter/thin###to be removed
		N0 <- niter <- NIter/thin
		delta0 <- -q+1
nus <- df <- delta0+n
B0 <- 0.1*matrix(runif(p*q),nrow=p,ncol=q)
L0 <- 1*diag(p)
S0 <- diag(runif(q),nrow=q,ncol=q)
PHI <- matrix(0, nrow=niter,ncol=p)
Phi <- matrix(NA,NIter,(p))
SIGMAE <- rep(0, niter)
THETA <- array(0,c(niter,p,q))
SIGMANU <- array(0,c(niter,q,q))
OMEGA <- array(0,c(niter,p,p))
temp0 <- Matrix(solve(L0))
temp1 <- Matrix(solve(L0,B0))
nug <- .02; nu <- 2; alpha <- 3/2
## Estimate of needed cutoff based on uniform sampling
mc <- find.tau(den = 1-sparsity, dim = p) * p   ###maximum cutoff
cf <- rep(mc/p/1.5,p)
SS <- (2.38/sqrt(p)) * chol(diag(cf^2/9, nrow = p))
Phi[1,] <- cf
acceptance_prob <- 1 #0.5
##############
listof_region <- c("EU_wide","Alpine","Northern","Atlantic","Continental","Southern")
sp_levels <- c(.99,.95,.90,.80)
