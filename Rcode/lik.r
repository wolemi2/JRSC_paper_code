#the likelihood and MCMC functions
lik <- function(cf){
	R <- mycor(Xtrain,cf,nthread,cor,nug,alpha,nu)
	Omegahat <- Matrix::solve(Matrix::crossprod(Xtrain,Matrix::solve(R,Xtrain)) + temp0)
	RZ  <- Matrix::solve(R,Ytrain)
	dett <- determinant(R)$modulus[1]
	rm(R); gc()
	mhat <- Matrix::crossprod(Xtrain,RZ)+temp1
	Mhat <- Matrix::crossprod(t(Omegahat),mhat)
	Shat <- Matrix::crossprod(Ytrain,RZ) + Matrix::crossprod(B0,temp1) + S0 - Matrix::crossprod(Mhat,solve(Omegahat,Mhat))
	post <- -0.5*q*dett - 0.5*q*determinant(Omegahat)$modulus[1] - 0.5*(q+df-1)*determinant(Shat)$modulus[1]
	return(list(post=post,Shat=Shat,Mhat=Mhat,Omegahat=Omegahat))
}

matern_lik <- function(phi){
    #R <- mycor(Xtrain,cf,nthread,cor,nug,alpha,nu)
    R0 <- stationary.taper.cov(Xtrain,V=diag(phi), Covariance="Matern",smoothness=5/2,Taper.args=list(k=2,aRange=theta0,dimension=p))
	R <- as.dgCMatrix.spam(R0) + Diagonal(n,nug)#
	Omegahat <- Matrix::solve(Matrix::crossprod(Xtrain,Matrix::solve(R,Xtrain)) + temp0)
	RZ  <- Matrix::solve(R,Ytrain)
	dett <- determinant(R)$modulus[1]
	rm(R) #gc()
	mhat <- Matrix::crossprod(Xtrain,RZ)+temp1
	Mhat <- Matrix::crossprod(t(Omegahat),mhat)
	Shat <- Matrix::crossprod(Ytrain,RZ) + Matrix::crossprod(B0,temp1) + S0 - Matrix::crossprod(Mhat,solve(Omegahat,Mhat))
	post <- -0.5*q*dett - 0.5*q*determinant(Omegahat)$modulus[1] - 0.5*(q+nus-1)*determinant(Shat)$modulus[1]
	return(list(post=post,Shat=as.matrix(Shat),Mhat=as.matrix(Mhat),Omegahat=as.matrix(Omegahat)))
	}

GP <- function(NIter,dat0,inp0,Phi){
     start.time <- proc.time()
	for(j in 2:NIter){
		ind1 <- j%%thin
		jj <- j[!ind1]
		if (j%%10 == 0) print(j)
    	met <- metropolis(j,SS,Phi,posterior,likelihood,acceptance_prob)
		OMEGA[jj/thin,,] <- as.matrix(met$likelihood$Omegahat)
		THETA[jj/thin,,] <- as.matrix(met$likelihood$Mhat)
		SIGMANU[jj/thin,,] <- as.matrix(met$likelihood$Shat)#Sigma_nu 
		Phi <- met$res
		PHI[jj/thin,] <- cf <- Phi[j,]
		SS <- met$SS
		posterior <- met$posterior
		likelihood <- met$likelihood
	acceptance_prob <- met$acceptance_prob
	}
	est.time <- (proc.time()[3] - start.time[3])/ 3600
      print(paste("Estimated running time:", round(est.time, 2), "hours"))
    result=list(PHI=PHI,SIGMANU=SIGMANU,THETA=THETA,OMEGA=OMEGA)
}