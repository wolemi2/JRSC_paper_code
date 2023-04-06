#this function predict the fitted MSGP models
pred <- function(xx0,cf,mhat,omegahat,shat){
	n <- nrow(Xtrain);n0=nrow(xx0)
	rr <- mycor2(rbind(Xtrain,xx0),n,n0,cf,nthread,cor,nug,alpha,nu)
	R <- rr$X1
	F1 <- Matrix::solve(R,Ytrain-Matrix::crossprod(t(Xtrain),mhat))
	F2 <- Matrix::solve(R,Xtrain)
	mm <- as.matrix(Matrix::crossprod(t(xx0),mhat) + Matrix::crossprod(rr$X01,F1))
	tq <- Matrix::t(xx0 - crossprod(rr$X01,F2))
	vv <- as.matrix(rr$X0 - Matrix::crossprod(rr$X01, Matrix::solve(R,rr$X01))+
	Matrix::crossprod(tq, Matrix::crossprod(omegahat,tq)))
	V0 <- diag(vv)%*%t(diag(shat))/df
	mu  <- mm + sqrt(V0)* matrix(rt(n=n0*q,df=nus),n0,q)
	return(list(mm=mm,V0=V0,mu=mu))
}