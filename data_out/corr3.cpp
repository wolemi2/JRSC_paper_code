//#define ARMA_USE_SUPERLU
//#define ARMA_64BIT_WORD 1
#define ARMA_USE_OPENMP
#define ARMA_USE_BLAS
#define ARMA_ALLOW_FAKE_GCC
#include <RcppArmadillo.h>
#include <R_ext/Boolean.h>

#define _UNIX03_SOURCE
#include <stdlib.h>
int setenv(const char *var_name, const char *new_value, int change_flag);

#include <cmath>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
using namespace arma;
using namespace Rcpp;
using std::string;
//using namespace RcppParallel;
using namespace std;
static double const log2pi = std::log(2.0 * M_PI);
const double pi = 3.141592653589793238463 ;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat mycor(const mat MAT_CORR, const rowvec cf,int const t, const string cor,double const nug,double const alpha,double const nu) {
	omp_set_num_threads(t);
// set stacksize, not really necessary unless you are running into memory issues
setenv("OMP_STACKSIZE","200M",1);
    int i,k,j;
	rowvec r(MAT_CORR.n_cols);
	rowvec ress(MAT_CORR.n_cols);
	sp_mat X(MAT_CORR.n_rows,MAT_CORR.n_rows);
	#if defined(_OPENMP) 
       #pragma omp parallel for private(k, i,r,ress)
#endif 
for (k = 0; k < MAT_CORR.n_rows-1; k++) {
	    for (i = k+1; i < MAT_CORR.n_rows; ++i) {
	r = (abs(MAT_CORR.row(k) - MAT_CORR.row(i)))/cf; 
	if(all(r <= 1)){
if(cor=="bohman"){
	ress = ((1-r) % cos(pi*r) + sin(pi*r)/pi);// % (r <= 1);// 0<=ress<=1
}else{
//double alpha=1.0, nu=1.0;
ress = pow(1 - pow(r,alpha),nu);// % (r <= 1);
}
ress.elem( find(ress < 0) ).zeros();
	X(k,i) = prod(ress);
}}}      
X.diag().ones();
X.diag() += nug;	
return(symmatu(X));
return(X);
 }

// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::export]]
List mycor2(const mat XX,int const n, int const n0, const rowvec cf,int const t, const string cor,double const nug,double const alpha,double const nu) {
	omp_set_num_threads(t);
// set stacksize, not really necessary unless you are running into memory issues
setenv("OMP_STACKSIZE","200M",1);
    int i,k;
	rowvec r(XX.n_cols), ress(XX.n_cols);
	sp_mat X(XX.n_rows,XX.n_rows);
#if defined(_OPENMP) 
    //#pragma omp parallel for shared(cf, XX)
#pragma omp parallel for private(k, i,r,ress)
#endif 
for (k = 0; k < XX.n_rows-1; k++) {
	    for (i = k+1; i < XX.n_rows; ++i) {
	r = (abs(XX.row(k) - XX.row(i)))/cf; 
	if(all(r <= 1)){
if(cor=="bohman"){
ress = ((1-r) % cos(pi*r) + sin(pi*r)/pi);// % (r <= 1);
}else{
ress = pow(1 - pow(r,alpha),nu);// % (r <= 1);
}
ress.elem( find(ress < 0) ).zeros();
X(k,i) = prod(ress);
		}}}      
X.diag().ones();
sp_mat X2 = symmatu(X);

sp_mat X1 = X2.submat(0,0, n-1, n-1);
X1.diag() += nug;
sp_mat X0 = X2.submat(n,n, n+n0-1, n+n0-1);
X0.diag() += nug;
sp_mat X01 = X2.submat(0,n,n-1,n+n0-1);
return(List::create(Named("X01") = X01,Named("X0")=X0,Named("X1")=X1));
}

// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::export]]
sp_mat mycor3(NumericMatrix MAT_CORR,NumericMatrix MAT_CORR2,NumericVector const &cf,int const t, const string cor,double const nug,double const alpha,double const nu) {
	omp_set_num_threads(t);
// set stacksize, not really necessary unless you are running into memory issues
setenv("OMP_STACKSIZE","64M",1);
    int i,k;
	NumericVector temp_SUM(MAT_CORR.ncol());
	vec r(MAT_CORR.ncol()), ress(MAT_CORR.ncol());
	sp_mat X(MAT_CORR.nrow(),MAT_CORR2.nrow());
#if defined(_OPENMP) 
   #pragma omp parallel for shared(cf, MAT_CORR2)
#endif
for (k = 0; k < MAT_CORR.nrow()-1; k++) {
	    for (i = 0; i < MAT_CORR2.nrow(); ++i) {
		temp_SUM = abs(MAT_CORR.row(k) - MAT_CORR2.row(i)); 
		if(is_false(any(temp_SUM >= cf))){
		r = temp_SUM/cf;
if(cor=="bohman"){
ress = ((1-r) % cos(pi*r) + sin(pi*r)/pi);// % (r <= 1);
}else{
ress = pow(1 - pow(r,alpha),nu);// % (r <= 1);
}
X(k,i) = prod(ress);
}}}
return(X);
}

