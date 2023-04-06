##########metropolis function
metropolis <- function(j,SS,Phi,posterior,likelihood,acceptance_prob){
u <- rnorm(p)
prop <- (Phi[j-1,]) + c(u %*%SS)
#if (any(prop < 0 | sum(prop) > mc)){
if (any(prop < 0) | sum(prop) > mc){
    Phi[j,] <- Phi[j-1,]
	} else{
	likelihood <- lik(prop)
posterior_proposal <- likelihood$post #+ prior
acceptance_prob <- min(1,exp(posterior_proposal - posterior))
        if (runif(1) < acceptance_prob) {
        Phi[j,] <- prop
   posterior <- posterior_proposal
      }else{
      Phi[j,]  <- Phi[j-1,]
  acceptance_prob <- 0
}
}
if(j < burn) {
      SS <- ramcmc::adapt_S(SS, u, acceptance_prob, n=j,target = 0.234, gamma = 2/3)
    }
          #}
      list(res=Phi,SS=SS,likelihood=likelihood,posterior=posterior,acceptance_prob=acceptance_prob)
   }

