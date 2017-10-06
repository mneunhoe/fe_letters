#################################
# Auxilliary Functions for Code #
#################################

power.function <- function(alpha=0.05,assumed.ate=NULL,var.assumed.diff=NULL,n.pairs=NULL){
  df <- n.pairs-1
  beta <- 1+pt(qt(alpha/2,df),df=df,ncp=(assumed.ate/sqrt(var.assumed.diff))*sqrt(n.pairs))-pt(-qt(alpha/2,df),df=df,ncp=(assumed.ate/sqrt(var.assumed.diff))*sqrt(n.pairs))
  return(beta)
}


power.function.neq <- function(alpha=0.05,assumed.ate=NULL,n.pairs=NULL,w.k=NULL,dif.k=NULL,n=NULL){ 
  df <- n.pairs-1
  beta <- 1+pt(qt(alpha/2,df),df=df,ncp=((n*assumed.ate)/sqrt(n.pairs*var(w.k*dif.k))))-pt(-qt(alpha/2,df),df=df,ncp=((n*assumed.ate)/sqrt(n.pairs*var(w.k*dif.k))))
  return(beta)
} 


get.z.star <- function(Z, ...){
  
  desmat.out <- desmat.sanitize(Z, blockvar, clustvar)
  
  
  permclus <- replicate(1, do.call(randfun.default, 
                                   list(desmat.out)))
  
  z.star <- as.integer(permclus)
  
  
  z.star
}


get.z.star1 <- function(Z, ...){
  
  desmat.out <- desmat.sanitize(Z, blockvar1, clustvar1)
  
  
  permclus <- replicate(1, do.call(randfun.default, 
                                   list(desmat.out)))
  
  z.star <- as.integer(permclus)
  
  
  z.star
}


