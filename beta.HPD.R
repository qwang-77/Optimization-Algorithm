# bisection to find the smaller x such that f(x) <- lam 
# where lam > 0 and 
# f(x) is the beta pdf with parameter alpha > 1 and beta > 1
rm(list = lm())
beta.lower <- function (lam, alpha, beta, max.iter=100, precision=1e-8) { 
  
  if (alpha <= 1)
    stop("alpha must be greater than one")
  if (beta <= 1)
    stop("beta must be greater than one")
  mod<-(alpha - 1)/(alpha + beta - 2 )
  if (lam > dbeta(mod, alpha, beta) || lam <= 0)
    stop("No solution") 
  
  # initialize 
  intv <- c(0, mod) 
  g <- dgamma(intv, alpha, beta) - lam 
  if(g[1] * g[2] > 0) 
    stop("error in beta.lower, bad initialize value")
  
  for (iter in 1:max.iter) {
    
    c <- 0.5 * (intv[1] + intv[2]) 
    gc <- dbeta(c, alpha, beta) - lam  
    
    if (gc * g[1] < 0) {
      intv[2] <- c
      g[2] <- gc
    }
    else {
      intv[1] <- c
      g[1] <- gc
    }
    
    if (intv[2] - intv[1] < precision)
      break
  }
  
  out <- list(c, gc, iter) 
  names(out) <- c("root", "function.value", "iteration.count")
  
  out
} 


#>---------------------------------Higer Solution--------------------------------------------
#>-------------------------------------------------------------------------------------------
#>-------------------------------------------------------------------------------------------

# bisection to find the larger x such that f(x) <- lam 
# where lam > 0 and 
# f(x) is the beta pdf with parameter alpha > 1 and beta > 1

beta.higher <- function (lam, alpha, beta, max.iter=100, precision=1e-8) { 
  
  if (alpha <= 1)
    stop("alpha must be greater than one")
  if (beta <= 1)
    stop("beta must be greater than one")
  mod<-(alpha - 1)/(alpha + beta - 2 )
  if (lam > dbeta(mod, alpha, beta) || lam <= 0)
    stop("No solution") 
  
  # initialize 
  intv <- c(mod,1) 
  g <- dgamma(intv, alpha, beta) - lam 
  if(g[1] * g[2] > 0) 
    stop("error in beta.lower, bad initialize value")
  
  for (iter in 1:max.iter) {
    
    c <- 0.5 * (intv[1] + intv[2]) 
    gc <- dbeta(c, alpha, beta) - lam  
    
    if (gc * g[1] < 0) {
      intv[2] <- c
      g[2] <- gc
    }
    else {
      intv[1] <- c
      g[1] <- gc
    }
    
    if (intv[2] - intv[1] < precision)
      break
  }
  
  out <- list(c, gc, iter) 
  names(out) <- c("root", "function.value", "iteration.count")
  
  out
} 


#-------------------------------------------HPD Calculation----------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

beta.HPD<-function(alpha, beta, percent = 0.95, precision.int = 1e-5, precision.lam = 1e-5)
  {
  max.iter <- 1/precision.lam
  mod <- (alpha - 1)/(alpha + beta - 2 )
  lam <- dbeta(mod,alpha, beta)/2
  i <- 1
  for (i in 1:max.iter){
    high <- beta.higher(lam, alpha <- alpha, beta <- beta)$root
    low <- beta.lower(lam, alpha <- alpha, beta <- beta)$root
    prob <- pbeta(high,alpha,beta)-pbeta(low, alpha, beta)
    if (abs(prob-percent)<= precision.int){break}
    if (prob-percent>0){lam <- lam + precision.lam}
    if (prob-percent<0){lam <- lam - precision.lam}
    i <- i + 1
  }
  out <- list(low, high,prob)
  names(out) <- c("Lower","Higher","Probability")
  if (i == max.iter){print("Warning: Reached maximum iteration, enlarge precision.lam or decrease precision.int")}
  out
 }