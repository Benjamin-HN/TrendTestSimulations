
################## Process simulating ##################

#rm(list=ls())


########### HPP ###########

hpp <- function(n, lambda, t.end){
  #Generate inter arrival times
  x.times <- rexp(n=n, rate = lambda)
  
  # Calculate arrival times
  t.times <- cumsum(x.times)
  
  while(max(t.times)<t.end){
    x.times <- c(x.times,rexp(rexp(n=n, rate = lambda)))
    t.times <- cumsum(x.times)
  }
  
  # Stop at t.end
  t.times <- t.times[t.times<= t.end]
  
  return(t.times)
}




########### NHPP Inversion ###########


# Define the cumulative intensity function
Lambda <- function(t) {
  return(t^2)
}

# Define the inverse of the cumulative intensity function
Lambda_inv <- function(u) {
  return(sqrt(u))
}

nhpp.inv <- function(n, Lambda_inv, t.end){
  #Generate inter arrival times, HPP
  x.times <- rexp(n=n, rate = 1)
  
  # Calculate arrival times, HPP
  t.times.HPP <- cumsum(x.times)
  
  # Transform the HPP into a NHPP, using the inverse function
  t.times <- Lambda_inv(t.times.HPP)
  
  while(max(t.times)<t.end){
    # Generate a new inter-arrival time for HPP
    x.new <- rexp(n=n, rate = 1)
    # Calculate the new arrival time for HPP
    x.times <- c(x.times,x.new)
    t.new.HPP <- cumsum(x.times)
    
    # Transform the new arrival time into a NHPP time
    t.new <- Lambda_inv(t.new.HPP)
    t.times <- t.new
  }
  
  # Stop at t.end
  t.times <- t.times[t.times <= t.end]
  
  return(t.times)
}

#print(nhpp.inv(n = 200, Lambda_inv, t.end = 10))





########### NHPP Thinning ###########


nhpp.thin <- function(n, lambda.max, t.end, lambda.t){
  # Generate HPP with a maximal rate
  t.times <- rexp(n = n, rate = lambda.max)
  
  # Calculate arrival times
  x.times <- cumsum(t.times)
  
  # Stop at t.end
  x.times <- x.times[x.times <= t.end]
  
  #Thinning
  x.times <- x.times[runif(length(x.times)) <= sapply(x.times, lambda.t) / lambda.max]
  
  
  return(x.times)
}

# Define lambda.t function
lambda.t <- function(t) {
  return(2*t)
}






########### RP (gamma) ###########



rp.gamma <- function(n, shape, rate, t.end){
  #Generate inter arrival times
  x.times <- rgamma(n=n, shape = shape, rate = rate)
  
  # Calculate arrival times
  t.times <- cumsum(x.times)
  
  while(max(t.times)<t.end){
    x.times <- c(x.times, rgamma(n=n, shape = shape, rate = rate))
    t.times <- cumsum(x.times)
  }
  # Stop at t.end
  t.times <- t.times[t.times <= t.end]
  
  
  return(t.times)
}
## Print
#print(rp.gamma(n = 1000, shape = 0.8, rate = 1, t.end = 10))




########### TRP Weibull ###########

# Define the inverse of the cumulative intensity function
Lambda_inv <- function(u) {
  return(sqrt(u))
}

# TRP using inversion method
trp.wei <- function(n, a, b, t.end){
  
  u <- runif(n)
  # Define the inverse of the CDF of the Weibull distribution
  x.times <-  ((-log(1 - u))/a)^(1/b) 
  
  # Calculate arrival times
  t.times.RP <- cumsum(x.times)
  
  # Transform the RP into a TRP, using the inverse function
  t.times <- Lambda_inv(t.times.RP)
  
  while(max(t.times)<t.end){
    
    # Generate a new inter-arrival time for RP
    x.new <- ((-log(1 - u))/a)^(1/b)
    
    # Calculate the new arrival time for RP
    x.times <- c(x.times,x.new)
    t.new.RP <- cumsum(x.times)
    # Transform the new arrival time into a TRP time
    t.new <- Lambda_inv(t.new.RP)
    t.times <- t.new
  }
  
  
  # Stop at t.end
  t.times <- t.times[t.times <= t.end]
  
  return(t.times)
}

### Print
#trp.wei(1000,a=1,b=4,t.end = 10)






########### TRP Gamma ###########


# Define the cumulative intensity function
Lambda <- function(t) {
  return(t^2)
}

# Define the inverse of the cumulative intensity function
Lambda_inv <- function(u) {
  return(sqrt(u))
}

#set.seed(123)
trp.gamma <- function(n, shape, rate, t.end, Lambda_inv){
  # Generate inter-arrival times, RP
  x.times <- rgamma(n=n, shape = shape, rate = rate)
  
  # Calculate arrival times, RP
  t.times.RP <- cumsum(x.times)
  
  # Transform the RP into a TRP using the inverse
  t.times <- Lambda_inv(t.times.RP)
  
  while(max(t.times)<t.end){
    # Generate a new inter-arrival time
    x.new <- rgamma(n=n, shape = shape, rate = rate)
    # Calculate the new arrival time
    t.new <- sum(x.times, x.new)
    # Transform the new arrival time into a TRP time
    t.new <- Lambda_inv(t.new)
    
    t.times <- t.new
  }
  
  # Stop at t.end
  t.times <- t.times[t.times <= t.end]
  
  return(t.times)
}

## Test the function
#print(trp.gamma(n = 200, shape = 1, rate = 1, t.end = 10, Lambda_inv = Lambda_inv))








######################## Failure truncated ########################



########### HPP ###########



hpp.fail <- function(n.fail, lambda){
  # Generate inter-arrival times
  x.times <- rexp(n=n.fail, rate = lambda)
  
  # Calculate arrival times
  t.times <- cumsum(x.times)
  
  return(t.times)
}

#print(hpp.fail(n.fail = 12, lambda = 2))



########### NHPP ###########

# Define the cumulative intensity function
Lambda <- function(t) {
  return(t^2)
}

# Define the inverse of the cumulative intensity function
Lambda_inv <- function(u) {
  return(sqrt(u))
}

nhpp.inv.fail <- function(n.fail, Lambda_inv){
  #Generate inter arrival times, HPP
  x.times <- rexp(n=n.fail, rate = 1)
  
  # Calculate arrival times, HPP
  t.times.HPP <- cumsum(x.times)
  
  # Transform the HPP into a NHPP, using the inverse function
  t.times <- Lambda_inv(t.times.HPP)
  
  # Stop at t.end
  # t.times <- t.times[t.times <= t.end]
  
  return(t.times)
}

#print(nhpp.inv.fail(n.fail = 20, Lambda_inv = Lambda_inv))







########### RP ###########

rp.gamma.fail <- function(n.fail, shape, rate){
  #Generate inter arrival times
  x.times <- rgamma(n=n.fail, shape = shape, rate = rate)
  
  # Calculate arrival times
  t.times <- cumsum(x.times)
  
  # Stop at t.end
  #t.times <- t.times[t.times <= t.end]
  
  
  return(t.times)
}

## Print
#print(rp.gamma.fail(n.fail = 20, shape = 0.8, rate = 1))




########### TRP Gamma ###########

# Define the cumulative intensity function
Lambda <- function(t) {
  return(t^2)
}

# Define the inverse of the cumulative intensity function
Lambda_inv <- function(u,a,b) {
  return((u/a)^(1/b))
#  return(u)
}

#set.seed(123)
trp.gamma.fail <- function(n.fail, shape, rate, Lambda_inv,a,b){
  # Generate inter-arrival times, RP
  x.times <- rgamma(n=n.fail, shape = shape, rate = rate)
  
  # Calculate arrival times, RP
  t.times.RP <- cumsum(x.times)
  
  # Transform the RP into a TRP using the inverse
  t.times <- Lambda_inv(t.times.RP,a,b)
  
  
  return(t.times)
}

### Test the function
#trp.gamma.fail(n.fail = 50, shape = 1, rate = 1, Lambda_inv = Lambda_inv)



########### TRP Weibull ###########



# TRP using inversion method
trp.wei.fail <- function(n.fail, alpha, beta,a,b){
  
  u <- runif(n.fail)
  # Define the inverse of the CDF of the Weibull distribution
  x.times <-  ((-log(1 - u))/alpha)^(1/beta) 
  
  # Calculate arrival times
  t.times.RP <- cumsum(x.times)
  
  # Transform the RP into a TRP, using the inverse function
  t.times <- Lambda_inv(t.times.RP,a,b)
  
  
  return(t.times)
}


## Print
#trp.wei.fail(n.fail=10,alpha=1,beta=4,a=0.1,b=1)


