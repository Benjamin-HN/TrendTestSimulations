
################## Code of each tests ##################



######  Laplace Test Time ###### 


La.time <- function(tvec,tau){
  n <- length(tvec)
  La <- sqrt(12*n) *(sum(tvec[1:n])-n*tau/2) /(n*tau)
  return(La)
}


######  Laplace Test Failure ###### 



La.fail <- function(tvec){
  n <- length(tvec)
  La <- sqrt(12*(n-1)) *(sum(tvec[1:(n-1)])-(n-1)*tvec[n]/2) /((n-1) * tvec[n])
  return(La)
}




######  Military Handbook Time ######


MH.time <- function(tvec,tau){
  
  # The MH test
  MH <- 2* sum(log(tau/tvec))
  
  return(MH)
}




######  Military Handbook Failure ###### 



MH.fail <- function(tvec){
  n <- length(tvec) # Find number of elements
  Tn <- tvec[n] # Time of the last event
  
  # The MH test
  MH <- 2* sum(log(Tn/tvec[-n]))
  
  
  return(MH)
}




###### Lewis-Robinson Time ###### 


LRtime <- function(tvec,tau,sigma="s"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  LR <- (mu/sdest) * sqrt(12)/(tau*sqrt(n)) * (sum(tvec) - n*tau/2)
  return(LR)
}




###### Lewis-Robinson Time sigma l ###### 


LRtime.l <- function(tvec,tau,sigma="l"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  LR <- (mu/sdest) * sqrt(12)/(tau*sqrt(n)) * (sum(tvec) - n*tau/2)
  return(LR)
}





###### Lewis-Robinson Failure ######


LRfail <- function(tvec,sigma="s"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  LR <- sqrt(12*mu^2/(n*sdest^2*tvec[n]^2))*(sum(tvec)-(n+1)*tvec[n]/2)
  return(LR)
}

###### Lewis-Robinson Failure, with sigma=l ######


LRfail.l <- function(tvec,sigma="l"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  LR <- sqrt(12*mu^2/(n*sdest^2*tvec[n]^2))*(sum(tvec)-(n+1)*tvec[n]/2)
  return(LR)
}





##### Mann test


Mann <- function(tvec){
  n <- length(tvec)
  xvec = diff(c(0,tvec))
  
  ## Calculate the Mann-Kendall formula
  S <- sum(sapply(1:(n-1),function(i){
    sum(ifelse(xvec[(i+1):n]-xvec[i]>0,1,0))
  }))
  
  ## Calculate the Expectation value and Variance
  E.S <- n*(n-1)/4
  VarS <- (2*n^3+3*n^2-5*n)/72
  
  ## Calculate Z statistic
  Z <- (S - E.S) / sqrt(VarS)

  return(Z)
}



###### Anderson-Darling type Time ######


ADtid <- function(tvec,tau,sigma="s"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  konst <- mu^3/(sdest^2*tau)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}



###### Anderson-Darling type Time Sigma l ######


ADtid.l <- function(tvec,tau,sigma="l"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  konst <- mu^3/(sdest^2*tau)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}



###### Anderson-Darling type Failure ######


ADfeil <- function(tvec,sigma="s"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  konst <- n*mu^2/sdest^2
  iseq <- 2:(n-1)
  imseq <- 1:(n-2)
  Niseq <- (n-1):2
  Nimseq <- (n-2):1
  qi <- tvec[2:(n-1)]/tvec[n]-iseq*xvec[2:(n-1)]/tvec[n]
  ri <- n*xvec[2:(n-1)]/tvec[n]-1
  qn <- 1-n*xvec[n]/tvec[n]
  r1 <- n*xvec[1]/tvec[n]-1
  rn <- n*xvec[n]/tvec[n]-1
  forsteledd <- r1^2*log(n/(n-1))-r1^2/n 
  sisteledd <- qn^2*log(n/(n-1))-rn^2/n
  sumledd <- sum(qi^2*log(iseq/imseq)+(qi+ri)^2*log(Niseq/Nimseq)-ri^2/n)
  AD <- konst*(forsteledd+sumledd+sisteledd)
  return(AD)
}


######### Anderson-Darling type Failure sigma l ######

ADfeil.l <- function(tvec,sigma="l"){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  }
  else
    sdest <- sqrt(var(xvec))
  konst <- n*mu^2/sdest^2
  iseq <- 2:(n-1)
  imseq <- 1:(n-2)
  Niseq <- (n-1):2
  Nimseq <- (n-2):1
  qi <- tvec[2:(n-1)]/tvec[n]-iseq*xvec[2:(n-1)]/tvec[n]
  ri <- n*xvec[2:(n-1)]/tvec[n]-1
  qn <- 1-n*xvec[n]/tvec[n]
  r1 <- n*xvec[1]/tvec[n]-1
  rn <- n*xvec[n]/tvec[n]-1
  forsteledd <- r1^2*log(n/(n-1))-r1^2/n 
  sisteledd <- qn^2*log(n/(n-1))-rn^2/n
  sumledd <- sum(qi^2*log(iseq/imseq)+(qi+ri)^2*log(Niseq/Nimseq)-ri^2/n)
  AD <- konst*(forsteledd+sumledd+sisteledd)
  return(AD)
}







############## Without using sigma in function call, for Time ##############
### This is used for Real_data_testing.R #####


###### Lewis-Robinson Time Sigma 1 ###### 

LRtime.1 <- function(tvec,tau){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  sdest <- sqrt(var(xvec))
  LR <- (mu/sdest) * sqrt(12)/(tau*sqrt(n)) * (sum(tvec) - n*tau/2)
  return(LR)
}




###### Lewis-Robinson Time Sigma 2 ###### 

LRtime.2 <- function(tvec,tau){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  xdiff <- xvec[2:n]-xvec[1:(n-1)]
  sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  LR <- (mu/sdest) * sqrt(12)/(tau*sqrt(n)) * (sum(tvec) - n*tau/2)
  return(LR)
}




###### Anderson-Darling type Time Sigma 1 ######

ADtid.1 <- function(tvec,tau){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  sdest <- sqrt(var(xvec))
  konst <- mu^3/(sdest^2*tau)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}



###### Anderson-Darling type Time Sigma 2 ######

ADtid.2 <- function(tvec,tau){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  mu <- mean(xvec)
  xdiff <- xvec[2:n]-xvec[1:(n-1)]
  sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
  sdest <- sqrt(var(xvec))
  konst <- mu^3/(sdest^2*tau)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}



