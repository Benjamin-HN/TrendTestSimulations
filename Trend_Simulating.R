
################## Simulating Trends ##################


################################################################################
###############################      Fig 1       ###############################
################################################################################
##########################   Different beta values   ########################### 
################################################################################
 
rm(list=ls())


source("Process_Simulating.R")
source("Trend_Tests_Codes.R")


antsim <- 20000
n.fail <- 5 # 5,15,30,60 ## Number of failures

##########################

# Define the test functions
tests <- list("Mann-Kendall" = Mann, "Laplace" = La.fail, "Military Handbook" = MH.fail, "Lewis-Robinson sigma 1" = LRfail,"Lewis-Robinson sigma 2" = LRfail.l, "Anderson-Darling type sigma 1" = ADfeil, "Anderson-Darling type sigma 2" = ADfeil.l)


# Weibull parameters
alpha <- 0.1
beta_values <- seq(0.01, 3, by=0.01)
a <- 0.1
b <- 0.75
#b <- 1.0
#b <- 1.5


# Store test statistics and power for each test
testobs <- list()
power <- list()
t.times.TRP.sim <- vector("list", length(beta_values))

################## Run the simulations ##################
for(s in 1:length(beta_values)){ 
  t.times.TRP.sim[[s]] <- vector("list", antsim)
  for(i in 1:antsim){
    t.times.TRP.sim[[s]][[i]] <- trp.wei.fail(n.fail = n.fail, alpha=alpha,beta=beta_values[s], a = a, b = b)
  }
}

# Run the simulations for each of the tests
for(test_name in names(tests)){
  
  # Get the current test function
  test_fun <- tests[[test_name]]
  
  # Store power estimates for each beta value
  power[[test_name]] <- numeric(length(beta_values))
  
  # Calculate for each test
  for(s in 1:length(beta_values)){
    
    #
    testobs <- sapply(t.times.TRP.sim[[s]], test_fun)
    
    # Calculate power for the current test                            
    if (test_name == "Anderson-Darling type sigma 1"||test_name == "Anderson-Darling type sigma 2") {
      power[[test_name]][s] <- mean(testobs > 2.492)
    } else if (test_name == "Military Handbook") {
      power[[test_name]][s] <- mean(testobs < qchisq(0.025, df = 2 * (n.fail - 1)) | testobs > qchisq(0.975, df = 2 * (n.fail - 1)))
    } else {
      power[[test_name]][s] <- mean(abs(testobs) > 1.96)
    }
    
    print(paste("Test:", test_name, "Alpha:", alpha, "Beta:", beta_values[s], "Power:", power[[test_name]][s]))
  }
}


### Plot
colors <- c("red", "blue", "green",  "purple1","purple4", "orange", "orange2")
lwd_values <- c(2, 2, 2, 2, 2, 2, 2)  # Line width
lty_values <- c(1, 1, 1, 1, 5, 1, 5)  # Line type
ymax <- max(unlist(power),na.rm = TRUE)
#ymax <- 0.2
#pdf(file = "Fig1beta_b0.75_Alp0.1_a0.1_nfail05_sim20k.pdf")
plot(beta_values, power[[names(tests)[1]]], type = "l", xlab = "beta", ylab = "Rejection Probability", 
     main = paste("TRP,", " Number of failures =",n.fail),xlim = c(0.01,3), ylim = c(0, ymax), col = colors[1], lwd = lwd_values[1], lty = lty_values[1])
for(t in 2:length(tests)){
  lines(beta_values, power[[names(tests)[t]]], type = "l", col = colors[t], lwd = lwd_values[t], lty = lty_values[t])
}
#abline(h=0.05, lty=2) # A dotted line for 0.05, where the test should be at for RP
### Text box
legend("topright", legend = names(tests),lty=lty_values, lwd=lwd_values, col = colors, cex = 0.8)
#dev.off()






################################################################################
###############################      Fig 2       ###############################
################################################################################
############################   Different b values   ############################ 
################################################################################


rm(list=ls())



source("Process_Simulating.R")
source("Trend_Tests_Codes.R")



antsim <- 20000
n.fail <- 5  # 5,15,30,60 ## Number of failures

##############################   Now with 2 "different" Lewis Robinson and AD

# Define the test functions
tests <- list("Mann-Kendall" = Mann, "Laplace" = La.fail, "Military Handbook" = MH.fail, "Lewis-Robinson sigma 1" = LRfail, "Lewis-Robinson sigma 2" = LRfail.l, "Anderson-Darling type sigma 1" = ADfeil, "Anderson-Darling type sigma 2" = ADfeil.l)


# Weibull parameters
alpha <- 100
#beta_values <- 0.75 #TRP
#beta_values <- 1    #NHPP
beta_values <- 1.5  #TRP
a <- 0.1
b_values <- seq(0.01, 3.5, by=0.01)


# Store test statistics and power for each test
testobs <- list()
power <- list()
t.times.TRP.sim <- vector("list", length(b_values))

# Run the simulations
for(s in 1:length(b_values)){ 
  t.times.TRP.sim[[s]] <- vector("list", antsim)
  for(i in 1:antsim){
    t.times.TRP.sim[[s]][[i]] <- trp.wei.fail(n.fail = n.fail, alpha=alpha,beta=beta_values, a = a, b = b_values[s])
  }
}

# Run the simulations for each of the tests
for(test_name in names(tests)){
  
  # Get the current test function
  test_fun <- tests[[test_name]]
  
  # Store power estimates for each b value
  power[[test_name]] <- numeric(length(b_values))
  
  # Calculate for each test
  for(s in 1:length(b_values)){
    
    #
    testobs <- sapply(t.times.TRP.sim[[s]], test_fun)
    
    # Calculate power for the current test                            
    if (test_name == "Anderson-Darling type sigma 1"||test_name == "Anderson-Darling type sigma 2") {
      power[[test_name]][s] <- mean(testobs > 2.492)
    } else if (test_name == "Military Handbook") {
      power[[test_name]][s] <- mean(testobs < qchisq(0.025, df = 2 * (n.fail - 1)) | testobs > qchisq(0.975, df = 2 * (n.fail - 1)))
    } else {
      power[[test_name]][s] <- mean(abs(testobs) > 1.96)
    }
    
    print(paste("Test:", test_name, "Alpha:", alpha, "Beta:", beta_values,"b:", b_values[s], "Power:", power[[test_name]][s]))
  }
}


### Plot
colors <- c("red", "blue", "green", "purple1","purple4", "orange", "orange2")
lwd_values <- c(2, 2, 2, 2, 2, 2, 2)  # Line width
lty_values <- c(1, 1, 1, 1, 5, 1, 5)  # Line type
ymax <- max(unlist(power),na.rm = TRUE)
#pdf(file = "Fig2b_beta1.5_Alp100_a0.1_nfail05_sim20k.pdf")
plot(b_values, power[[names(tests)[1]]], type = "l", xlab = "b", ylab = "Rejection probability", xlim = c(0.01,3.5),
     main = paste("TRP,", " Number of failures =",n.fail), ylim = c(0, ymax), col = colors[1], lwd = lwd_values[1], lty = lty_values[1]) # NHPP or TRP
for(t in 2:length(tests)){
  lines(b_values, power[[names(tests)[t]]], type = "l", col = colors[t], lwd = lwd_values[t], lty = lty_values[t])
}
# Text box
legend("bottomright", legend = names(tests),lty=lty_values, lwd=lwd_values, col = colors, cex = 0.8)
#dev.off()




################################################################################
###############################      Fig 3       ###############################
################################################################################
####################  Different number of failures values   ####################
################################################################################



rm(list=ls())

source("Process_Simulating.R")
source("Trend_Tests_Codes.R")


antsim <- 20000  ## Number of simulations 
n.fail_values <- seq(5,60, by=1) ## Number of failures
  
##########################

# Define the test functions
tests <- list("Mann-Kendall" = Mann, "Laplace" = La.fail, "Military Handbook" = MH.fail, "Lewis-Robinson sigma 1" = LRfail,"Lewis-Robinson sigma 2" = LRfail.l, "Anderson-Darling type sigma 1" = ADfeil, "Anderson-Darling type sigma 2" = ADfeil.l)


# Weibull parameters
alpha <- 0.1
#beta <- 0.75 ## RP Overdispersed
#beta <- 1.0  ## HPP
beta <- 1.5   ## RP Underdispersed
a <- 0.1
b <- 1.0


# Store test statistics and power for each test
testobs <- list()
power <- list()
t.times.TRP.sim <- vector("list", length(n.fail_values))

###############################################     Run the simulations 
for(s in 1:length(n.fail_values)){ 
  t.times.TRP.sim[[s]] <- vector("list", antsim)
  for(i in 1:antsim){
    t.times.TRP.sim[[s]][[i]] <- trp.wei.fail(n.fail = n.fail_values[s], alpha=alpha,beta=beta, a = a, b = b)
  }
}

# Run the simulations for each of the tests
for(test_name in names(tests)){
  
  # Get the current test function
  test_fun <- tests[[test_name]]
  
  # Store power estimates for each b value
  power[[test_name]] <- numeric(length(n.fail_values))
  
  # Calculate for each test
  for(s in 1:length(n.fail_values)){
    
    #
    testobs <- sapply(t.times.TRP.sim[[s]], test_fun)
    
    # Calculate power for the current test                            
    if (test_name == "Anderson-Darling type sigma 1"||test_name == "Anderson-Darling type sigma 2") {
      power[[test_name]][s] <- mean(testobs > 2.492)
    } else if (test_name == "Military Handbook") {
      power[[test_name]][s] <- mean(testobs < qchisq(0.025, df = 2 * (n.fail_values[s] - 1)) | testobs > qchisq(0.975, df = 2 * (n.fail_values[s] - 1)))
    } else {
      power[[test_name]][s] <- mean(abs(testobs) > 1.96)
    }
    
    print(paste("Test:", test_name, "Alpha:", alpha, "Beta:", beta, "n.fail:",n.fail_values[s], "Power:", power[[test_name]][s]))
  }
}


### Plot
colors <- c("red", "blue", "green",  "purple1","purple4", "orange", "orange2")
lwd_values <- c(2, 2, 2, 2, 2, 2, 2)  # Line width
lty_values <- c(1, 1, 1, 1, 5, 1, 5)  # Line type
ymax <- max(unlist(power),na.rm = TRUE)
##ymax <- 0.2
#pdf(file = "Fig3nfail60_b1_beta1.5_Alp0.1_a0.1_sim20k.pdf")
plot(n.fail_values, power[[names(tests)[1]]], type = "l", xlab = "Number of failures", ylab = "Rejection probability", 
     main = paste("RP (underdispersed)"), ylim = c(0, ymax), col = colors[1], lwd = lwd_values[1], lty = lty_values[1])
for(t in 2:length(tests)){
  lines(n.fail_values, power[[names(tests)[t]]], type = "l", col = colors[t], lwd = lwd_values[t], lty = lty_values[t])
}
abline(h=0.05, lty=2) ## A dotted line for 0.05, where the test should be at for b=1
####### Text box
legend("topright", legend = names(tests),lty=lty_values, lwd=lwd_values, col = colors, cex = 0.8)
#dev.off()



