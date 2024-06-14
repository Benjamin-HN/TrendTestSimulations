
################## Real Data Testing ##################

rm(list=ls())


source("Process_Simulating.R")
source("Trend_Tests_Codes.R")
### goftest needed for p-values for Anderson-Darling
if (!require(goftest)) {
  install.packages("goftest")
}
library(goftest)


### USSH times
USSH_ftimes <- c(1.382, 2.990, 4.124, 6.827, 7.472, 
                 7.567, 8.845, 9.450, 9.794, 10.848, 
                 11.993, 12.300, 15.413, 16.497, 17.352, 
                 17.632, 18.122, 19.067, 19.172, 19.299, 
                 19.360, 19.686, 19.940, 19.944)
USSH_tau <- 20

USSH.n <- length(USSH_ftimes)
event_number_USSH <- 1:length(USSH_ftimes) # Sequence from 1 to the length

#### Plot
#pdf(file = "Fig4_USSH_times.pdf")
plot(USSH_ftimes, event_number_USSH, type = "p", main = "Plot of USSH times", xlab = "Time (hours)", ylab = "Cumulative number of failures",col= "black", pch = 19, cex=1.1, xlim = c(0,USSH_tau))
lines(c(0,USSH_ftimes[USSH.n]), c(0,USSH.n), lty=3)
#dev.off()



###### LHD times
LHD_ftimes <- c(16, 39, 71, 95, 98, 110, 114, 
                226, 294, 344, 555, 599, 757, 
                822, 963, 1077, 1167, 1202, 1257, 
                1317, 1345, 1372, 1402, 1536, 
                1625, 1643, 1675, 1726, 1736,
                1772, 1796, 1799, 1814, 1868, 
                1894, 1970)
LHD_tau <- 2000

LHD.n <- length(LHD_ftimes)
event_number_LHD <- 1:length(LHD_ftimes) # Sequence from 1 to the length

#### Plot
#pdf(file = "Fig5_LHD_times.pdf")
plot(LHD_ftimes, event_number_LHD, type = "p", main = "Plot of LHD times", xlab = "Time (hours)", ylab = "Cumulative number of failures",col= "black", pch = 19, cex=1.1, xlim = c(0,LHD_tau))
lines(c(0,LHD_ftimes[LHD.n]), c(0,LHD.n), lty=3)
#dev.off()




### List up the tests
tests <- list(
  "Mann-Kendall" = Mann,
  "Laplace" = La.time, 
  "Military Handbook" = MH.time, 
  "Lewis-Robinson sigma 1" = LRtime.1,
  "Lewis-Robinson sigma 2" = LRtime.2, 
  "Anderson-Darling type sigma 1" = ADtid.1, 
  "Anderson-Darling type sigma 2" = ADtid.2
)


# Function to apply all tests and to calculate power
apply_tests <- function(tvec, tau){
  results <- list()
  power <- list()
  pvalue <- list()
  
  for (test_name in names(tests)){
    test_fun <- tests[[test_name]]
    if(test_name == "Mann-Kendall"){
      testobs <- test_fun(tvec)
    } else{
      testobs <- test_fun(tvec,tau)
    }
    #testobs <- test_fun(tvec, tau) # Can be used if all test have tau
    
    # Calculate power for the current test
    if (test_name == "Anderson-Darling type sigma 1" || test_name == "Anderson-Darling type sigma 2"){
      power[[test_name]] <- mean(testobs > 2.492)
    } else if (test_name == "Military Handbook"){
      n <- length(tvec)
      power[[test_name]] <- mean(testobs < qchisq(0.025, df = 2 * n) | testobs > qchisq(0.975, df = 2 * n))
    } else {
      power[[test_name]] <- mean(abs(testobs) > 1.96)
    }
    
    # Calculate the p-value
    if (test_name == "Anderson-Darling type sigma 1" || test_name == "Anderson-Darling type sigma 2"){
      pvalue[[test_name]] <- 1 - pAD(q = testobs)
    }
      else if (test_name == "Military Handbook"){
      n <- length(tvec)
      pvalue[[test_name]] <- 2 * min(pchisq(testobs, df = 2 * n), 1 - pchisq(testobs, df = 2 * n))
    }
    # else if(test_name == "Mann-Kendall"){}
      
      else{
      pvalue[[test_name]] <- 2 * pnorm(-abs(testobs))
    }
    
    results[[test_name]] <- list("TestResult" = testobs, "Power" = power[[test_name]], "pvalue"= pvalue[[test_name]])
  }
  return(results)
}


### Calculate

# Function to apply all tests to the data and calculate power
USSH_results <- apply_tests(USSH_ftimes, USSH_tau)
LHD_results <- apply_tests(LHD_ftimes, LHD_tau)

USSH_results ### Print the results
LHD_results  ### Print the results



#### Used for the plots to work
# Extract out the power results
USSH_power <- sapply(USSH_results, function(x) x$Power)
LHD_power <- sapply(LHD_results, function(x) x$Power)

# Convert the power results to a data frame for plot
USSH_power_df <- data.frame(Test = names(USSH_power), Power = unlist(USSH_power))
LHD_power_df <- data.frame(Test = names(LHD_power), Power = unlist(LHD_power))

# Extract out the p-value
USSH_p <- sapply(USSH_results, function(x) x$pvalue)
LHD_p <- sapply(LHD_results, function(x) x$pvalue)

# Convert the power results to a data frame for plot
USSH_p_df <- data.frame(Test = names(USSH_p), pvalue = unlist(USSH_p))
LHD_p_df <- data.frame(Test = names(LHD_p), pvalue = unlist(LHD_p))

# Extract the test result for all tests
USSH_test_results <- sapply(USSH_results, function(x) x$TestResult)
LHD_test_results <- sapply(LHD_results, function(x) x$TestResult)


### Used for all test
short_names <- c("Mann","Laplace", "MH", "LR1", "LR2", "AD1", "AD2")
colors <- c("red","blue", "green",  "purple1","purple4", "orange", "orange2")




########### Plot Power

#### Plot the power results for USSH
#pdf(file = "Fig6_USSH_times3.pdf")
plot(USSH_power_df$Power, type = "p", xaxt = "n", xlab = "", ylab = "Test Outcome", main = "USSH Test Outcome", ylim = c(-0.5,1.5),col= colors, pch = 19, cex = 2,yaxt="n")
axis(1, at = 1:length(USSH_power_df$Test),  labels = short_names)
axis(2, at = c(0, 1), labels = c(expression(paste("Not Reject "," ", H[0])), expression(paste("Reject "," ", H[0]))))
#dev.off()


#### Plot the power results for LHD
#pdf(file = "Fig7_LHD_times3.pdf")
plot(LHD_power_df$Power, type = "p", xaxt = "n", xlab = "", ylab = "Test Outcome", main = "LHD Test Outcome",ylim = c(-0.5,1.5),col= colors, pch = 19, cex = 2, yaxt="n")
axis(1, at = 1:length(LHD_power_df$Test),  labels = short_names)
axis(2, at = c(0, 1), labels = c(expression(paste("Not Reject "," ", H[0])), expression(paste("Reject "," ", H[0]))))
#dev.off()




################### Plot Test Result, these 2 plots are not used in thesis

###### USSH

# Create a numeric vector for the x-axis
x_values_USSH <- 1:length(USSH_test_results)

# Plot the test results
#pdf(file = "Fig8_USSH_times3.pdf")
plot(x_values_USSH, USSH_test_results, type = "p", xlab = "Test", ylab = "Test Result", main = "Test Results for USSH", xaxt = "n", ylim = c(-5,30), col= colors, pch = 19, cex =1.5)
axis(1, at = x_values_USSH, labels = short_names)
text(x = x_values_USSH, y = USSH_test_results, labels = round(USSH_test_results, 2), pos = 3, cex = 0.9, col = colors)
#dev.off()


###### LHD


# Create a numeric vector for the x-axis
x_values_LHD <- 1:length(LHD_test_results)
ymax <- max(unlist(LHD_test_results),na.rm = TRUE)


# Plot the test results
#pdf(file = "Fig9_LHD_times3.pdf")
plot(x_values_LHD, LHD_test_results, type = "p", xlab = "Test", ylab = "Test Result", main = "Test Results for LHD", xaxt = "n", ylim = c(-5,ymax+5), col= colors, pch = 19, cex=1.5)
axis(1, at = x_values_LHD, labels = short_names)
text(x = x_values_LHD, y = LHD_test_results, labels = round(LHD_test_results, 2), pos = 3, cex = 0.9, col = colors)
#dev.off()




################### Plot for p-values

###### USSH

# Create a numeric vector for the x-axis
x_values_USSH <- 1:length(USSH_test_results)

# Plot the test results
#pdf(file = "Fig10_USSH_pvalue3.pdf")
plot(USSH_p_df$pvalue, type = "p", xaxt = "n", xlab = "", ylab = "p-value", main = "USSH p-value Results", ylim = c(-0.02,1.02), col = colors, pch = 19, cex = 1.5, xlim = c(0.9, length(USSH_p_df$pvalue) + 0.1))
axis(1, at = 1:length(USSH_p_df$pvalue),  labels = short_names)
text(x = x_values_USSH, y = USSH_p, labels = round(USSH_p, 3), pos = 3, cex = 1.2, col = colors)
#dev.off()



###### LHD

# Create a numeric vector for the x-axis
x_values_LHD <- 1:length(LHD_test_results)

# Plot the test results
#pdf(file = "Fig11_LHD_pvalue3.pdf")
plot(LHD_p_df$pvalue, type = "p", xaxt = "n", xlab = "", ylab = "p-value", main = "LHD p-value Results", ylim = c(-0.02,1.02), col = colors, pch = 19, cex = 1.5, xlim = c(0.9, length(LHD_p_df$pvalue) + 0.1))
axis(1, at = 1:length(LHD_p_df$pvalue),  labels = short_names)
text(x = x_values_LHD, y = LHD_p, labels = round(LHD_p, 3), pos = 3, cex = 1.2, col = colors)
#dev.off()


