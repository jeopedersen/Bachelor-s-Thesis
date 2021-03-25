# Pair Trading S&P 500 - Bachelor's Project Spring 2021
# Author: Jens Pedersen
# Acknowledgements: Anton Rask Lundborg Hansen

library(urca)
library(vars)
library(VARtests)
library(tsDyn)
library(forecast)
library(tidyverse)
library(tseries)
library(Rfast)
setwd("C:/Users/jop.msc/Google Drev/KU/BSc i Økonomi/6. Semester/Bachelorprojekt/Datasæt")

# Load dataset
data <- read.csv("SP500 Log.csv", sep = ";", stringsAsFactors = F, dec = ",")
dates <- data[, 1]
data <- data[, -1]
rownames(data) <- dates
data <- data[1:2513,] #Omitting last three rows due to data corruption

# Find all pairs
pairs <- combn(1:(dim(data)[2]), 2)



misspecification_autocorrelation_detected <- numeric(ncol(pairs))
misspecification_arch_detected <- numeric(ncol(pairs))
misspecification_normality_detected <- numeric(ncol(pairs))
no_cointegration_rejected <- numeric(ncol(pairs)) 
ca_jo_trace_test_stats <- numeric(ncol(pairs))
ca_jo_trace_crit_vals <- numeric(ncol(pairs))
selected_lag <- numeric(ncol(pairs))
used_lag <- numeric(ncol(pairs))

# Loop over pairs
for (n_pair in 1:ncol(pairs)) {
  
  # Declare Time Series Objects
  pair <- pairs[, n_pair]
  
  ts1 <- ts(data[,pair[1]], start = c(1,1), end = c(2265,1), frequency = 1) # end = 2265, as we are modelling until end of 2018
  ts2 <- ts(data[,pair[2]], start = c(1,1), end = c(2265,1), frequency = 1)
  
  # Bind into a system
  
  dset <- cbind(ts1,ts2)
  
  # Test lag length for pairs
  
  lagselect <- VARselect(dset, lag.max = 20, type = "trend")
  selected_lag[n_pair] <- lagselect$selection[1]
  if(selected_lag[n_pair] == 1){
    used_lag[n_pair] <- which.min(lagselect$criteria[1, -1]) + 1
  } else {
    used_lag[n_pair] <- selected_lag[n_pair]
  }
  
  print(pair)
  print(selected_lag[n_pair])
  
  fitted_var_model <- VAR(dset, p = used_lag[n_pair], type ="trend")
  
  # Misspecifications tests
  
  misspecification_autocorrelation_detected[n_pair] <- serial.test(fitted_var_model, lags.bg = used_lag[n_pair], type = "BG")$serial$p.value
  misspecification_arch_detected[n_pair] <- 1 # arch.test(fitted_var_model)$arch$p.value 
  misspecification_normality_detected[n_pair] <- 1 # normality.test(fitted_var_model)$jb.mul$JB$p.value
  

  if(min(misspecification_autocorrelation_detected[n_pair], 
         misspecification_arch_detected[n_pair],
         misspecification_normality_detected[n_pair])  < 0.05) {
    next
  } 
  # If the residuals suffer from autocorrelation discard and move on to the next model
  
  #Wild Bootstrap
  ctrace <- cointBootTest(dset, r = "sequence", p = used_lag[n_pair], model = 3, signif = 0.05, B = 399, boot_type = "WB", WB_dist = "rademacher")
  ca_jo_trace_test_stats[n_pair] <- ctrace$WB.pv[1] 
  no_cointegration_rejected[n_pair] <- (ca_jo_trace_test_stats[n_pair] < 0.05 )
}

pairs_with_cointegration <- pairs[,which(no_cointegration_rejected==1 )]

named_pairs_with_cointegration <- apply(pairs_with_cointegration, 2, function(r) colnames(data)[r])

named_pairs_with_cointegration

# Getting summary for VAR

var_with_cointegration <- list()

for (n_pair in 1:ncol(pairs_with_cointegration)) {
  pair <- pairs_with_cointegration[, n_pair]
  

  ts1 <- ts(data[,pair[1]], start = c(1,1), end = c(2265,1), frequency = 1)
  ts2 <- ts(data[,pair[2]], start = c(1,1), end = c(2265,1), frequency = 1)
  
  # Bind into a system
  dset <- cbind(ts1,ts2)
  
  lagselect <- VARselect(dset, lag.max = 20, type = "trend")
  selected_lag[n_pair] <- lagselect$selection[1]
  if(selected_lag[n_pair] == 1){
    used_lag[n_pair] <- which.min(lagselect$criteria[1, -1]) + 1
  } else {
    used_lag[n_pair] <- selected_lag[n_pair]
  }
  
  fitted_var_model <- VAR(dset, p = used_lag[n_pair], type ="trend")
  
  var_with_cointegration[[n_pair]] <- fitted_var_model
}


# Getting summary for Bootstrap Cointegration

var_with_cointegration_boot <- list()

for (n_pair in 1:ncol(pairs_with_cointegration)) {
  pair <- pairs_with_cointegration[, n_pair]
  
  
  ts1 <- ts(data[,pair[1]], start = c(1,1), end = c(2265,1), frequency = 1)
  ts2 <- ts(data[,pair[2]], start = c(1,1), end = c(2265,1), frequency = 1)
  
  # Bind into a system
  dset <- cbind(ts1,ts2)
  
  ctrace <- cointBootTest(dset, r = "sequence", p = used_lag[n_pair], model = 3, signif = 0.05, B = 399, boot_type = "WB", WB_dist = "rademacher")

  
  var_with_cointegration_boot[[n_pair]] <- ctrace
}


summary_file <- "summary.txt"
write("Summaries for S&P 500 modelling:", file=summary_file)
# Summaries 
count_number_of_pairs_after_beta_adjust <- 0
equation_of_selected_pairs <- numeric(0)
z_scores_selected <- list()
roots_selected <- list()

for (n_pair in 1:ncol(pairs_with_cointegration)) {
  pair <- pairs_with_cointegration[, n_pair]
  
  name_1 <- named_pairs_with_cointegration[1, n_pair]
  name_2 <- named_pairs_with_cointegration[2, n_pair]
  
  beta <- var_with_cointegration_boot[[n_pair]]$beta[2,1]
  
  rho <- var_with_cointegration_boot[[n_pair]]$rho[1,1]
  
  alpha <- var_with_cointegration_boot[[n_pair]]$alpha[1,1]
  
  roots <- Rfast::nth(Re(sort(var_with_cointegration_boot[[n_pair]]$companion_eigen[[2]])),2, descending = T)
  
  # Removing beta-values not between0 -0.5 and 2. 
  
  if (!between(beta, -0.5, 2)) {
    next
  }
  # Speed of convergence chosen due to Half-Life
  if ((roots > 0.975)  ) {
    next
  }

  cointegrated_residual <- data[,pair[1]] - data[,pair[2]] - rho
  
  cointegrated_residual_mean <- mean(cointegrated_residual[1:2265])
  cointegrated_residual_sd <- sd(cointegrated_residual[1:2265])
  
  z_scores_selected[[n_pair]] <- (cointegrated_residual[2265:2513] - cointegrated_residual_mean)/cointegrated_residual_sd
  roots_selected[[n_pair]] <- Rfast::nth(Re(sort(var_with_cointegration_boot[[n_pair]]$companion_eigen[[2]])),2, descending = T)
  
  count_number_of_pairs_after_beta_adjust <- count_number_of_pairs_after_beta_adjust+1
  
  eq <- paste(n_pair, ":", name_1, "=", round(beta, 3), "*", name_2, "+", round(rho, 3))
  
  equation_of_selected_pairs[count_number_of_pairs_after_beta_adjust] <- eq
  
  cat(eq, file=summary_file, append=TRUE, sep = "\n")
  
}

cat(paste("Number of pairs:", count_number_of_pairs_after_beta_adjust), file = summary_file, append = TRUE, sep ="\n")

z_score_matrix <- matrix(unlist(z_scores_selected), byrow=TRUE, nrow=count_number_of_pairs_after_beta_adjust)


for(i in 1:nrow(z_score_matrix)) {
  if (i!=5) {
   next
  }
  plot(2265:2513, z_score_matrix[i,], ylim=c(-5, 5), type="l", xlab="Time",  ylab="z-score", main=equation_of_selected_pairs[i])
  abline(h=0, lty=2, col=4)
  abline(h=1.25, col=2)
  abline(h=-1.25, col=2)
  readline(prompt="Press Enter for next plot.")
}



# Trading
write.table(z_score_matrix, "z_scores.txt", row.names=FALSE, col.names=FALSE)
z_loaded <- unname(as.matrix(read.table("z_scores.txt")))
