# Statistical Arbitrage in the U.S. Equities Market - Pairs Trading Standard \& Poor's 500
# Author: Jens Pedersen


library(urca)
library(vars)
library(VARtests)
library(tsDyn)
library(forecast)
library(tidyverse)
library(tseries)
library(Rfast)
setwd("C:/Users/jop.msc/Google Drev/KU/BSc i Økonomi/6. Semester/Bachelorprojekt/Datasæt")

# Loading dataset
data <- read.csv("SP500 Log.csv", sep = ";", stringsAsFactors = F, dec = ",")
dates <- data[, 1]
data <- data[, -1]
rownames(data) <- dates
# Omitting last three rows due to data corruption
data <- data[1:2513,]

# Finding all combinations of pairs
pairs <- combn(1:(dim(data)[2]), 2)

# 
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
  
  # Misspecifications tests (only testing for serial correlation)
  misspecification_autocorrelation_detected[n_pair] <- serial.test(fitted_var_model, lags.bg = used_lag[n_pair], type = "BG")$serial$p.value
  misspecification_arch_detected[n_pair] <- 1 # arch.test(fitted_var_model)$arch$p.value 
  misspecification_normality_detected[n_pair] <- 1 # normality.test(fitted_var_model)$jb.mul$JB$p.value
  
  # If the residuals suffer from serial correlation discard and move on to the next model

  if(min(misspecification_autocorrelation_detected[n_pair], 
         misspecification_arch_detected[n_pair],
         misspecification_normality_detected[n_pair])  < 0.05) {
    next
  } 
  
  
  # Bootstrap testing the co-integration rank
  ctrace <- cointBootTest(dset, r = "sequence", p = used_lag[n_pair], model = 3, signif = 0.05, B = 399, boot_type = "WB", WB_dist = "rademacher")
  ca_jo_trace_test_stats[n_pair] <- ctrace$WB.pv[1] 
  no_cointegration_rejected[n_pair] <- (ca_jo_trace_test_stats[n_pair] < 0.05 )
}

pairs_with_cointegration <- pairs[,which(no_cointegration_rejected==1 )]

named_pairs_with_cointegration <- apply(pairs_with_cointegration, 2, function(r) colnames(data)[r])

named_pairs_with_cointegration

# Getting summaries for CVAR models

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

# Getting summaries for the bootstrapped models

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

# Creating a text file for all pairs that hold certain criterias defined beneath

summary_file <- "summary.txt"
write("Summaries for S&P 500 modelling:", file=summary_file)

count_number_of_pairs_after_beta_adjust <- 0
equation_of_selected_pairs <- numeric(0)
z_scores_selected <- list()
roots_selected <- list()
pairs_selected <- numeric(0)

for (n_pair in 1:ncol(pairs_with_cointegration)) {
  pair <- pairs_with_cointegration[, n_pair]
  
  name_1 <- named_pairs_with_cointegration[1, n_pair]
  name_2 <- named_pairs_with_cointegration[2, n_pair]
  
  beta <- var_with_cointegration_boot[[n_pair]]$beta[2,1]
  
  rho <- var_with_cointegration_boot[[n_pair]]$rho[1,1]
  
  alpha <- var_with_cointegration_boot[[n_pair]]$alpha[1,1]
  
  roots <- Rfast::nth(Re(sort(var_with_cointegration_boot[[n_pair]]$companion_eigen[[2]])),2, descending = T)
  
  # Removing beta-values not between -0.5 and 2. 
  
  if (!between(beta, -0.5, 2)) {
    next
  }
  # Finding the second largest unrestricted eigenvalue in the companion matrix. (Speed of convergence to equilibrium)
  if ((roots > 0.97)  ) {
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
  pairs_selected[count_number_of_pairs_after_beta_adjust] <- n_pair
  
  cat(eq, file=summary_file, append=TRUE, sep = "\n")
  
}

cat(paste("Number of pairs:", count_number_of_pairs_after_beta_adjust), file = summary_file, append = TRUE, sep ="\n")

z_score_matrix <- matrix(unlist(z_scores_selected), byrow=TRUE, nrow=count_number_of_pairs_after_beta_adjust)


# Plotting z-scores for each pair

for(i in 1:nrow(z_score_matrix)) {
  par(bg=(col=rgb(254,230,206, maxColorValue = 255)))
  par("usr")
  plot(2265:2513, z_score_matrix[i,], ylim=c(-4, 4), type="l", xlab="Time",  ylab="z-score", main=equation_of_selected_pairs[i], col="red", xaxt ='t, t+1, ... , t+n')
  abline(h=1.25, col=rgb(230,85,13, maxColorValue = 255))
  abline(h=-1.25, col=rgb(230,85,13, maxColorValue = 255))
  abline(h=0.5, col=rgb(253,174,107, maxColorValue = 255), lty =2)
  abline(h=-0.5, col=rgb(253,174,107, maxColorValue = 255), lty =2)
  abline(v = 2305, lty=3, col=rgb(127,39,4, maxColorValue = 255))
  abline(v = 2335, lty=3, col=rgb(217,72,1, maxColorValue = 255))
  abline(v = 2360, lty=3, col=rgb(127,39,4, maxColorValue = 255))
  abline(v = 2440, lty=3, col=rgb(217,72,1, maxColorValue = 255))
  abline(v = 2492, lty=3, col=rgb(127,39,4, maxColorValue = 255))
  readline(prompt="Press Enter for next plot.")
}

par(bg=(col=rgb(254,230,206, maxColorValue = 255)))
par("usr")
plot(2265:2513, z_score_matrix[11,], ylim=c(-4, 4), type="l", xlab="Time",  ylab="z-score", main=equation_of_selected_pairs[n_pair] , col="red", xaxt ='t, t+1, ... , t+n')
abline(h=1.25, col=rgb(230,85,13, maxColorValue = 255))
abline(h=-1.25, col=rgb(230,85,13, maxColorValue = 255))
abline(h=0.5, col=rgb(253,174,107, maxColorValue = 255), lty =2)
abline(h=-0.5, col=rgb(253,174,107, maxColorValue = 255), lty =2)
abline(v = 2305, lty=3, col=rgb(127,39,4, maxColorValue = 255))
abline(v = 2335, lty=3, col=rgb(217,72,1, maxColorValue = 255))
abline(v = 2360, lty=3, col=rgb(127,39,4, maxColorValue = 255))
abline(v = 2440, lty=3, col=rgb(217,72,1, maxColorValue = 255))
abline(v = 2492, lty=3, col=rgb(127,39,4, maxColorValue = 255))


# Trading strategy
write.table(z_score_matrix, "z_scores.txt", row.names=FALSE, col.names=FALSE)
z_loaded <- unname(as.matrix(read.table("z_scores.txt")))

balance_func <- function(z_scores, prices_1, prices_2, beta) {
  if(length(prices_1) != length(z_scores)) {
    stop()
  }
  if(length(prices_2) != length(z_scores)) {
    stop()
  }
  t <- length(z_scores)
  balance <- 0
  stock_status <- 0 # 0 hverken købt eller solgt, 1 købt stock_1, -1 solgt stock_1
  for(i in 1:t) {
    current_z <- z_scores[i]
    
    if((current_z < -1.25) & (stock_status == 0)) {
      balance <- balance - ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] + ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
      stock_status <- 1
    }
    
    if((current_z > 1.25) & (stock_status == 0)) {
      balance <- balance + ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] - ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
      stock_status <- -1
    }
    
    if((abs(z_scores[i]) > 3.25) & (stock_status != 0)) {
      if(stock_status == 1) {
        balance <- balance + ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] - ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
        stock_status <- 0
      }
      if(stock_status == -1) {
        balance <- balance - ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] + ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
        stock_status <- 0
      }
    }
    
    if((current_z > -0.5) & (stock_status == 1)) {
      balance <- balance + ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] - ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
      stock_status <- 0
    }
    
    if((current_z < 0.5) & (stock_status == -1)) {
      balance <- balance - ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] + ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
      stock_status <- 0
    }
  }
  
  if((stock_status == 1)) {
    balance <- balance + ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] - ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
    stock_status <- 0
  }
  
  if((stock_status == -1)) {
    balance <- balance - ((prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_1[i] + ((beta*prices_1[i]/(prices_1[i]+beta*prices_1[i])))*prices_2[i]
    stock_status <- 0
  }
  
  return(balance)
}

balances <- numeric(length(pairs_selected))
returns <- numeric(length(pairs_selected))
beta <- numeric(length(pairs_selected))
price_3 <- numeric(length(pairs_selected))
price_4 <- numeric(length(pairs_selected))
for(i in 1:length(pairs_selected) ) {
  z_scores <- z_score_matrix[i, ]
  pair <- pairs_with_cointegration[, pairs_selected[i]]
  prices_1 <- data[2265:2513, pair[1]]
  prices_2 <- data[2265:2513, pair[2]]
  prices_3[i] <- data[2265, pair[1]] 
  prices_4[i] <- data[2265, pair[2]] 
  beta[i] <- var_with_cointegration_boot[[pairs_selected[i]]]$beta[2,1]
  balances[i] <- balance_func(z_scores, prices_1, prices_2, beta)
  returns[i] <- balances[i]/(prices_3[i]+beta*prices_4[i])
}



