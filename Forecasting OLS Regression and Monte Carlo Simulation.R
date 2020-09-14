# Empirical Asset Pricing - Assignment1
#### Load data and packages ####
setwd('/Users/regina/Documents/UWO/Asset Pricing')
EAP1 = read.table("A1_data.csv", header = TRUE, sep = ",")
EAP1 = subset(EAP1, EAP1$caldt <= 20041231)

library(fBasics)
library(dplyr)
library(timeDate)
library(tseries)
library(stats)
library(dynlm)


#### Forecasting Regressions ####
Ret = EAP1$Return
ExRet = EAP1$Returnx
DP = EAP1$DP
Div_growth = EAP1$D_growth

# OLS Regression of Real Return on D/P Ratio
n = nrow(EAP1)
reg1 = lm(Ret[2:n] ~ DP[1:n-1])
summary(reg1)

b1 = summary(reg1)$coefficients[2,1]
tval1 = summary(reg1)$coefficients[2,3]
Rsq1 = summary(reg1)$r.squared
sd1 = sd(reg1$fitted.values)

# OLS Regression of Excess Return on D/P Ratio
reg2 = lm(ExRet[2:n] ~ DP[1:n-1])
b2 = summary(reg2)$coefficients[2,1]
tval2 = summary(reg2)$coefficients[2,3]
Rsq2 = summary(reg2)$r.squared
sd2 = sd(reg2$fitted.values)

# OLS Regressions of real dividend growth on D/P ratio
reg3 = lm(Div_growth[2:n] ~ DP[2:n])
b3 = summary(reg3)$coefficients[2,1]
tval3 = summary(reg3)$coefficients[2,3]
Rsq3 = summary(reg3)$r.squared
sd3 = sd(reg3$fitted.values)

# Results table
Results = matrix(c(b1, tval1, Rsq1, sd1,
                   b2, tval2, Rsq2, sd2,
                   b3, tval3, Rsq3, sd3), nrow = 3, ncol = 4, byrow = TRUE)
dimnames(Results) = list(c("Real Return on D/P Ratio", "Excess Return on D/P Ratio", "Real Dividend Growth on D/P ratio"), 
                         c("b", "t(b)", "R-squared", "SD"))
Results


#### Monte Carlo ####
# Parameters
phi = 0.941
rho = 0.9638
vec_sd = c(0.196, 0.14, 0.153)
nSim = 50000 + 1

# function - Simulated data
Sim_data = function(nSim = 50000 + 1, ...){
  
  shocks_dp = rnorm(nSim, sd = vec_sd[3])
  shocks_d = rnorm(nSim-1, sd = vec_sd[2])
  vec_dp = numeric(nSim)
  
  vec_dp[1] <- rnorm(1, sd = vec_sd[3]^2/(1-phi^2))
  for (i in 2:nSim) {
    vec_dp[i] = phi*vec_dp[i-1] + shocks_dp[i]
  }
  vec_d = (phi*rho - 1)*vec_dp[1:nSim-1] + shocks_d  ## br = 0
  vec_r = shocks_d - rho*shocks_dp[2:nSim]
  
  return(data.frame(dp = vec_dp[2:nSim],
                    div_growth = vec_d,
                    ret = vec_r))
}

# function - Regression with simulated data
MonteCarloReg = function(nSim = 50000+1, ...){
  
  Stats = matrix(nrow = nSim, ncol = 6)
  
  for (i in 2:nSim) {
    dat = Sim_data(...)
    list_reg = list(lm(dp ~ lag(dp), data = dat),
                    lm(div_growth ~ lag(dp), data = dat),
                    lm(ret ~ lag(dp), data = dat))
    Stats[i,1:3] = unlist(lapply(list_reg, FUN = function(reg) summary(reg)$coef[2,1]))
    Stats[i,4:6] = unlist(lapply(list_reg, FUN = function(reg) summary(reg)$coef[2,3]))
  }
  
  Stats = na.omit(as.data.frame(Stats))
  reg_names = c("dp", "div_growth", "ret")
  names(Stats) = c(paste0("Coef_", reg_names), paste0("tval_", reg_names))
  
  return(Stats)
}

# Regression results
Sim_results = MonteCarloReg(nSim = 50000+1)

# Plot
library(ggplot2)
par(mfrow=c(2,2))

# plot maximum 5000 points
ggplot(Sim_results[1:min(nSim,5000),], aes(x = Coef_ret, y = Coef_div_growth)) + 
  geom_point() +
  geom_vline(xintercept = 0.097) + 
  geom_hline(yintercept = 0.008) +
  ggtitle("Coefficients, phi = 0.94") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(expression(b[r])) +
  ylab(expression(b[d]))
ggplot(Sim_results[1:min(nSim,5000),], aes(x = tval_ret, y = tval_div_growth)) + 
  geom_point() +
  geom_vline(xintercept = 1.92) + 
  geom_hline(yintercept = 0.18) +
  ggtitle("t-stats, phi = 0.94") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(expression(t(b[r]))) +
  ylab(expression(t(b[d])))





