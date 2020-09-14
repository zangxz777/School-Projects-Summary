##### Load Data and Packages ####
setwd('/Users/regina/Documents/UWO/ECON9505 Econometrics')
max.print = 100000
ProjectData = read.table("data.csv", header = TRUE, sep = ",")

library(tseries)
library(forecast)
library(vars)
library(urca)


#### ARIMA Forecasting ####
# Test for Stationary - Not Stationary
WeeklyP = rev(ProjectData[,1])
adf.test(WeeklyP)

# First Difference Stationary
Fd_weeklyP = diff(WeeklyP, lag = 1, differences = 1)
adf.test(Fd_weeklyP)

# Compare and Choose optimal p&q order
auto.arima(WeeklyP, max.p = 10, max.q = 10, ic = "aic")
#ARIMA(1,1,0)
AIC(arima(WeeklyP, order = c(1,1,0)))
#ARIMA(0,1,1)
AIC(arima(WeeklyP, order = c(0,1,1)))
#ARIMA(0,1,0)
AIC(arima(WeeklyP, order = c(0,1,0)))

# In sample forecasting and Plot
summary(arima(WeeklyP, order = c(1,1,0)))
summary(arima(WeeklyP, order = c(0,1,1)))
summary(arima(WeeklyP, order = c(0,1,0)))

par(mfrow = c(3,1))
#ARIMA(1,1,0) forecast
n = length(WeeklyP)
F_P = c(WeeklyP[1],WeeklyP[2])
e = rnorm(n,0,sqrt(0.0002194))
for (i in 3:n){
  F_P[i] = 1.2449*F_P[i-1] - 0.2449*F_P[i-2] + e[i]
}
plot(F_P, type = 'l', 
          xlim = c(0,800),
          ylim = c(0,2), 
          col = 'red', 
          main = "ARIMA(1,1,0) in-sample forecast")
par(new = TRUE)
plot(WeeklyP, type = 'l', xlim = c(0,800), ylim = c(0,2))

#ARIMA(0,1,1) forecast
n = length(WeeklyP)
F_P = c(WeeklyP[1])
e = rnorm(n,0,sqrt(0.0002191))
for (i in 2:n){
  F_P[i] = F_P[i-1] + e[i]+ 0.2521*e[i-1]
}
plot(F_P, type = 'l',
          xlim = c(0,800),
          ylim = c(0,2), 
          col = 'red', 
          main = "ARIMA(0,1,1) in-sample forecast")
par(new = TRUE)
plot(WeeklyP, type = 'l', xlim = c(0,800), ylim = c(0,2))

#ARIMA(1,1,1) forecast
arima(WeeklyP, order = c(1,1,1))
n = length(WeeklyP)
F_P = c(WeeklyP[1], WeeklyP[2])
e = rnorm(n,0,sqrt(0.000219))
for (i in 3:n){
  F_P[i] = 1.072*F_P[i-1] - 0.072*F_P[i-2] + e[i]+ 0.1832*e[i-1]
}
plot(F_P, type = 'l',
          xlim = c(0,800),
          ylim = c(0,2), 
          col = 'red', 
          main = "ARIMA(1,1,1) in-sample forecast")
par(new = TRUE)
plot(WeeklyP, type = 'l', xlim = c(0,800), ylim = c(0,2))

# Out of sample forecasting and forecasting accuracy
# (0,1,1) vs (0,1,0)
# 50:50 split
n = length(WeeklyP)/2
historical = WeeklyP[1:n]
forecast1 = forecast(arima(historical, order = c(0,1,1)), h=n)$mean
forecast2 = forecast(arima(historical, order = c(0,1,0)), h=n)$mean
actual = WeeklyP[(n+1):length(WeeklyP)]

accuracy(forecast1, actual) #RMSE = 0.1409149, MAE = 0.1261144
accuracy(forecast2, actual) #RMSE = 0.143801, MAE = 0.1282423


#### VAR Forecasting ####
Euro = read.table("EURO.csv", header = FALSE, sep = ",")
US = read.table("US.csv", header = TRUE, sep = ",")
InRate = cbind(Euro, US)
Euro = InRate[,1]
US = InRate[,2]
InDiff = Euro - US
InDiff = InDiff[-1]

# Test for interest rate differential stationary
adf.test(InDiff)
# First difference stationary
Fd_InDiff = diff(InDiff)
adf.test(Fd_InDiff)
varmodel = cbind(Fd_weeklyP, Fd_InDiff)

# Johansen Test for Conintegration
JCT_in = Euro - US
JCT_wp = Fd_weeklyP
JCT_comb = cbind(JCT_in, JCT_wp)
Jcontest = ca.jo(JCT_comb, type = "trace", ecdet = "const", K = 5, spec = "longrun")
summary(Jcontest) #at least one cointegration

# Fit VAR model decide optimal lag length as per AIC
VARselect(varmodel, lag.max = 10, type = "both")$selection
ExRate_VAR = VAR(varmodel, type = "both", lag.max = 10, ic = "AIC")
summary(ExRate_VAR)

# Causality Test
causality(ExRate_VAR)

# Normality Test
normality.test(ExRate_VAR)

# Impulse Response Function
var_irf = irf(ExRate_VAR, impulse = "Fd_InDiff", response = "Fd_weeklyP")
plot(var_irf)

# Forecasting Accuracy 50:50 split
nhalf = nrow(varmodel)/2
historical_VAR = varmodel[1:nhalf]
actual_VAR = varmodel[(nhalf+1):nrow(varmodel)]

VARselect(historical_VAR, lag.max = 8, type = "both")$selection
VAR_fh = VAR(historical_VAR, type = "both", lag.max = 8, ic = "AIC")
summary(VAR_fh)
Fvar_fh <- predict(VAR_fh, n.ahead = nhalf, ci = 0.95)
summary(Fvar_fh)

accuracy(Fvar_fh$fcst$Fd_weeklyP[,1], actual_VAR[,1]) #RMSE = 0.01871266, MAE = 0.01464852






