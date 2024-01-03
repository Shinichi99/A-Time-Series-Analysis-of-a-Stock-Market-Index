rm(list=ls())

######################
# Univariate Analysis#
######################

# Import data
dataset <- read.table(file=file.choose(),header=TRUE, sep=",")
sandp500 <- ts(dataset$sandp500, frequency = 12, start = c(1980, 1) )

#Plot time series in level
par(mfrow=c(1,1))
ts.plot(sandp500, main = "S&P 500 index", xlab = "Time", ylab = "S&P 500")
#transform time series
log_sandp500 <- log(sandp500)
dlog_sandp500 <- diff(log_sandp500)
#Plot log-transformed time series in first differences
ts.plot(dlog_sandp500, main = "Log-transformed S&P 500 in First Differences", xlab = "Time", ylab = "Log-transformed S&P 500 index in First Differences")
monthplot(dlog_sandp500, main = "Monthly Plot of Log-transformed S&P 500 in First Differences", xlab = "Month", ylab = "Log-transformed S&P 500 in First Differences")

#unit-root
library(CADFtest)
max.lag<-round(sqrt(length(log_sandp500))) 
CADFtest(dlog_sandp500, type= "drift", criterion= "BIC", max.lag.y=max.lag)
par(mfrow=c(2,1))
acf(dlog_sandp500, main = "Log-transformed S&P 500 in First Differences")
pacf(dlog_sandp500, main = "Log-transformed S&P 500 in First Differences")
Box.test(dlog_sandp500, lag = max.lag, type = "Ljung-Box")

### Model Selection ###
#function for MAE
Forecast <- function(y, h, order, seasonal) {
  error_h <- c()
  S <- round(0.75 * length(y))
  for (i in S:(length(y) - h)) {
      mymodel_sub <- arima(y[1:i], order = order, seasonal = seasonal)
      predict_h <- predict(mymodel_sub, n.ahead = h)$pred[h]
      error_h <- c(error_h, y[i + h] - predict_h)
  }
  return(error_h)
}

#MA(1)
fit_ma1 <- arima(log_sandp500, order = c(0,1,1) ,seasonal=c(0,0,0))
fit_ma1
par(mfrow=c(1,1))
ts.plot(fit_ma1$residuals)
acf(fit_ma1$residuals)
Box.test(fit_ma1$residuals, lag = max.lag, type = "Ljung-Box")
AIC(fit_ma1,k=log(length(dlog_sandp500)))
error_ma1 = Forecast(log_sandp500,1,c(0,1,1),c(0,0,0))
mean(abs(error_ma1))
error4_ma1 = Forecast(log_sandp500,4,c(0,1,1),c(0,0,0))
mean(abs(error4_ma1))
acf(fit_ma1$residuals^2)

#AR(1)
fit_ar1 <- arima(log_sandp500, order = c(1,1,0), seasonal =c(0,0,0))
fit_ar1
par(mfrow=c(1,1))
ts.plot(fit_ar1$residuals)
acf(fit_ar1$residuals)
Box.test(fit_ar1$residuals, lag = max.lag, type = "Ljung-Box")
AIC(fit_ar1,k=log(length(dlog_sandp500)))
error_ar1 = Forecast(log_sandp500,1,c(1,1,0),c(0,0,0))
mean(abs(error_ar1))
error4_ar1 =Forecast(log_sandp500,4,c(1,1,0),c(0,0,0))
mean(abs(error4_ar1))
    
#ARIMA(1,1,1) 
fit_arima1 <- arima(log_sandp500, order = c(1,1,1), seasonal =c(0,0,0))
fit_arima1
par(mfrow=c(1,1))
ts.plot(fit_arima1$residuals)
acf(fit_arima1$residuals)
Box.test(fit_arima1$residuals, lag = max.lag, type = "Ljung-Box")
AIC(fit_arima1,k=log(length(dlog_sandp500)))
error_arima1 = Forecast(log_sandp500,1,c(1,1,1),c(0,0,0))
mean(abs(error_arima1))
error4_arima1 =Forecast(log_sandp500,4,c(1,1,1),c(0,0,0))
mean(abs(error4_arima1))
acf(fit_arima1$residuals^2)

#GARCH
library(fGarch)
fit_garch<-garchFit(~arma(0,1)+garch(1,1),log_sandp500)
summary(fit_garch)
fit_garch<-garchFit(~arma(0,1)+garch(1,1),cond.dist="QMLE",log_sandp500)
summary(fit_garch)
fit_garch<-garchFit(~arma(1,0)+garch(1,1),log_sandp500)
summary(fit_garch)
fit_garch<-garchFit(~arma(1,0)+garch(1,1),cond.dist="QMLE",log_sandp500)
summary(fit_garch)
fit_garch<-garchFit(~arma(1,1)+garch(1,1),log_sandp500)
summary(fit_garch)
fit_garch<-garchFit(~arma(1,1)+garch(1,1),cond.dist="QMLE",log_sandp500)
summary(fit_garch)


#Diebold-Mariano test
library(forecast)
dm.test(error4_ma1,error4_ar1,h=4,power=1) #For absolute value loss, use power=1.
dm.test(error4_ma1,error4_arima1,h=4,power=1) 
dm.test(error4_ar1,error4_arima1,h=4,power=1) 



######################
# Multivariate #
######################

# Import treasury data
treasury <- ts(dataset$treasury, frequency = 12, start = c(1980, 1) )
# Plot treasury data
ts.plot(treasury, main = "US 10-year Treasury Rate", xlab = "Time", ylab = "US 10-year Treasury Rate")
log_treasury <- log(treasury)
dlog_treasury <- diff(log_treasury)
ts.plot(dlog_treasury, main = "Log-transformed US 10-year Treasury Rate in First Differences", xlab = "Time", ylab = "Log-transformed US 10-year Treasury Rate in First Differences")
#unit-root test
max.lag<-round(sqrt(length(log_treasury))) 
CADFtest(dlog_treasury, type= "drift", criterion= "BIC", max.lag.y=max.lag)


### VAR model ###
library(vars)
dlogdata<-data.frame(dlog_sandp500, dlog_treasury)
names(dlogdata)<-c("dlog_sandp500","dlog_treasury")
VARselect(dlogdata,lag.max=10,type="const")
fit_var1<-VAR(dlogdata,type="const",p=1)
summary(fit_var1)
var1_residuals<-resid(fit_var1)
par(mfrow=c(2,2))
acf(var1_residuals[,1])
acf(var1_residuals[,2])
ccf(var1_residuals[,1],var1_residuals[,2])

###  Impulse response ### 
par(mfrow=c(1,1))
irf_var<-irf(fit_var1,ortho=F,boot=T)
plot(irf_var)
#Enlarge Figure
irf_var <- irf(fit_var1, impulse = "dlog_treasury", response = "dlog_sandp500", ortho = F,boot=T)
plot(irf_var)


### VECM model ###
data<-data.frame(log_sandp500,log_treasury)
VARselect(data,lag.max=10,type="const")
library(urca)
trace_test<-ca.jo(data,type="trace",K=2,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(data,type="eigen",K=2,ecdet="const",spec="transitory")
summary(maxeigen_test)
fit_vecm1<-cajorls(trace_test,r=1)
fit_vecm1
fit_vecm2<-cajorls(maxeigen_test,r=1)
fit_vecm2

### Granger causality ###
lag <- 1
n <- length(dlog_sandp500)
dlog_sandp500.0 <- dlog_sandp500[(lag+1):n]
dlog_sandp500.1 <- dlog_sandp500[lag:(n-1)]
dlog_treasury.1 <- dlog_treasury[lag:(n-1)]
fit_dlm <- lm(dlog_sandp500.0 ~ dlog_sandp500.1+dlog_treasury.1)
acf(fit_dlm$residuals)
sqrt(length(dlog_sandp500))
Box.test(fit_dlm$residuals, lag = max.lag, type = "Ljung-Box")
fit_dlm_nox <- lm(dlog_sandp500.0 ~ dlog_sandp500.1)
anova(fit_dlm,fit_dlm_nox)


######################
# Forecasting #
######################

#Visualize forecasts
fit_var<-vec2var(trace_test,r=1)
myforecast<-predict(fit_var,n.ahead=12)
par(mfrow=c(1,1))
ts.plot(sandp500,xlim=c(2020,2025), ylim=c(0,6000))
log_sandp500_forecast<-ts(exp(myforecast$fcst$log_sandp500[,1]),frequency=12,start=c(2023,11))
log_sandp500_lower<-ts(exp(myforecast$fcst$log_sandp500[,2]),frequency=12,start=c(2023,11))
log_sandp500_upper<-ts(exp(myforecast$fcst$log_sandp500[,3]),frequency=12,start=c(2023,11))
lines(log_sandp500_forecast,col="red", lwd = 2)
lines(log_sandp500_lower,col="blue", lwd = 2)
lines(log_sandp500_upper,col="blue", lwd = 2)
title(main = "12-step-ahead Forecast of S&P500")
myforecast_ma1<-predict(fit_ma1,n.ahead=12)
expected_ma1<-exp(myforecast_ma1$pred)
lower_ma1<-exp(myforecast_ma1$pred-qnorm(0.975)*myforecast_ma1$se);
upper_ma1<-exp(myforecast_ma1$pred+qnorm(0.975)*myforecast_ma1$se);
cbind(lower_ma1,expected_ma1,upper_ma1)
lines(expected_ma1, col = "pink", lwd = 2, lty = "dashed")
lines(lower_ma1, col = "skyblue", lwd = 2, lty = "dashed")
lines(upper_ma1, col = "skyblue", lwd = 2, lty = "dashed")

legend("bottomleft", legend = c("VECM_Expected", "VECM_Bounds", "MA(1)_Expected", "MA(1)_Bounds"), 
       col = c("red", "blue", "pink", 'skyblue'), lty = 2)



fit_var<-vec2var(trace_test,r=1)
myforecast<-predict(fit_var,n.ahead=12)
par(mfrow=c(1,1))
ts.plot(treasury,xlim=c(2020,2025))
log_treasury_forecast<-ts(exp(myforecast$fcst$log_treasury[,1]),frequency=12,start=c(2023,11))
log_treasury_lower<-ts(exp(myforecast$fcst$log_treasury[,2]),frequency=12,start=c(2023,11))
log_treasury_upper<-ts(exp(myforecast$fcst$log_treasury[,3]),frequency=12,start=c(2023,11))
lines(log_treasury_forecast,col="red", lwd = 2)
lines(log_treasury_lower,col="blue", lwd = 2)
lines(log_treasury_upper,col="blue", lwd = 2)
title(main = "12-step-ahead Forecast of treasury")
legend("topleft", legend = c("VECM_Expected", "VECM_Bounds"), col = c("red", "blue"), lty = 1)

