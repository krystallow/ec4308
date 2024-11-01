rm(list=ls())

#Auxiliary function to compute root MSE (same as MSE before, but with square root):
RMSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(sqrt(mean((truth - pred)^2)))
} 

#install.packages("githubinstall") #this package is needed to install packages from GitHub (a popular code repository)
library(githubinstall)

#Here, some of you may encounter the following issue:
#You may see an error like this when trying to load the "githubinstall" library:
#Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) :  namespace 'cli' 3.x.x is already loaded, but >= 3.x.x is required

#If you see this, the following should fix it (reinstall the "cli" package):
#remove.packages("cli")
#install.packages("cli")

#install Medeiros et al's package, you will be prompted to say "Yes" to confirm the name of the installed package:

#githubinstall("HDeconometrics")

#You will see something like:
#Suggestion:
#  - gabrielrvsc/HDeconometrics  Set of R functions for high-dimensional econometrics
#Do you want to install the package (Y/n)?  
#Type "Y" in the command line and press Enter

library(HDeconometrics)
#install.packages("sandwich")
library(sandwich) #library to estimate variance for DM test regression using NeweyWest()
library(randomForest)
library(hdm)
library(readr)
library(xts)
library(dplyr)
library(tibble)
###########################
###########################
###########################
#Set working directory and load data:
#setwd("D:/data") #set your working directory here - or load data directly via RStudio interface

# Move 'date' column to rownames and remove it from the dataframe
Final_Transformed_Dataset <- read_csv("Final_Transformed_Dataset.csv", 
                                      col_types = cols(sasdate = col_date(format = "%m/%d/%Y")))
df <- Final_Transformed_Dataset %>%
  column_to_rownames(var = "sasdate")

### Some preliminary data manipulation
Y=df #this is the matrix of all variables, Y and X. It is in a format of a matrix with added column names.
dum=rep(0,nrow(Y)) # 1st step in creating the 2020 outlier dummy - a vector of zeros of conformable length to rest of the data
dum[which.min(Y[,1])]=1 #create dummy for outlier in 2020 - corresponds to the minimum of the first column, which is INDPRO
Y=cbind(Y,dum=dum) #add dummy to data matrix

yy=Y[,1] #get the y variable - INDPRO

nprev=230 #number of out-of-sample observations (test window - 30 % of the observations)

oosy=tail(yy,nprev) #auxiliary:get the out-of-sample true values (last 230 obs. using tail())

#create lags
rwtemp=embed(yy,13)
# this creates a matrix where each row contains the current observation of yy along with the 12 preceding observations

#Simple Random Walk forecast:
rw1c=tail(rwtemp[,2],nprev)
rw3c=tail(rwtemp[,4],nprev) 
rw6c=tail(rwtemp[,7],nprev)
rw12c=tail(rwtemp[,13],nprev)
# rwtemp[,13] refers to the fourth column of the rwtemp matrix, which contains the values lagged by 12 periods
# 12-step ahead random walk forecast, where the value at time t is predicted to be the same as the value at time t−12

#Collect RMSE's for randomw walk:
rw.rmse1=RMSE(oosy,rw1c)
rw.rmse3=RMSE(oosy,rw3c)
rw.rmse6=RMSE(oosy,rw6c)
rw.rmse12=RMSE(oosy,rw12c)

######################################
### Getting optimal lags using BIC ###
######################################
# Load necessary library
#install.packages("flexmix")
library(flexmix)

optimal_lags_bic <- function(Y, max_lag = 10) {
  
  bic_values <- numeric(max_lag)  # Initialize vector to store BIC values
  
  # Loop through 1 to max_lag and fit AR models
  for (lag in 1:max_lag) {
    
    # Fit AR model using arima (with no differencing and no MA terms)
    ar_model <- arima(Y, order = c(lag, 0, 0))  # order = (AR lags, differencing, MA lags)
    
    # Compute BIC using the BIC() function from flexmix
    bic_values[lag] <- BIC(ar_model)
  }
  
  # Find the lag that gives the minimum BIC
  optimal_lag <- which.min(bic_values)
  
  # Return the optimal lag and the BIC values for all lags
  return(list("optimal_lag" = optimal_lag, "bic_values" = bic_values))
}

result <- optimal_lags_bic(yy, max_lag = 12)
optimal_lag <- result$optimal_lag  # Optimal number of lags
print(optimal_lag)
print(result$bic_values)   # BIC values for all lags

######################################
#Benchmark 2: AR(p) forecast
######################################

#Perform AR(p) forecasts using rolling window.
#See the file func-ar.R for the forecast construction details there

#Add the functions  in func-ar.R (must be in your working directory)
#Or simply open up func-ar.R and execute the function commands there
source("func-ar.R")
par(mfrow = c(2, 2))
bar1c=ar.rolling.window(Y,nprev,1,1,type="bic") #1-step AR forecast
title(main = 'AR(1) 1-Step Forecast')
bar3c=ar.rolling.window(Y,nprev,1,3,type="bic") #3-step AR forecast
title(main = 'AR(1) 3-Step Forecast')
bar6c=ar.rolling.window(Y,nprev,1,6,type="bic") #6-step AR forecast
title(main = 'AR(1) 6-Step Forecast')
bar12c=ar.rolling.window(Y,nprev,1,12,type="bic") #12-step AR forecast
title(main = 'AR(1) 12-Step Forecast')

#Benchmark forecast graphics:

#Plot benchmark coefficients
#Here I use plot.ts(), which plots time series objects that are dated

#First, I make a coefficient matrix a time-series object by using ts() function
#Main reason: I want dates on the X-axis

#Syntax: ts(object, start=startdate, end=enddate, freq=frequency (periods per year))
arcoef.ts=ts(bar1c$coef, start=c(2010,1), end=c(2024,7), freq=12)
colnames(arcoef.ts)=c("Constant","1st Lag","2nd Lag") #name columns to distinguish plots

#Plot all the coefficients over time (plot.ts() same as plot, but for tme series objects):
#quartz()
plot.ts(arcoef.ts, main="AR regression coefficients", cex.axis=1.5)

#Similarly, I create ts objects out of 1-step to 12-step benchmark forecasts
bench1.ts=ts(cbind(rw1c,bar1c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench1.ts)=c("RW","AR(1)","True Value")

bench3.ts=ts(cbind(rw3c,bar3c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench3.ts)=c("RW","AR(1)","True Value")

bench6.ts=ts(cbind(rw6c,bar6c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench6.ts)=c("RW","AR(1)","True Value")

bench12.ts=ts(cbind(rw12c,bar12c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench12.ts)=c("RW","AR(1)","True Value")

#Plot 1-step forecasts:
par(mfrow = c(2, 2))
#quartz()
# plot.ts(bench1.ts[,1], main="1-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="blue", ylab="INDRPO")
# points(bench1.ts[,2], type="l", col="red",lwd=1.8) # AR(1) Model
# points(bench1.ts[,3], type="l", col="black",lwd=2) #True Value
# legend("bottomleft",legend=c("RW","AR(1)","INDPRO"))

plot.ts(bench1.ts[,3], main="1-step AR Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="black", ylab="INDRPO") #True Value
points(bench1.ts[,2], type="l", col="red",lwd=1.8) # AR(1) Model
legend("bottomright",legend=c("AR(1)","INDPRO"),col = c("red", "black"),lwd = 1.8,lty = 1,cex = 0.4)

#Plot 3-step forecasts:
plot.ts(bench3.ts[,3], main="3-step AR Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="black", ylab="INDRPO") #True Value
points(bench3.ts[,2], type="l", col="red",lwd=1.8) # AR(1) Model
legend("bottomright",legend=c("AR(1)","INDPRO"),col = c("red", "black"),lwd = 1.8,lty = 1,cex = 0.4)

#Plot 6-step forecasts:
plot.ts(bench6.ts[,3], main="6-step AR Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="black", ylab="INDRPO") #True Value
points(bench6.ts[,2], type="l", col="red",lwd=1.8) # AR(1) Model
legend("bottomright",legend=c("AR(1)","INDPRO"),col = c("red", "black"),lwd = 1.8,lty = 1,cex = 0.4)

#Plot 12-step forecasts:

#quartz()
# plot.ts(bench12.ts[,1], main="12-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="blue", ylab="INDPRO")
# points(bench12.ts[,2], type="l", col="red",lwd=1.8)
# points(bench12.ts[,3], type="l", col="black",lwd=2)
# legend("bottomleft",legend=c("RW","AR(1)","INDPRO"))
plot.ts(bench12.ts[,3], main="12-step AR Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="black", ylab="INDRPO") #True Value
points(bench12.ts[,2], type="l", col="red",lwd=1.8) # AR(1) Model
legend("bottomright",legend=c("AR(1)","INDPRO"),col = c("red", "black"),lwd = 1.8,lty = 1,cex = 0.4)

#AR forecasts RMSE:

ar.rmse1=bar1c$errors[1]
ar.rmse3=bar3c$errors[1]
ar.rmse6=bar6c$errors[1]
ar.rmse12=bar12c$errors[1]

# For DM test:
predicted_values = bar12c$pred