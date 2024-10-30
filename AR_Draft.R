rm(list=ls())

#Auxiliary function to compute root MSE (same as MSE before, but with square root):
RMSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(sqrt(mean((truth - pred)^2)))
} #end function with a return(output) statement. Here we can go straight to return because the object of interest is a simple function of inputs


#Install the HDeconometrics package used in Medeiros et al. (2019) for convenient estimation
#of LASSO and ElNet using information criteria (basically uses glmnet, and selects on criterion)

install.packages("githubinstall") #this package is needed to install packages from GitHub (a popular code repository)
library(githubinstall)

#Here, some of you may encounter the following issue:
#You may see an error like this when trying to load the "githubinstall" library:
#Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) :  namespace 'cli' 3.x.x is already loaded, but >= 3.x.x is required

#If you see this, the following should fix it (reinstall the "cli" package):
#remove.packages("cli")
#install.packages("cli")


#install Medeiros et al's package, you will be prompted to say "Yes" to confirm the name of the installed package:

githubinstall("HDeconometrics")

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
# rwtemp[,4] refers to the fourth column of the rwtemp matrix, which contains the values lagged by 3 periods
# 3-step ahead random walk forecast, where the value at time t is predicted to be the same as the value at time t−3
rw6c=tail(rwtemp[,7],nprev)
rw12c=tail(rwtemp[,13],nprev)

#Collect RMSE's for randomw walk:
rw.rmse1=RMSE(oosy,rw1c)
rw.rmse3=RMSE(oosy,rw3c)
rw.rmse6=RMSE(oosy,rw6c)
rw.rmse12=RMSE(oosy,rw12c)

######################################
#Benchmark 2: AR(p) forecast
######################################

#Perform AR(p) forecasts using rolling window.
#See the file func-ar.R for the forecast construction details there

#Add the functions  in func-ar.R (must be in your working directory)
#Or simply open up func-ar.R and execute the function commands there
source("func-ar.R")


bar1c=ar.rolling.window(Y,nprev,1,1,type="fixed") #1-step AR forecast
bar3c=ar.rolling.window(Y,nprev,1,3,type="fixed") #3-step AR forecast
bar6c=ar.rolling.window(Y,nprev,1,6,type="fixed") #6-step AR forecast
bar12c=ar.rolling.window(Y,nprev,1,12,type="fixed") #12-step AR forecast
print(bar12c$pred)
#Benchmark forecast graphics:

#Plot benchmark coefficients
#Here I use plot.ts(), which plots time series objects that are dated

#First, I make a coefficient matrix a time-series object by using ts() function
#Main reason: I want dates on the X-axis

#Syntax: ts(object, start=startdate, end=enddate, freq=frequency (periods per year))
arcoef.ts=ts(bar1c$coef, start=c(2010,1), end=c(2024,7), freq=12)
colnames(arcoef.ts)=c("Constant","1st Lag","2nd Lag","3rd Lag","4th Lag") #name columns to distinguish plots

#Plot all the coefficients over time (plot.ts() same as plot, but for tme series objects):
quartz()
plot.ts(arcoef.ts, main="AR regression coefficients", cex.axis=1.5)

#Similarly, I create ts objects out of 1-step and 12-step benchmark forecasts
bench1.ts=ts(cbind(rw1c,bar1c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench1.ts)=c("RW","AR(4)","True Value")

bench12.ts=ts(cbind(rw12c,bar12c$pred,oosy), start=c(2010,1), end=c(2024,7), freq=12)
colnames(bench12.ts)=c("RW","AR(4)","True Value")

#Plot 1-step forecasts:

quartz()
plot.ts(bench1.ts[,1], main="1-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="blue", ylab="INDRPO")
points(bench1.ts[,2], type="l", col="red",lwd=1.8) # AR(4) Model
points(bench1.ts[,3], type="l", col="black",lwd=2) #True Value
legend("bottomleft",legend=c("RW","AR(4)","INDPRO"))

#Plot 12-step forecasts:

quartz()
plot.ts(bench12.ts[,1], main="12-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="blue", ylab="INDPRO")
points(bench12.ts[,2], type="l", col="red",lwd=1.8)
points(bench12.ts[,3], type="l", col="black",lwd=2)
legend("bottomleft",legend=c("RW","AR(4)","INDPRO"))

#AR forecasts RMSE:

ar.rmse1=bar1c$errors[1]
ar.rmse3=bar3c$errors[1]
ar.rmse6=bar6c$errors[1]
ar.rmse12=bar12c$errors[1]

