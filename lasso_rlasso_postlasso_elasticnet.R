##################################
# Using only Selected Variables #
#################################
###########################################################
### Forcasting for LASSO, rLASSO, postLASSO, Elastic Net ###
############################################################

library(readxl)
library(sandwich) 
library(randomForest)
library(hdm)
library(glmnet)
library(githubinstall)
githubinstall("HDeconometrics")
library(HDeconometrics)

RMSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(sqrt(mean((truth - pred)^2)))
}

df <- read_excel("Final_Transformed_Dataset.xlsx", 
                                        col_types = c("date", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric"))

df <- df[-1, ]   
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("sasdate", "INDPRO", colnames(selected_variables))
df <- df[, selected_cols]


## Create dummy variable for 2020 Covid
df$sasdate <- as.Date(df$sasdate)
dum <- ifelse(format(df$sasdate, "%Y") == "2020", 1, 0)

Y <- cbind(df, dum = dum)
Y <- Y[,-1]
Y <- as.matrix(Y)

yy=Y[,1]


source("func-lasso.R")
source("func-rlasso.R")

nprev=230 #number of out-of-sample observations (test window ) # 30%

oosy=tail(yy,nprev) #auxiliary:get the out-of-sample true values 


############################################################################
#Penalized regression: LASSO forecasts (BIC, AIC, AICc, rlasso, post-LASSO)
############################################################################
alpha=1 #set alpha=1 for LASSO

#Run forecasts for LASSO (BIC)
lasso1c=lasso.rolling.window(Y,nprev,1,1,alpha,IC="bic")
lasso3c=lasso.rolling.window(Y,nprev,1,3,alpha,IC="bic")
lasso6c=lasso.rolling.window(Y,nprev,1,6,alpha,IC="bic")
lasso12c=lasso.rolling.window(Y,nprev,1,12,alpha,IC="bic")

#LASSO(BIC) RMSE's
lasso.rmse1=lasso1c$errors[1]
lasso.rmse3=lasso3c$errors[1]
lasso.rmse6=lasso6c$errors[1]
lasso.rmse12=lasso12c$errors[1]


# rLASSO #
#Note the final input are the coefficients from results of rLASSO:
rlasso1c=rlasso.rolling.window(Y,nprev,1,1,FALSE)
rlasso3c=rlasso.rolling.window(Y,nprev,1,3,FALSE)
rlasso6c=rlasso.rolling.window(Y,nprev,1,6,FALSE)
rlasso12c=rlasso.rolling.window(Y,nprev,1,12,FALSE)

# rLASSO RMSE's
rlasso.rmse1=rlasso1c$errors[1]
rlasso.rmse3=rlasso3c$errors[1]
rlasso.rmse6=rlasso6c$errors[1]
rlasso.rmse12=rlasso12c$errors[1]



# post-LASSO #
#Note the final input are the coefficients from results of rLASSO:
prlasso1c=rlasso.rolling.window(Y,nprev,1,1,TRUE)
prlasso3c=rlasso.rolling.window(Y,nprev,1,3,TRUE)
prlasso6c=rlasso.rolling.window(Y,nprev,1,6,TRUE)
prlasso12c=rlasso.rolling.window(Y,nprev,1,12,TRUE)

# post-LASSO RMSE's
prlasso.rmse1=rlasso1c$errors[1]
prlasso.rmse3=rlasso3c$errors[1]
prlasso.rmse6=rlasso6c$errors[1]
prlasso.rmse12=rlasso12c$errors[1]



########################################################################
#Penalized regression: Elastic Net (BIC, alpha=0.5)
########################################################################

alpha=0.5 #set alpha to 0.5 for elastic net

elasticnet1c=lasso.rolling.window(Y,nprev,1,1,alpha,IC="bic")
elasticnet3c=lasso.rolling.window(Y,nprev,1,3,alpha,IC="bic")
elasticnet6c=lasso.rolling.window(Y,nprev,1,6,alpha,IC="bic")
elasticnet12c=lasso.rolling.window(Y,nprev,1,12,alpha,IC="bic")

#See the RMSE
elnet.rmse1=elasticnet1c$errors[1]
elnet.rmse3=elasticnet3c$errors[1]
elnet.rmse6=elasticnet6c$errors[1]
elnet.rmse12=elasticnet12c$errors[1]



#######################################################
##Sparsity analysis over time plots for LASSO and ElNet
#######################################################

#Get nonzero coefficient numbers for different horizons (LASSO(BIC))
c1c=rowSums(prlasso1c$coef != 0)
c3c=rowSums(prlasso3c$coef != 0)
c6c=rowSums(prlasso6c$coef != 0)
c12c=rowSums(prlasso12c$coef != 0)

#Create a ts object for the plot
lcoef.ts=ts(cbind(c1c,c3c,c6c,c12c), start=c(1960,3), end=c(2024,07), freq=12)
colnames(lcoef.ts)=c("1-step","3-step","6-step","12-step")
#Plot numbers of nonzero coefficients across the test window
windows()
plot.ts(lcoef.ts, main="Sparsity Analysis for Post-LASSO",cex.axis=1.5)

#Get nonzero coefficient numbers for different horizons (Elastic Net)
ce1c=rowSums(elasticnet1c$coef != 0)
ce3c=rowSums(elasticnet3c$coef != 0)
ce6c=rowSums(elasticnet6c$coef != 0)
ce12c=rowSums(elasticnet12c$coef != 0)

#Create a respective ts object for the plot
elcoef.ts=ts(cbind(ce1c,ce3c,ce6c,ce12c), start=c(1960,3), end=c(2024,07), freq=12)
colnames(elcoef.ts)=c("1-step","3-step","6-step","12-step")
#Plot numbers of nonzero coefficients across the test window:
windows()
plot.ts(elcoef.ts, main="Sparsity Analysis for ElNet",cex.axis=1.5)










#######################
# Using all Variables #
#######################
###########################################################
### Forcasting for LASSO, rLASSO, postLASSO, Elastic Net ###
############################################################

library(readxl)
library(sandwich) 
library(randomForest)
library(hdm)
library(glmnet)
library(githubinstall)
githubinstall("HDeconometrics")
library(HDeconometrics)

RMSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(sqrt(mean((truth - pred)^2)))
}

df <- read_excel("Final_Transformed_Dataset.xlsx", 
                 col_types = c("date", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric"))

df <- df[-1, ]   

## Create dummy variable for 2020 Covid
df$sasdate <- as.Date(df$sasdate)
dum <- ifelse(format(df$sasdate, "%Y") == "2020", 1, 0)

Y <- cbind(df, dum = dum)
Y <- Y[,-1]
Y <- Y[apply(Y, 1, function(row) all(is.finite(row) & !is.na(row))), ]
Y <- as.matrix(Y)


yy=Y[,1]

source("func-lasso.R")
source("func-rlasso.R")

nprev=230 #number of out-of-sample observations (test window ) # 30%

oosy=tail(yy,nprev) #auxiliary:get the out-of-sample true values (last 180 obs. using tail())


############################################################################
#Penalized regression: LASSO forecasts (BIC, AIC, AICc, rlasso, post-LASSO)#
############################################################################
alpha=1 #set alpha=1 for LASSO


#Run forecasts for LASSO (BIC)
lasso1c=lasso.rolling.window.all(Y,nprev,1,1,alpha,IC="bic")
lasso3c=lasso.rolling.window.all(Y,nprev,1,3,alpha,IC="bic")
lasso6c=lasso.rolling.window.all(Y,nprev,1,6,alpha,IC="bic")
lasso12c=lasso.rolling.window.all(Y,nprev,1,12,alpha,IC="bic")

#LASSO(BIC) RMSE's
lasso.rmse1=lasso1c$errors[1]
lasso.rmse3=lasso3c$errors[1]
lasso.rmse6=lasso6c$errors[1]
lasso.rmse12=lasso12c$errors[1]

# rLASSO #
#Note the final input are the coefficients from results of rLASSO:
rlasso1c=rlasso.rolling.window.all(Y,nprev,1,1,FALSE)
rlasso3c=rlasso.rolling.window.all(Y,nprev,1,3,FALSE)
rlasso6c=rlasso.rolling.window.all(Y,nprev,1,6,FALSE)
rlasso12c=rlasso.rolling.window.all(Y,nprev,1,12,FALSE)

# rLASSO RMSE's
rlasso.rmse1=rlasso1c$errors[1]
rlasso.rmse3=rlasso3c$errors[1]
rlasso.rmse6=rlasso6c$errors[1]
rlasso.rmse12=rlasso12c$errors[1]



# post-LASSO #
#Note the final input are the coefficients from results of rLASSO:
prlasso1c=rlasso.rolling.window.all(Y,nprev,1,1,TRUE)
prlasso3c=rlasso.rolling.window.all(Y,nprev,1,3,TRUE)
prlasso6c=rlasso.rolling.window.all(Y,nprev,1,6,TRUE)
prlasso12c=rlasso.rolling.window.all(Y,nprev,1,12,TRUE)

# post-LASSO RMSE's
prlasso.rmse1=rlasso1c$errors[1]
prlasso.rmse3=rlasso3c$errors[1]
prlasso.rmse6=rlasso6c$errors[1]
prlasso.rmse12=rlasso12c$errors[1]



########################################################################
#Penalized regression: Elastic Net (BIC, alpha=0.5)
########################################################################

alpha=0.5 #set alpha to 0.5 for elastic net

elasticnet1c=lasso.rolling.window.all(Y,nprev,1,1,alpha,IC="bic")
elasticnet3c=lasso.rolling.window.all(Y,nprev,1,3,alpha,IC="bic")
elasticnet6c=lasso.rolling.window.all(Y,nprev,1,6,alpha,IC="bic")
elasticnet12c=lasso.rolling.window.all(Y,nprev,1,12,alpha,IC="bic")

#See the RMSE
elnet.rmse1=elasticnet1c$errors[1]
elnet.rmse3=elasticnet3c$errors[1]
elnet.rmse6=elasticnet6c$errors[1]
elnet.rmse12=elasticnet12c$errors[1]



#######################################################
##Sparsity analysis over time plots for LASSO and ElNet
#######################################################

#Get nonzero coefficient numbers for different horizons (LASSO(BIC))
c1c=rowSums(prlasso1c$coef != 0)
c3c=rowSums(prlasso3c$coef != 0)
c6c=rowSums(prlasso6c$coef != 0)
c12c=rowSums(prlasso12c$coef != 0)

#Create a ts object for the plot
lcoef.ts=ts(cbind(c1c,c3c,c6c,c12c), start=c(1960,3), end=c(2024,07), freq=12)
colnames(lcoef.ts)=c("1-step","3-step","6-step","12-step")
#Plot numbers of nonzero coefficients across the test window
windows()
plot.ts(lcoef.ts, main="Sparsity Analysis for Post-LASSO",cex.axis=1.5)

#Get nonzero coefficient numbers for different horizons (Elastic Net)
ce1c=rowSums(elasticnet1c$coef != 0)
ce3c=rowSums(elasticnet3c$coef != 0)
ce6c=rowSums(elasticnet6c$coef != 0)
ce12c=rowSums(elasticnet12c$coef != 0)

#Create a respective ts object for the plot
elcoef.ts=ts(cbind(ce1c,ce3c,ce6c,ce12c), start=c(1960,3), end=c(2024,07), freq=12)
colnames(elcoef.ts)=c("1-step","3-step","6-step","12-step")
#Plot numbers of nonzero coefficients across the test window:
windows()
plot.ts(elcoef.ts, main="Sparsity Analysis for ElNet",cex.axis=1.5)

