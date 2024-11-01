library(glmnet)
library(hdm)
library(pls)
library(readxl)

set.seed(2457829)
library(readxl)
data <- read_excel("Final_Transformed_Dataset.xlsx", range = "B1:DK775")
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
                                                      "numeric"), range = "A1:DK775" )
data = data[-1,]
df = df[-1,]
data_clean <- na.omit(data)
df = na.omit(df)
# Now create model matrix
x = model.matrix(INDPRO~., data = data_clean)
y = data_clean$INDPRO
ntrain= floor(0.7* 773)  #number of training data -> 70% of total obs
train = sample(1:nrow(x),ntrain)
test = (-train)
training_set = data_clean[train,]
testing_set = data_clean[-train,]
MSE <- function(pred, truth){
  return(mean((truth - pred)^2)) 
}

  
################################
# LASSO
################################
lasso.mod <- glmnet(x[train,], y[train], alpha = 1)
plot(lasso.mod) #this allows you to plot coefficients as a function of your L1-budget (and hence lambda)

#This is may be useful to see how and where shrinkage kicks in.
cv.out10u = cv.glmnet(x[train,], y[train], alpha = 1)
plot(cv.out10u)
plot(cv.out10u$lambda,cv.out10u$cvm, main="10-fold CV", xlab="Lambda", ylab="CV MSE")
bestlam10u = cv.out10u$lambda.min
lasso.pred <- predict(lasso.mod, s = bestlam10u, newx = x[test,])
mse_lasso_cv = MSE(lasso.pred, y[test])

lasso.selected <- glmnet(x[train, ], y[train], alpha = 1, lambda = bestlam10u)
lasso_coef <- coef(lasso.selected)
# Convert to a data frame for easier interpretation
selected_vars <- as.data.frame(as.matrix(lasso_coef))
selected_vars <- selected_vars[selected_vars[, 1] != 0, , drop = FALSE]

# Rename the columns for clarity
colnames(selected_vars) <- c("Coefficient")

# View the selected variables
print(selected_vars)

# Plug in lambda choice
rlasso.fit = rlasso(y[train]~x[train,],  post=FALSE)
yhat.rlasso<- predict(rlasso.fit, newdata=x[test,])
mse_plugin_LASSO = MSE(yhat.rlasso, y[test])

variables_selected <- colnames(data[, colnames(data) %in% rownames(selected_vars)[-1]])
variables_selected_df <- data[, variables_selected]
write.csv(variables_selected_df, file = "lasso_variables_selected.csv", row.names = FALSE)

################################
# POST LASSO
################################

# Set-up for cross-validation
Kf = 10  # 10 folds
split = runif(ntrain) #uniform random number for assigning folds
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for 10-fold CV

cv.lambda = cv.out10u$lambda #extract lambdas used for our 10-fold CV previously
nlam = length(cv.lambda) #count the number of lambdas

#Set up training and test sets in matrices for convenience (as we will loop over the elements)
trainfs.mat=x[train,]
testfs.mat=x[test,]
ytr=y[train]
yte=y[test]
cv_pl = matrix(0,Kf,nlam)
for(j in 1:Kf) {
  ii = cvgroup == j
  nii = cvgroup != j
  temp.lasso = glmnet(trainfs.mat[nii,],ytr[nii],lambda = cv.lambda) #NB: run the usual LASSO for the grid of lambdas
  for(i in 1:nlam) {
    temp.coef = abs(predict(temp.lasso, s = cv.lambda[i], type = "coef")) #extract coefficients
    temp.use = which(temp.coef > 0)-1
    temp.df = data.frame(out = ytr[nii], W = trainfs.mat[nii,temp.use]) #extract data for included coefficients
    temp.ols = lm(out~., data = temp.df) #run post-LASSO OLS with the retained predictors 
    temp.testdf = data.frame(out = ytr[ii], W = trainfs.mat[ii,temp.use]) #temporary test data for retained predictors
    cv_pl[j,i] = sum((ytr[ii] - predict(temp.ols, newdata = temp.testdf))^2) #compute CV, but note it is compute for post-LASSO prediction each time - this is the key point!
  }
}
lambdapl=cv.lambda[which.min(colSums(cv_pl))] #extract the lambda minimizing CV criterion - note much smaller than the lambda in the usual LASSO CV!
plasso1 = glmnet(trainfs.mat,ytr,lambda = cv.lambda)
temp.use = which(abs(predict(plasso1, s = lambdapl, type = "coef")) > 0)-1 #find nonzero coefficients and create variable indices
# (note -1 here reduces indices to adjust for the absence of constant in trainfs.mat)
temp.df = data.frame(out = ytr, W = trainfs.mat[,temp.use]) #create training dataframe with retained predictors
temp.ols = lm(out~., data = temp.df) #run post-LASSO OLS
temp.testdf = data.frame(out = yte, W = testfs.mat[,temp.use]) #create the test dataframe with retained predictors
MSEpl = mean((yte - predict(temp.ols, newdata = temp.testdf))^2) 
MSEpl

################################
# Elastic Net
################################

# Function to run Elastic Net on multiple alphas
run_elastic_net <- function(alpha_val) {
  elnet.mod <- glmnet(x[train,], y[train], alpha = alpha_val)
  cv.out10el = cv.glmnet(x[train,], y[train], alpha = alpha_val)
  
  plot(cv.out10el$lambda, cv.out10el$cvm, main = paste("CV for alpha =", alpha_val), xlab = "Lambda", ylab = "CV MSE")
  
  best_lambda = cv.out10el$lambda.min
  elnet.pred <- predict(elnet.mod, s = best_lambda, newx = x[test,])
  
  mse_alpha = MSE(elnet.pred, y[test])
  return(list(best_lambda = best_lambda, mse = mse_alpha))
}

# Run Elastic Net for various alpha values
alpha_values = seq(0.1, 0.9, by = 0.1)
results = lapply(alpha_values, run_elastic_net)

for (i in seq_along(results)) {
  cat(paste("Elastic Net with alpha =", alpha_values[i], "-> Best Lambda:", results[[i]]$best_lambda, ", MSE:", results[[i]]$mse, "\n"))
}

# Elastic Net with alpha = 0.9 -> Best Lambda: 2.40465794564183e-05 , MSE: 1.7853529881859e-07 
mse_elasticnet = 1.7853529881859e-07 
#################################
# Principal components regression
##################################

#Examine the fit syntax:
pcr.fit=pcr(y~.,data=data_clean, subset=train, scale=TRUE, validation="CV")

#Summary of CV results and variance shares explained by components
summary(pcr.fit)

#Plot loadings for the first 2 components to attempt interpretation
plot(pcr.fit, "loadings", comps = 1:2, legendpos = "topleft")
abline(h = 0) #add the zero line for reference

#Built-in plot command to display CV results
validationplot(pcr.fit, val.type="MSEP", main="CV",legendpos = "topright") 

#Use built-in selection plot for M based on min CV MSE as well as the 1SE rule of thumb
ncomp.onesigma = selectNcomp(pcr.fit, method = "onesigma", plot = TRUE, main="1 SE rule selection")

pcr.pred=predict(pcr.fit, newdata=x[test,])
mse_pcr = MSE(pcr.pred,y[test])

##########################################
# Principal components regression + LASSO
##########################################
comps.mat=pcr.fit$scores
cv.out = cv.glmnet(comps.mat, y[train], alpha = 1) #determine lambda.min and lambda.1SE
lasso.mod <- glmnet(comps.mat, y[train], alpha = 1, lambda = cv.out$lambda.1se)
lasso.mod$beta #check the coefficients

# predict with the LASSO model using PCR components as input
# Attempt to extract components, ensuring it maintains a matrix structure
comps.test <- pcr.pred[, , 1]  # If this gives NULL, try the following:
comps.test <- matrix(pcr.pred[, , 1], nrow = dim(pcr.pred)[1], ncol = dim(pcr.pred)[3])  # Convert it to a 2D matrix
yhat.lasso_pcr <- predict(lasso.mod, newx = comps.test)
mse_lasso_pcr_cv = MSE(yhat.lasso_pcr, y[test])

# Append results to the data frame
mse_results <- data.frame(
  Model = character(),
  MSE = numeric(),
  stringsAsFactors = FALSE
)

mse_results <- rbind(mse_results, 
                     data.frame(Model = "LASSO CV", MSE = mse_lasso_cv),
                     data.frame(Model = "LASSO plugin", MSE = mse_plugin_LASSO),
                     data.frame(Model = "LASSO + PCR", MSE = mse_lasso_pcr_cv),
                     data.frame(Model = "POST LASSO", MSE = MSEpl),
                     data.frame(Model = "PCR", MSE = mse_pcr))

# View the results
print(mse_results)

###############
# FORECASTING
###############
library(githubinstall)
githubinstall("HDeconometrics")
library(HDeconometrics)
library(sandwich)
df$sasdate <- as.Date(df$sasdate) 
dum <- ifelse(format(df$sasdate, "%Y") == "2020", 1, 0)
df <- cbind(df, dum = dum)

nprev=230 #number of out-of-sample observations (test window )

oosy=tail(y,nprev) #auxiliary:get the out-of-sample true values (last 180 obs. using tail())

RMSE <- function(pred, truth){
  return(sqrt(mean((truth - pred)^2)))
}

#create lags
rwtemp=embed(y,13)
#Simple RW forecast:
rw1c=tail(rwtemp[,2],nprev)
rw3c=tail(rwtemp[,4],nprev)
rw6c=tail(rwtemp[,7],nprev)
rw12c=tail(rwtemp[,13],nprev)

#Collect RMSE's for randomw walk:
rw.rmse1=RMSE(oosy,rw1c)
rw.rmse3=RMSE(oosy,rw3c)
rw.rmse6=RMSE(oosy,rw6c)
rw.rmse12=RMSE(oosy,rw12c)

################
# Random Forest
################

library(randomForest)
runrf=function(Y,indice,lag){
  
  dum=Y[,ncol(Y)] # extract dummy from data
  Y=Y[,-ncol(Y)]  #data without the dummy
  comp=princomp(scale(Y,scale=FALSE)) # compute principal components to add as predictors
  Y2=cbind(Y,comp$scores[,1:35]) #augment predictors by the first 35 principal components
  aux=embed(Y2,6+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))]  # lags of Y (predictors) corresponding to forecast horizon 
  
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)]   #retrieve the last  observations if one-step forecast
  }else{
    X.out=aux[,-c(1:(ncol(Y2)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)]  #last observations: y_T,y_t-1...y_t-h
  }
  
  dum=tail(dum,length(y)) #cut the dummy to size to account for lost observations due to lags
  
  model=randomForest(cbind(X,dum),y,importance = TRUE) #fit the random forest on default settings
  pred=predict(model,c(X.out,0)) #generate forecast
  
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}


#This function will repeatedly call the previous function in the rolling window h-step forecasting

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) nprev - number of out-of-sample observations (at the end of the sample)

#3) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#4) lag - the forecast horizon

rf.rolling.window=function(Y,nprev,indice=1,lag=1){
  
  save.importance=list() #blank for saving variable importance
  save.pred=matrix(NA,nprev,1) ##blank for forecasts
  for(i in nprev:1){#NB: backwards FOR loop: going from 180 down to 1
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 491, then 2 to 492 etc.)
    rf=runrf(Y.window,indice,lag)#call the function to fit the Random Forest and generate h-step forecast
    save.pred[(1+nprev-i),]=rf$pred #save the forecast
    save.importance[[i]]=importance(rf$model) #save variable importance
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  #Some helpful stuff:
  real=Y[,indice]#get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  
  return(list("pred"=save.pred,"errors"=errors,"save.importance"=save.importance)) #return forecasts, history of variable importance, and RMSE and MAE for the period.
}


df = df[,-1]
df = as.matrix(df)
rf1c=rf.rolling.window(df,nprev,1,1)
rf3c=rf.rolling.window(df,nprev,1,3)
rf6c=rf.rolling.window(df,nprev,1,6)
rf12c=rf.rolling.window(df,nprev,1,12)

par(mfrow = c(2, 2)) # Arrange plots in a 2x2 layout

# Define actual values
real_values <- df[,1] # or whichever variable corresponds to the actual values in df

# Plot for 1-step forecast
plot(real_values, type = "l", main = "1-Step Forecast", ylab = "Real", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf1c$pred)), rf1c$pred), col = "red")

# Plot for 3-step forecast
plot(real_values, type = "l", main = "3-Step Forecast", ylab = "Real", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf3c$pred)), rf3c$pred), col = "red")

# Plot for 6-step forecast
plot(real_values, type = "l", main = "6-Step Forecast", ylab = "Real", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf6c$pred)), rf6c$pred), col = "red")

# Plot for 12-step forecast
plot(real_values, type = "l", main = "12-Step Forecast", ylab = "Real", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf12c$pred)), rf12c$pred), col = "red")

# Restore plotting layout
par(mfrow = c(1, 1))


#See the RMSE:
rf.rmse1=rf1c$errors[1]
rf.rmse3=rf3c$errors[1]
rf.rmse6=rf6c$errors[1]
rf.rmse12=rf12c$errors[1]

###########
# Bagging
###########

library(readr)
variables_selected <- read_csv("variables_selected.csv")
# Simulate future feature data for the next 12 months (e.g., based on historical data or trends)
future_IPFPNSS <- rnorm(12, mean = mean(variables_selected$IPFPNSS), sd = sd(variables_selected$IPFPNSS))
future_IPMAT <- rnorm(12, mean = mean(variables_selected$IPMAT), sd = sd(variables_selected$IPMAT))
future_IPDMAT <- rnorm(12, mean = mean(variables_selected$IPDMAT), sd = sd(variables_selected$IPDMAT))
future_IPMANSICS <- rnorm(12, mean = mean(variables_selected$IPMANSICS), sd = sd(variables_selected$IPMANSICS))
future_UNRATE <- rnorm(12, mean = mean(variables_selected$UNRATE), sd = sd(variables_selected$UNRATE))
future_CES1021000001 <- rnorm(12, mean = mean(variables_selected$CES1021000001), sd = sd(variables_selected$CES1021000001))
future_AWOTMAN <- rnorm(12, mean = mean(variables_selected$AWOTMAN), sd = sd(variables_selected$AWOTMAN))
future_CES0600000008 <- rnorm(12, mean = mean(variables_selected$CES0600000008), sd = sd(variables_selected$CES0600000008))
future_CES3000000008 <- rnorm(12, mean = mean(variables_selected$CES3000000008), sd = sd(variables_selected$CES3000000008))
future_CIVPART <- rnorm(12, mean = mean(variables_selected$CIVPART), sd = sd(variables_selected$CIVPART))
future_LNS11300036 <- rnorm(12, mean = mean(variables_selected$LNS11300036), sd = sd(variables_selected$LNS11300036))
future_FEDFUNDS <- rnorm(12, mean = mean(variables_selected$FEDFUNDS), sd = sd(variables_selected$FEDFUNDS))

# Combine into a future data frame
future_variables_selected <- data.frame(
  IPFPNSS = future_IPFPNSS,
  IPMAT = future_IPMAT,
  IPDMAT = future_IPDMAT,
  IPMANSICS = future_IPMANSICS,
  UNRATE = future_UNRATE,
  CES1021000001 = future_CES1021000001,
  AWOTMAN = future_AWOTMAN,
  CES0600000008 = future_CES0600000008,
  CES3000000008 = future_CES3000000008,
  CIVPART = future_CIVPART,
  LNS11300036 = future_LNS11300036,
  FEDFUNDS = future_FEDFUNDS
)

# Create a data frame for 12 months
months <- paste(1:230)
future_values <- rep(NA, 230)  # Placeholder for future values
baggingfit = randomForest(y~.,data=variables_selected,ntree=5000, mtry=230)
plot(baggingfit) #plot the last fitted (largest) OOB error
bagging_prediction = predict(baggingfit, newdata = future_variables_selected)

prediction_df_bagging <- data.frame(Month = months, Predicted_Values = bagging_prediction)


################
# Random Forest
################
library(fanplot)
rffit = randomForest(y~.,data=variables_selected,ntree=5000)
plot(rffit) #plot the last fitted (largest) OOB error
rf_prediction = predict(rffit, newdata = future_variables_selected)

# Random Forest Fan Chart
# Simulate multiple scenarios for each month by adding random noise
set.seed(123)
num_scenarios <- 100  # Number of scenarios to simulate
INDPRO = rf_prediction
predicted_values_matrix <- replicate(num_scenarios, rf_prediction + rnorm(12, 0, 0.002))
# Convert matrix into time series format
predicted_values_ts <- ts(predicted_values_matrix, start = c(1, 1), frequency = 12)
plot(INDPRO, main="Fan Chart for Random Forest", xlim = c(1, 12), ylim = c(-0.02, 0.02))
fan(predicted_values_ts, data.type = "simulations", 
    probs = seq(0.1, 0.9, by = 0.1), 
    fan.col = colorRampPalette(c("blue", "white")))

# Bagging Fan Chart
# Simulate multiple scenarios for each month by adding random noise
set.seed(123)
num_scenarios <- 100  # Number of scenarios to simulate
INDPRO = bagging_prediction
predicted_values_matrix <- replicate(num_scenarios, bagging_prediction + rnorm(12, 0, 0.002))
# Convert matrix into time series format
predicted_values_ts <- ts(predicted_values_matrix, start = c(1, 1), frequency = 12)
plot(INDPRO, main="Fan Chart for Random Forest", xlim = c(1, 12), ylim = c(-0.02, 0.02))
fan(predicted_values_ts, data.type = "simulations", 
    probs = seq(0.1, 0.9, by = 0.1), 
    fan.col = colorRampPalette(c("blue", "white")))



