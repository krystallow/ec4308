library(glmnet)
library(hdm)
library(pls)
library(readxl)

set.seed(2457829)
data <- read_excel("Final_Transformed_Dataset.xlsx", range = "B1:DK775")
data = data[-1,]
data_clean <- na.omit(data)
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

###########
# Bagging
###########
library(randomForest)
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
months <- paste(1:12)
future_values <- rep(NA, 12)  # Placeholder for future values
baggingfit = randomForest(y~.,data=variables_selected,ntree=5000, mtry=12)
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



