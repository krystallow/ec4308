library(leaps)
library(tictoc)
library(glmnet)

set.seed(2457829)
ntrain <- floor(0.7 * nrow(df))  # take 70% train
tr <- sample(seq_len(nrow(df)), size = ntrain)
train <- df[tr, ] 
test <- df[-tr, ]  

formula <- target_column ~ . - others_columns  # change 


###########################################################
######## Best Subset Selection >> method = exhaustive 
###########################################################
tic("Best Subset Selection")
best_subset_model <- regsubsets(formula, data = train, nvmax = ncol(train) - 1, method = "exhaustive")
toc()

summary_bss <- summary(best_subset_model)


best_k_bic <- which.min(summary_bss$bic)
best_k_aic <- which.min(summary_bss$cp) 
best_k_adjr2 <- which.max(summary_bss$adjr2) # see model with highest adjusted R-squared >> check if overfit

plot((summary_bss$bic - mean(summary_bss$bic)) / sd(summary_bss$bic), type = "l", col = "blue",
     main = "AIC and BIC", xlab = "Number of predictors (k)", ylab = "Standardized values")
points((summary_bss$cp - mean(summary_bss$cp)) / sd(summary_bss$cp), type = "l", col = "red")
legend("topright", c("BIC", "AIC"), lty = 1, col = c("blue", "red"))


#AIC and BIC using error variance estimate from the largest model
varest=summary_bss$rss[num_params]/(ntrain-num_params+1) #estimate error variance of the model with k 

#construct the IC with this estimate >> gen sequence of 1 to k to plug into formula
BICL = summary_bss$rss/ntrain + log(ntrain)*varest*((seq(1,num_params,1))/ntrain)
AICL = summary_bss$rss/ntrain + 2*varest*((seq(1,num_params,1))/ntrain)

kbicl=which.min(BICL) #BIC choice

kaicl=which.min(AICL)  #AIC choice, same

#AIC and BIC using iterative procedure
## Calculate OOS MSE using AIC and BIC Minimizing models:
#Get the X-matrix for the test set:
test.mat = model.matrix( formula , data = test)

#extract coefficients from the best model on BIC
temp.coef = coef(summary_bss, id = kbic)
MSEBIC = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselbic = names(temp.coef)  # Selected variables on BIC


#Repeat this for AIC >> check if all agree/same
temp.coef = coef(summary_bss, id = kaic)
MSEAIC = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaic = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(summary_bss, id = kbicl)
MSEBICL = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = kaicl)
MSEAICL = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = k1bic)
MSEBIC1 = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = k1aic)
MSEAIC1 = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)


############################
###10-fold Cross Validation for Best Subset
############################

Kf=10 

split = runif(ntrain) 
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for K-fold cv

train.mat = model.matrix(target_column~.-other_columns, data=train) 

cv_bss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(target_column~.-other_columns, data = train[nii,], nvmax = num_params, method = "exhaustive") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_bss[j,i] = sum((train$target_column[ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestK = which.min(colSums(cv_bss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvbss_min = min(colSums(cv_bss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_bss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bssel, id = bestK) #extract the coefficient vector of the best model
MSECV = mean((test$Balance-test.mat[,names(temp.coef)]%*%temp.coef)^2)



###########################################################
######## Forward Selection >> method = forward 
###########################################################
tic("Forward Selection")
bfssel <- regsubsets(formula, data = train, nvmax = ncol(train) - 1, method = "forward")
toc()

sumbfs <- summary(bfssel)


best_k_bic <- which.min(sumbfs$bic)
best_k_aic <- which.min(sumbfs$cp) 
best_k_adjr2 <- which.max(sumbfs$adjr2) # see model with highest adjusted R-squared >> check if overfit

plot((sumbfs$bic - mean(sumbfs$bic)) / sd(sumbfs$bic), type = "l", col = "blue",
     main = "AIC and BIC", xlab = "Number of predictors (k)", ylab = "Standardized values")
points((sumbfs$cp - mean(sumbfs$cp)) / sd(sumbfs$cp), type = "l", col = "red")
legend("topright", c("BIC", "AIC"), lty = 1, col = c("blue", "red"))


#AIC and BIC using error variance estimate from the largest model
varestf=sumbfs$rss[num_params]/(ntrain-num_params+1) #estimate error variance of the model with k 

#construct the IC with this estimate >> gen sequence of 1 to k to plug into formula
BICLf = sumbfs$rss/ntrain + log(ntrain)*varestf*((seq(1,num_params,1))/ntrain)
AICLf = sumbfs$rss/ntrain + 2*varestf*((seq(1,num_params,1))/ntrain)

kbiclf=which.min(BICLf) #BIC choice

kaiclf=which.min(AICLf)  #AIC choice, same

#AIC and BIC using iterative procedure
## Calculate OOS MSE using AIC and BIC Minimizing models:
#Get the X-matrix for the test set:
test.mat = model.matrix( formula , data = test)

#extract coefficients from the best model on BIC
temp.coef = coef(summary_bss, id = kbicf)
MSEBICf = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselbic = names(temp.coef)  # Selected variables on BIC


#Repeat this for AIC >> check if all agree/same
temp.coef = coef(summary_bss, id = kaicf)
MSEAICf = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaic = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(summary_bss, id = kbiclf)
MSEBICLf = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = kaiclf)
MSEAICfL = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = k1bicf)
MSEBIC1f = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(summary_bss, id = k1aicf)
MSEAIC1f = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)


############################
###10-fold Cross Validation for Forward
############################

Kf=10 

split = runif(ntrain) 
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for K-fold cv

train.mat = model.matrix(target_column~.-other_columns, data=train) 

cv_fss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(target_column~.-other_columns, data = train[nii,], nvmax = num_params, method = "forward") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_fss[j,i] = sum((train$target_column[ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestK = which.min(colSums(cv_fss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvfss_min = min(colSums(cv_fss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_fss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bssel, id = bestK) #extract the coefficient vector of the best model
MSECVf = mean((test$Balance-test.mat[,names(temp.coef)]%*%temp.coef)^2)



###########################################################
######## Backward Selection >> method = backward 
###########################################################
tic("Backward Selection")
bbssel <- regsubsets(formula, data = train, nvmax = ncol(train) - 1, method = "backward")
toc()

sumbbs <- summary(bbssel)


best_k_bic <- which.min(sumbbs$bic)
best_k_aic <- which.min(sumbbs$cp) 
best_k_adjr2 <- which.max(sumbbs$adjr2) # see model with highest adjusted R-squared >> check if overfit

plot((sumbbs$bic - mean(sumbbs$bic)) / sd(sumbbs$bic), type = "l", col = "blue",
     main = "AIC and BIC", xlab = "Number of predictors (k)", ylab = "Standardized values")
points((sumbbs$cp - mean(sumbbs$cp)) / sd(sumbbs$cp), type = "l", col = "red")
legend("topright", c("BIC", "AIC"), lty = 1, col = c("blue", "red"))


#AIC and BIC using error variance estimate from the largest model
varestb=sumbbs$rss[num_params]/(ntrain-num_params+1) #estimate error variance of the model with k 

#construct the IC with this estimate >> gen sequence of 1 to k to plug into formula
BICLb = sumbbs$rss/ntrain + log(ntrain)*varestb*((seq(1,num_params,1))/ntrain)
AICLb = sumbbs$rss/ntrain + 2*varestb*((seq(1,num_params,1))/ntrain)

kbiclb=which.min(BICLb) #BIC choice

kaiclb=which.min(AICLb)  #AIC choice, same

#AIC and BIC using iterative procedure
## Calculate OOS MSE using AIC and BIC Minimizing models:
#Get the X-matrix for the test set:
test.mat = model.matrix( formula , data = test)

#extract coefficients from the best model on BIC
temp.coef = coef(sumbbs, id = kbicb)
MSEBICb = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselbicb = names(temp.coef)  # Selected variables on BIC


#Repeat this for AIC >> check if all agree/same
temp.coef = coef(sumbbs, id = kaicb)
MSEAICb = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaicb = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(sumbbs, id = kbiclb)
MSEBICLb = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(sumbbs, id = kaiclb)
MSEAICLb = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(sumbbs, id = k1bicb)
MSEBIC1b = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(sumbbs, id = k1aicb)
MSEAIC1b = mean((test$target_column-test.mat[,names(temp.coef)]%*%temp.coef)^2)


############################
###10-fold Cross Validation for Backward
############################

Kf=10 

split = runif(ntrain) 
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for K-fold cv

train.mat = model.matrix(target_column~.-other_columns, data=train) 

cv_bss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(target_column~.-other_columns, data = train[nii,], nvmax = num_params, method = "backward") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_bss[j,i] = sum((train$target_column[ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestKb = which.min(colSums(cv_bss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvbss_min = min(colSums(cv_bss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_bss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bssel, id = bestK) #extract the coefficient vector of the best model
MSECVb = mean((test$Balance-test.mat[,names(temp.coef)]%*%temp.coef)^2)



###########################################################
######## Ridge
###########################################################

x = model.matrix(target_column ~ .-other_columns, data = df)
x = x[, -1]  #omit the intercept 
y = df$target_column

MSE <- function(pred, truth){ 
  return(mean((truth - pred)^2))
}
grid = 10^seq(10, -2, length = 100)

ridge.mod <- glmnet(x[train,], y[train], alpha = 0, lambda = grid, thresh = 1e-12) # check test


## Choose Lambda by 10-fold CV with auto defined grid:
cv.out10d = cv.glmnet(x[train,], y[train], alpha = 0) 
plot(cv.out10d) 
plot(cv.out10d$lambda,cv.out10d$cvm, main="10-fold CV default settings", xlab="Lambda", ylab="CV MSE")
bestlam10d = cv.out10d$lambda.min
bestlam10d

# RMSE associated with the value of lambda chosen by cross-validation >> AUTO GRID
ridge.pred <- predict(ridge.mod, s = bestlam10d, newx = x[test,])
MSE(ridge.pred, y[test])

## Choose Lambda with 10-fold CV with user defined grid:
cv.out10u = cv.glmnet(x[train,], y[train], alpha = 0, lambda=grid)
plot(cv.out10u)
plot(cv.out10d$lambda,cv.out10d$cvm, main="10-fold CV user-defined grid", xlab="Lambda", ylab="CV MSE") 
bestlam10u = cv.out10u$lambda.min
bestlam10u

# RMSE associated with the value of lambda chosen by cross-validation >> USER GRID
ridge.pred <- predict(ridge.mod, s = bestlam10u, newx = x[test,])
MSE(ridge.pred, y[test])


# Check LOOCV
cv.outl = cv.glmnet(x[train,], y[train], alpha = 0, lambda=grid, nfolds=300)
plot(cv.outl$lambda,cv.outl$cvm, main="LOOCV user-defined grid", xlab="Lambda", ylab="CV MSE")
bestlaml = cv.outl$lambda.min
bestlaml

# RMSE associated with the value of lambda chosen by cross-validation >> LOOCV USER GRID
ridge.pred <- predict(ridge.mod, s = bestlaml, newx = x[test,])
MSE(ridge.pred, y[test])


# define a finer grid just to check values closer to 0 (100 equally spaced point from 0.001 to 5)
grid2 = seq(0.001, 5, length=100) ## change 

# Re-Check LOOCV with the new grid
cv.outl = cv.glmnet(x[train,], y[train], alpha = 0, lambda=grid2, nfolds=300)
plot(cv.outl$lambda,cv.outl$cvm, main="LOOCV user-defined grid", xlab="Lambda", ylab="CV MSE")
bestlaml = cv.outl$lambda.min
bestlaml

# RMSE associated with the value of lambda chosen by cross-validation >> LOOCV with NEW grid
ridge.pred <- predict(ridge.mod, s = bestlaml, newx = x[test,])
MSE(ridge.pred, y[test])




###########################################################
######## Ridge-PCA Hybrid Approach
###########################################################

x = model.matrix(target_column ~ .-other_columns, data = df)  # Create matrix of predictors
x = x[, -1] 
y = df$target_column

train_index <- sample(1:nrow(x), nrow(x) * 0.7) 
train <- train_index
test <- -train_index

# standardise before PCA
pca_model <- preProcess(x[train,], method = c("center", "scale", "pca"))

x_train_pca <- predict(pca_model, x[train,])
x_test_pca <- predict(pca_model, x[test,])

grid = 10^seq(10, -2, length = 100) ## define again further on
ridge_mod <- glmnet(as.matrix(x_train_pca), y[train], alpha = 0, lambda = grid, thresh = 1e-12)

# choose lambda based on 10 fold CV
cv_out = cv.glmnet(as.matrix(x_train_pca), y[train], alpha = 0, lambda = grid, nfolds = 10)
best_lambda = cv_out$lambda.min

plot(cv_out)
plot(cv_out$lambda, cv_out$cvm, main = "10-fold CV Ridge-PCA", xlab = "Lambda", ylab = "CV MSE")

ridge_pred <- predict(ridge_mod, s = best_lambda, newx = as.matrix(x_test_pca))

MSE_ridge_pca <- MSE(ridge_pred, y[test])

rmse_ridge_pca <- sqrt(test_mse)


