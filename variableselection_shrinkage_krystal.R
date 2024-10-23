library(leaps)
library(tictoc)
library(glmnet)
library(tree)
library(rpart)
library(caret)
library(readxl)
library(gbm)

Transformed_Dataset <- read_excel("Final_Transformed_Dataset.xlsx")
Transformed_Dataset <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))


df <- Transformed_Dataset
df <- df[-1, ]

set.seed(2457829)
ntrain <- floor(0.7 * nrow(df))  # take 70% train
tr = sample(1:nrow(df),ntrain)  
train = df[tr,]   
test = df[-tr,]  

## CHANGE
target_column <- "INDPRO"
others_columns <- c("sasdate")
formula <- as.formula(paste(target_column, "~.-", paste(others_columns, collapse = "-")))
##

num_params = ncol(train)-1

##########################################################
######## Best Subset Selection >> method = exhaustive ### >> Not Ran, too long run time
##########################################################
tic("Best Subset Selection")
best_subset_model <- regsubsets(formula, data = train, nvmax = ncol(train)-1, method = "exhaustive", really.big = TRUE)
toc()

summary_bss <- summary(best_subset_model)

kbic <- which.min(summary_bss$bic)
kaic <- which.min(summary_bss$cp) 
which.max(summary_bss$adjr2) 

plot((summary_bss$bic - mean(summary_bss$bic)) / sd(summary_bss$bic), type = "l", col = "blue",
     main = "AIC and BIC", xlab = "Number of predictors (k)", ylab = "Standardized values")
points((summary_bss$cp - mean(summary_bss$cp)) / sd(summary_bss$cp), type = "l", col = "red")
legend("topright", c("BIC", "AIC"), lty = 1, col = c("blue", "red"))


### AIC and BIC using error variance estimate from the largest model ### 
varest=summary_bss$rss[num_params]/(ntrain-num_params+1) #estimate error variance of the model with k 

#construct the IC with this estimate >> gen sequence of 1 to k to plug into formula
BICL = summary_bss$rss/ntrain + log(ntrain)*varest*((seq(1,num_params,1))/ntrain)
AICL = summary_bss$rss/ntrain + 2*varest*((seq(1,num_params,1))/ntrain)

kbicl=which.min(BICL) #BIC choice
kaicl=which.min(AICL)  #AIC choice, same


### AIC and BIC using iterative procedure ###
## Calculate OOS MSE using AIC and BIC Minimizing models:
#Get the X-matrix for the test set:
test.mat = model.matrix(formula, data = test)

#extract coefficients from the best model on BIC
temp.coef = coef(best_subset_model, id = kbicl)
MSEBIC = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)
varselbic = names(temp.coef)  # Selected variables on BIC

#Repeat this for AIC >> check if all agree/same
temp.coef = coef(best_subset_model, id = kaicl)
MSEAIC = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaic = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(best_subset_model, id = kbicl)
MSEBICL = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(best_subset_model, id = kaicl)
MSEAICL = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(best_subset_model, id = k1bic)
MSEBIC1 = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(best_subset_model, id = k1aic)
MSEAIC1 = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)


###############################################
###10-fold Cross Validation for Best Subset ###
###############################################

Kf=10 

split = runif(ntrain) 
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for K-fold cv

train.mat = model.matrix(formula, data=train) 

cv_bss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(formula, data = train[nii,], nvmax = num_params, method = "exhaustive") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_bss[j,i] = sum((train[[target_column]][ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestK = which.min(colSums(cv_bss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvbss_min = min(colSums(cv_bss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_bss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bssel, id = bestK) #extract the coefficient vector of the best model
MSECV = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)



#################################################
######## Forward Selection >> method = forward### 
#################################################
tic("Forward Selection")
bfssel <- regsubsets(formula, data = train, nvmax = num_params, method = "forward")
toc()

sumbfs <- summary(bfssel)

kbicf <- which.min(sumbfs$bic)
kaicf <- which.min(sumbfs$cp) 
best_k_adjr2 <- which.max(sumbfs$adjr2) # see model with highest adjusted R-squared >> check if overfit

plot((sumbfs$bic-mean(sumbfs$bic))/sd(sumbfs$bic), main = "AIC and BIC (Standardized)", xlab = "k",
     ylab = "IC", type = "l", col = "blue")
points((sumbfs$cp-mean(sumbfs$cp))/sd(sumbfs$cp), type = "l", col = "red")
legend("topright", c("BIC","AIC"),lty=c(1,1) ,col=c("blue","red"))

#AIC and BIC using error variance estimate from the largest model
### Estimate error variance of the largest model selected
varestf <- sumbfs$rss[length(sumbfs$rss)] / (ntrain - length(sumbfs$rss) - 1)

### Construct the BIC using error variance estimate and the correct model sizes
BICLf <- sumbfs$rss / ntrain + log(ntrain) * varestf * ((seq(1, length(sumbfs$rss), 1)) / ntrain)
### Plot the standardized BIC for comparison
plot(BICLf, type = "l", col = "green", main = "BIC with Error Variance Estimate", 
     xlab = "Number of Predictors", ylab = "BIC")
### Construct the AIC using error variance estimate and the correct model sizes
AICLf <- sumbfs$rss / ntrain + 2 * varestf * ((seq(1, length(sumbfs$rss), 1)) / ntrain)
### Plot the standardized AIC for comparison
plot(AICLf, type = "l", col = "orange", main = "AIC with Error Variance Estimate", 
     xlab = "Number of Predictors", ylab = "AIC")

kbiclf=which.min(BICLf) #BIC choice

kaiclf=which.min(AICLf)  #AIC choice, same


#######################################
#AIC and BIC using iterative procedure
#######################################
sig0=var(train$INDPRO)
BIC0f = sumbfs$rss/ntrain + log(ntrain)*sig0*((seq(1,11,1))/ntrain)
AIC0f = sumbfs$rss/ntrain + 2*sig0*((seq(1,11,1))/ntrain)

k0bicf = which.min(BIC0f)
k0aicf = which.min(AIC0f)

BIC1f = sumbfs$rss/ntrain + log(ntrain)*(sumbfs$rss[k0bicf]/(ntrain-k0bicf-1))*((seq(1,11,1))/ntrain)
AIC1f = sumbfs$rss/ntrain + 2*(sumbfs$rss[k0aicf]/(ntrain-k0aicf-1))*((seq(1,11,1))/ntrain)

k1bicf = which.min(BIC1f)
k1aicf = which.min(AIC1f)

BIC2f = sumbfs$rss/ntrain + log(ntrain)*(sumbfs$rss[k1bicf]/(ntrain-k1bicf-1))*((seq(1,11,1))/ntrain)
AIC2f = sumbfs$rss/ntrain + 2*(sumbfs$rss[k1aicf]/(ntrain-k1aicf-1))*((seq(1,11,1))/ntrain)

k2bicf = which.min(BIC2f)
k2aicf = which.min(AIC2f)
#See here we can stop after two steps (i.e., with BIC1/AIC1 - the choice converges to 4 parameters)


## Calculate OOS MSE using AIC and BIC Minimizing models
test.mat = model.matrix(formula, data = test)
temp.coef = coef(bfssel, id = kbicf)
MSEBICf = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselbicf = names(temp.coef)  # Selected variables on BIC

temp.coef = coef(bfssel, id = kaicf)
MSEAICf = mean((test$INDPRO-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaic = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(bfssel, id = kbiclf)
MSEBICLf = mean((test$INDPRO-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(bfssel, id = kaiclf)
MSEAICLf = mean((test$INDPRO-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(bfssel, id = k1bicf)
MSEBIC1f = mean((test$INDPRO-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(bfssel, id = k1aicf)
MSEAIC1f = mean((test$INDPRO-test.mat[,names(temp.coef)]%*%temp.coef)^2)



##########################################
###10-fold Cross Validation for Forward###
##########################################

Kf=10 

split = runif(ntrain) 
cvgroup = as.numeric(cut(split,quantile(split, probs = seq(0,1,.1)),include.lowest = TRUE))  # groups for K-fold cv

train.mat = model.matrix(formula, data=train) 

cv_fss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(formula, data = train[nii,], nvmax = num_params, method = "forward") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_fss[j,i] = sum((train[[target_column]][ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestKf = which.min(colSums(cv_fss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvfss_min = min(colSums(cv_fss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_fss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bfssel, id = bestKf) #extract the coefficient vector of the best model
MSECVf = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

kbicf <- which.min(sumbfs$bic)  # Best model according to BIC
temp.coef_bic <- coef(bfssel, id = kbicf)  # Extract coefficients
varselbic <- names(temp.coef_bic)  # Selected variables for BIC

kaicf <- which.min(sumbfs$cp)  # Best model according to AIC
temp.coef_aic <- coef(bfssel, id = kaicf)  # Extract coefficients
varselaic <- names(temp.coef_aic)  # Selected variables for AIC

bestKf <- which.min(colSums(cv_fss))  # Best model based on cross-validation
temp.coef_cv <- coef(bfssel, id = bestKf)  # Extract coefficients
varselcv <- names(temp.coef_cv)  # Selected variables for CV


# BIC-selected model
cat("Variables selected by BIC:\n", varselbic, "\n")
cat("Coefficients for BIC-selected model:\n", temp.coef_bic, "\n")

# AIC-selected model
cat("Variables selected by AIC:\n", varselaic, "\n")
cat("Coefficients for AIC-selected model:\n", temp.coef_aic, "\n")

# CV-selected model
cat("Variables selected by CV:\n", varselcv, "\n")
cat("Coefficients for CV-selected model:\n", temp.coef_cv, "\n")



##################################################
######## Backward Selection >> method = backward #
##################################################
tic("Backward Selection")
bbssel <- regsubsets(formula, data = train, nvmax = ncol(train) - 1, method = "backward")
toc()

sumbbs <- summary(bbssel)

kbicb <- which.min(sumbbs$bic)
kaicb <- which.min(sumbbs$cp) 
which.max(sumbbs$adjr2) # see model with highest adjusted R-squared >> check if overfit

plot((sumbbs$bic - mean(sumbbs$bic)) / sd(sumbbs$bic), type = "l", col = "blue",
     main = "AIC and BIC", xlab = "Number of predictors (k)", ylab = "Standardized values")
points((sumbbs$cp - mean(sumbbs$cp)) / sd(sumbbs$cp), type = "l", col = "red")
legend("topright", c("BIC", "AIC"), lty = 1, col = c("blue", "red"))


#AIC and BIC using error variance estimate from the largest model
### Estimate error variance of the largest model selected
varestb <- sumbbs$rss[length(sumbbs$rss)] / (ntrain - length(sumbbs$rss) - 1)

### Construct the BIC using error variance estimate and the correct model sizes
BICLb <- sumbbs$rss / ntrain + log(ntrain) * varestb * ((seq(1, length(sumbbs$rss), 1)) / ntrain)
### Plot the standardized BIC for comparison
plot(BICLb, type = "l", col = "green", main = "BIC with Error Variance Estimate", 
     xlab = "Number of Predictors", ylab = "BIC")
### Construct the AIC using error variance estimate and the correct model sizes
AICLb <- sumbbs$rss / ntrain + 2 * varestb * ((seq(1, length(sumbbs$rss), 1)) / ntrain)
### Plot the standardized AIC for comparison
plot(AICLb, type = "l", col = "orange", main = "AIC with Error Variance Estimate", 
     xlab = "Number of Predictors", ylab = "AIC")

kbiclb=which.min(BICLb) #BIC choice

kaiclb=which.min(AICLb)  #AIC choice, same

#AIC and BIC using iterative procedure
## Calculate OOS MSE using AIC and BIC Minimizing models:
test.mat = model.matrix(formula, data = test)
temp.coef = coef(bbssel, id = kbicb)
MSEBICb = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselbicb = names(temp.coef)  # Selected variables on BIC


#Repeat this for AIC >> check if all agree/same
temp.coef = coef(bbssel, id = kaicb)
MSEAICb = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

varselaicb = names(temp.coef)  # Selected variables on AIC

temp.coef = coef(bbssel, id = kbiclb)
MSEBICLb = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)

temp.coef = coef(bbssel, id = kaiclb)
MSEAICLb = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)



############################################
###10-fold Cross Validation for Backward ###
############################################

cv_bss = matrix(0,Kf,num_params) #blank for results with dimensions: number of folds by number of models
for(j in 1:Kf) {
  ii = cvgroup == j #create a logical array with TRUE if cvgroup value belongs to fold number j - a way to identify observations belonging to fold j (test set)
  nii = cvgroup != j #create a logical array with TRUE if cvgroup value DOES NOT belong to (!= denotes "not equal") fold number j - a way to identify observations outside fold j (train set)
  temp.bss = regsubsets(formula, data = train[nii,], nvmax = num_params, method = "backward") #run model search excluding fold j
  for(i in 1:11) {
    temp.coef = coef(temp.bss, id = i)
    cv_bss[j,i] = sum((train[[target_column]][ii]-train.mat[ii,names(temp.coef)]%*%temp.coef)^2) #compute fold MSE for each model (k=1 to k=11) 
  }
}
bestKb = which.min(colSums(cv_bss)) #figure out preferred number of regressors based on some of squared errors (colSums() sums columns of the matrix - the result is a vector of sums)
cvbss_min = min(colSums(cv_bss)/ntrain) #in-sample MSE of CV-selected model (just need to divide by N - squared errors already summed up in the previous line)


plot((colSums(cv_bss)/ntrain), main = "10-fold CV selection", xlab = "k",
     ylab = "10-fold MSE", type = "l", col = "blue")


# get test set MSE of the model chosen by CV:
temp.coef = coef(bbssel, id = bestKb) #extract the coefficient vector of the best model
MSECVb = mean((test[[target_column]]-test.mat[,names(temp.coef)]%*%temp.coef)^2)


# BIC-selected model
kbicb <- which.min(sumbbs$bic)
temp.coef_bicb <- coef(bbssel, id = kbicb)  # Extract coefficients for the BIC-selected model
varsel_bicb <- names(temp.coef_bicb)  # Selected variables for BIC

# Print selected variables and coefficients
cat("Variables selected by BIC:\n", varsel_bicb, "\n")
cat("Coefficients for BIC-selected model:\n", temp.coef_bicb, "\n")

# AIC-selected model
kaicb <- which.min(sumbbs$cp)
temp.coef_aicb <- coef(bbssel, id = kaicb)  # Extract coefficients for the AIC-selected model
varsel_aicb <- names(temp.coef_aicb)  # Selected variables for AIC

# Print selected variables and coefficients
cat("Variables selected by AIC:\n", varsel_aicb, "\n")
cat("Coefficients for AIC-selected model:\n", temp.coef_aicb, "\n")

# CV-selected model
bestKb <- which.min(colSums(cv_bss))  # Best model based on CV
temp.coef_cv <- coef(bbssel, id = bestKb)  # Extract coefficients for the CV-selected model
varsel_cv <- names(temp.coef_cv)  # Selected variables for CV

# Print selected variables and coefficients
cat("Variables selected by CV:\n", varsel_cv, "\n")
cat("Coefficients for CV-selected model:\n", temp.coef_cv, "\n")


#############
### Ridge ###
#############
x <- model.matrix(formula, data = df)
x <- x[, -1]  # Omit the intercept
y <- df[[target_column]]

print(dim(df))
print(dim(x))
print(length(y))

set.seed(2457829)
ntrain <- floor(0.7 * nrow(x))
tr <- sample(seq_len(nrow(x)), ntrain)  
tst <- setdiff(seq_len(nrow(x)), tr)
grid = 10^seq(10, -2, length = 100)

if (all(tr <= nrow(x))) {  
  # Fit the ridge regression model
  ridge.mod <- glmnet(x[tr, ], y[tr], alpha = 0, lambda = grid, thresh = 1e-12)
} else {  
  stop("Some indices in tr are out of bounds for the model matrix x.")
}

MSE <- function(pred, truth){ 
  return(mean((truth - pred)^2))
}

## Choose Lambda by 10-fold CV with auto defined grid:
cv.out10d = cv.glmnet(x[tr,], y[tr], alpha = 0) 
plot(cv.out10d) 
plot(cv.out10d$lambda,cv.out10d$cvm, main="10-fold CV default settings", xlab="Lambda", ylab="CV MSE")
bestlam10d = cv.out10d$lambda.min
bestlam10d

# RMSE associated with the value of lambda chosen by cross-validation >> AUTO GRID
ridge.pred <- predict(ridge.mod, s = bestlam10d, newx = x[tst,])
MSE(ridge.pred, y[tst])

## Choose Lambda with 10-fold CV with user defined grid:
cv.out10u = cv.glmnet(x[tr,], y[tr], alpha = 0, lambda=grid)
plot(cv.out10u)
plot(cv.out10d$lambda,cv.out10d$cvm, main="10-fold CV user-defined grid", xlab="Lambda", ylab="CV MSE") 
bestlam10u = cv.out10u$lambda.min
bestlam10u

# RMSE associated with the value of lambda chosen by cross-validation >> USER GRID
ridge.pred <- predict(ridge.mod, s = bestlam10u, newx = x[tst,])
MSE(ridge.pred, y[tst])


# Check LOOCV
cv.outl = cv.glmnet(x[tr,], y[tr], alpha = 0, lambda=grid, nfolds=300)
plot(cv.outl$lambda,cv.outl$cvm, main="LOOCV user-defined grid", xlab="Lambda", ylab="CV MSE")
bestlaml = cv.outl$lambda.min
bestlaml

# RMSE associated with the value of lambda chosen by cross-validation >> LOOCV USER GRID
ridge.pred <- predict(ridge.mod, s = bestlaml, newx = x[tst,])
MSE(ridge.pred, y[tst])


# define a finer grid just to check values closer to 0 (100 equally spaced point from 0.001 to 5)
grid2 = seq(0.001, 5, length=100) ## change 

# Re-Check LOOCV with the new grid
cv.outl = cv.glmnet(x[tr,], y[tr], alpha = 0, lambda=grid2, nfolds=300)
plot(cv.outl$lambda,cv.outl$cvm, main="LOOCV user-defined grid", xlab="Lambda", ylab="CV MSE")
bestlaml = cv.outl$lambda.min
bestlaml

# RMSE associated with the value of lambda chosen by cross-validation >> LOOCV with NEW grid
ridge.pred <- predict(ridge.mod, s = bestlaml, newx = x[tst,])
MSE(ridge.pred, y[tst])


#################################
### Ridge-PCA Hybrid Approach ###
#################################
# standardise before PCA
pca_model <- preProcess(x[tr,], method = c("center", "scale", "pca"))

x_train_pca <- predict(pca_model, x[tr,])
x_test_pca <- predict(pca_model, x[tst,])

grid = 10^seq(10, -2, length = 100) ## define again further on
ridge_mod <- glmnet(as.matrix(x_train_pca), y[tr], alpha = 0, lambda = grid, thresh = 1e-12)

# choose lambda based on 10 fold CV
cv_out = cv.glmnet(as.matrix(x_train_pca), y[tr], alpha = 0, lambda = grid, nfolds = 10)
best_lambda = cv_out$lambda.min

plot(cv_out)
plot(cv_out$lambda, cv_out$cvm, main = "10-fold CV Ridge-PCA", xlab = "Lambda", ylab = "CV MSE")

ridge_pred <- predict(ridge_mod, s = best_lambda, newx = as.matrix(x_test_pca))

MSE_ridge_pca <- MSE(ridge_pred, y[tst])


########################
### Regression Tree  ###
########################

# Grow a large tree
big_tree <- tree(formula, data=train, mindev=0.0001) ## change mindev for depth of tree

# Check the number of leaves in the big tree
length(unique(big_tree$where))

# Perform cross-validation to find the optimal tree size
cv_result <- cv.tree(big_tree, FUN = prune.tree)

# Plot the cross-validated error as a function of the tree size (number of terminal nodes)
plot(cv_result$size, cv_result$dev, type = "b", col = "blue", 
     xlab = "Number of Terminal Nodes", ylab = "Cross-Validated Deviance", 
     main = "Cross-Validation for Tree Pruning")

# Find the optimal tree size based on the minimum cross-validated error
optimal_size <- cv_result$size[which.min(cv_result$dev)]


# Prune the tree to the optimal size
pruned_tree <- prune.tree(big_tree, best=optimal_size) ## change best, see the shape of tree first to prune

# Check the pruned tree size
length(unique(pruned_tree$where))

# Plot the big tree and the pruned tree
#windows()  if can
#dev.new() 
par(mfrow=c(1,2))  # Create a 1x2 grid for plotting

# Plot the big tree
plot(big_tree, type="uniform")
text(big_tree, col="blue", label=c("yval"), cex=.8)

# Plot the pruned tree
plot(pruned_tree, type="uniform")
text(pruned_tree, col="blue", label=c("yval"), cex=.8)

# Predict on the test set and compute Mean Squared Error (MSE)
pred_big_tree <- predict(big_tree, newdata=test)
pred_pruned_tree <- predict(pruned_tree, newdata=test)

mse <- function(pred, truth) {
  mean((truth - pred)^2)
}

mse_big_tree <- mse(pred_big_tree, test[[target_column]])
mse_pruned_tree <- mse(pred_pruned_tree, test[[target_column]])



## Saving Dataframe for variables selected by Best Model > Backward Selection:
selected_columns <- c(varselcv[-1])  # Combine the target variable with selected variables
variables_selected <- df[, selected_columns]  # Select the relevant columns
write.csv(variables_selected, file = "variables_selected.csv", row.names = FALSE)


