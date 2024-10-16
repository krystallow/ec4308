library(glmnet)
library(hdm)
library(pls)

set.seed(2457829)

ntrain= #number of training data
  train = sample(1:nrow(x),ntrain)
test = (-train)
MSE <- function(pred, truth){
  return(mean((truth - pred)^2)) 
}

#model.matrix creates a matrix of predictors and transforms qualitative variables into dummy variables which is necessary for glmnet()
x = model.matrix(, data = xx) #placeholder data

x=x[,2:ncol(x)] #omit the intercept -run on demeaned data without one

y =  #define Y
  
################################
# LASSO
################################

grid = 10^seq(10, -2, length = 500) # we can just remove this to use their default grid, this is a logarithmic grid
lasso.mod <- glmnet(x[train,], y[train], alpha = 1, lambda = grid)
plot(lasso.mod) #this allows you to plot coefficients as a function of your L1-budget (and hence lambda)

#This is may be useful to see how and where shrinkage kicks in.
cv.out10u = cv.glmnet(x[train,], y[train], alpha = 1, lambda=grid)
plot(cv.out10u)
plot(cv.out10u$lambda,cv.out10u$cvm, main="10-fold CV user-defined grid", xlab="Lambda", ylab="CV MSE")

bestlam10u = cv.out10u$lambda.min
lasso.pred <- predict(lasso.mod, s = bestlam10u, newx = x[test,])
mse_lasso_cv = MSE(lasso.pred, y[test])

# Plug in lambda choice
rlasso.fit = rlasso(y[train]~x[train,],  post=FALSE)
yhat.rlasso<- predict(rlasso.fit, newdata=x[test,])
mse_plugin_LASSO = MSE(yhat.rlasso, y[test])

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
  elnet.mod <- glmnet(x[train,], y[train], alpha = alpha_val, lambda = grid)
  cv.out10el = cv.glmnet(x[train,], y[train], alpha = alpha_val, lambda = grid)
  
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

#################################
# Principal components regression
##################################

#Examine the fit syntax:
pcr.fit=pcr(y~.,data=, subset=train, scale=TRUE, validation="CV")

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
comps.test <- predict(pcr.pred, newdata = x[test, ]) 
yhat.lasso_pcr <- predict(lasso.mod, newx = comps.test)
mse_lasso_pcr_cv = MSE(yhat.lasso_pcr, y[test])
