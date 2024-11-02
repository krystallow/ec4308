#Load the libraries:
library(glmnet)
library(gbm)
library(randomForest)
library(hdm)
library(rpart)
library(tictoc)
#install.packages("lsei") #package for finding weights using quadratic programming
library(lsei)

##############################
#Preliminary data manipulation
##############################

set.seed(2457829)
library(readxl)
data <- read_excel("Final_Transformed_Dataset.xlsx", range = "B1:DK775")
data = data[-1,]
data_clean <- na.omit(data)
# Now create model matrix
x = model.matrix(INDPRO~., data = data_clean)
y = data_clean$INDPRO
ntrain= floor(0.7* 773)  #number of training data -> 70% of total obs
train = sample(1:nrow(x),ntrain)
training_set = data_clean[train,]
testing_set = data_clean[-train,]
vald=testing_set[1:130,] #reserve the first 250 in the test data as validation set for ensemble weights
testing_set=testing_set[130:nrow(testing_set),] #use the rest as the test data

#Create model matrices, omitting the constant (-1) for the different data sets (in doing this, we expand some categorical predictors into dummies):
x1 = model.matrix(INDPRO ~ ., data = training_set) #training data
x2 = model.matrix(INDPRO ~ ., data = vald) #validation data
x3 = model.matrix(INDPRO ~ ., data = testing_set) #test data
y1=training_set$INDPRO #y for training data
y2 = vald$INDPRO
y3=testing_set$INDPRO  # y for test data

#A bit more advanced MSE function: computes both MSE and its standard error via linear regression of squared errors on a constant
MSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(summary(lm((truth-pred)^2~1))$coef[1:2]) #end function with a return(output) statement. Here we can go straight to return because the object of interest is a simple function of inputs
}


####################################
##LASSO/Ridge/Elastic Net results
##################################

################################
# LASSO
################################
lasso.mod <- glmnet(x[train,], y[train], alpha = 1)

#This is may be useful to see how and where shrinkage kicks in.
cv.out10l = cv.glmnet(x[train,], y[train], alpha = 1)
plot(cv.out10l)
plot(cv.out10l$lambda,cv.out10l$cvm, main="10-fold CV", xlab="Lambda", ylab="CV MSE")
bestlam10l = cv.out10l$lambda.min
lasso.pred1 <- predict(lasso.mod, s = bestlam10l, newx = x3)
mse_lasso_cv = MSE(lasso.pred1, y3)

lasso.selected <- glmnet(x[train, ], y1, alpha = 1, lambda = bestlam10u)

# Plug in lambda choice
rlasso.fit = rlasso(y[train]~x[train,],  post=FALSE)
rlasso.pred <- predict(rlasso.fit, x3)

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
rlassop.fit = rlasso(y[train]~x[train,],  post=TRUE)
yhat.rlassop<- predict(rlassop.fit, newdata=x3)
################################
# Elastic Net
################################

elnet.mod <- glmnet(x[train,], y[train], alpha = 0.5)
cv.out10el = cv.glmnet(x[train,], y[train], alpha = 0.5)
bestlam10e = cv.out10el$lambda.min
elnet.pred1 <- predict(elnet.mod, s = bestlam10e, newx = x3)

#Predict the test data

tic()
rftune=tuneRF(x1, y1, mtryStart=floor(sqrt(ncol(x1))), stepFactor=2, improve=0.05, nodesize=5, ntree=5000, doBest=TRUE, plot=FALSE, trace=FALSE)
toc()
rftune$mtry #tunde number of predictors

rft.pred = predict(rftune, newdata=x3)
MSE(rft.pred,y3)

############################
#Tree-based method results
############################

#Fit a single pruned tree first using default settings:
big.tree = rpart(INDPRO~.,method="anova",data=training_set) #big tree

bestcp=big.tree$cptable[which.min(big.tree$cptable[,"xerror"]),"CP"] #get best penalty

best.tree = prune(big.tree,cp=bestcp) #get tree for best cp on CV 

tree.pred = predict(best.tree, newdata=testing_set) #predict on test data
MSE(tree.pred,y3)

#Fit boosted tree, with d=5 (takes a long time!!):
tic()
boost.fit = gbm(INDPRO~.,data=training_set,distribution='gaussian',interaction.depth=5,n.trees=10000,shrinkage=.01,cv.folds=10)
bestd5cv=gbm.perf(boost.fit, method="cv")
boost.pred = predict(boost.fit,newdata = testing_set,n.trees = bestd5cv)
toc()
MSE(boost.pred, y3)

#Fit boosted tree, with d=2 (takes a long time!!):
tic()
boost.fit2 = gbm(INDPRO~.,data=training_set,distribution='gaussian',interaction.depth=2,n.trees=10000,shrinkage=.01,cv.folds=10)
bestd5cv2=gbm.perf(boost.fit2, method="cv")
boost.pred2 = predict(boost.fit2,newdata = testing_set,n.trees = bestd5cv2)
toc()
MSE(boost.pred2, y3)


#Fit a random forest using the tuneRF() function to pick the optimal number of predictors (attention - takes a long time!):
#(takes a long time!)
tic()
rftune=tuneRF(x1, y1, mtryStart=floor(sqrt(ncol(x1))), stepFactor=2, improve=0.05, nodesize=5, ntree=5000, doBest=TRUE, plot=FALSE, trace=FALSE)
toc()
rftune$mtry #tunde number of predictors

rft.pred = predict(rftune, newdata=x3)
MSE(rft.pred,y3)

## Ridge
## Choose Lambda by 10-fold CV with grid:
grid = 10^seq(10, -2, length = 100)
if (all(train <= nrow(x))) {  
  # Fit the ridge regression model
  ridge.mod <- glmnet(x[train, ], y[train], alpha = 0, lambda = grid, thresh = 1e-12)
} else {  
  stop("Some indices in tr are out of bounds for the model matrix x.")
}

cv.out10r = cv.glmnet(x[train,], y[train], alpha = 0) 
bestlam10r = cv.out10r$lambda.min
ridge.pred1 <- predict(ridge.mod, s = bestlam10r, newx = x3)

############################################
####Hybrid method: alternating LASSO and RF
############################################

#Initial step: we already fitted rLASSO - take its residuals:

rhat.rlasso=rlasso.fit$residuals

#Second stage: fit random forest to rLASSO residuals 
#Using tuneRF here as well:
#(NB: takes a long time!)
tic()
rffit.ht = tuneRF(x1, rhat.rlasso, mtryStart=floor(sqrt(ncol(x1))), stepFactor=2, improve=0.05, nodesize=5, ntree=5000, doBest=TRUE, plot=FALSE, trace=FALSE)
toc()
yhat.rlasso<- predict(rlasso.fit, newdata=x3)
rfht.pred = predict(rffit.ht, newdata=x3)
hybridt=rfht.pred+yhat.rlasso
mse_lasso_rf = MSE(hybridt, y3)
#######################################
###Forecast combinations
#######################################

#In order to get weights, we need to fit all the models on the validation data first:
lasso.predc <- predict(lasso.mod, s = bestlam10l, newx = x2)
yhat.rlassoc<- predict(rlasso.fit, newdata=x2)
yhat.rlassopc<- predict(rlassop.fit, newdata=x2)
ridge.predc <- predict(ridge.mod, s = bestlam10r, newx = x2)
elnet.predc <- predict(elnet.mod, s = bestlam10e, newx = x2)
rft.predc=predict(rftune, newdata=x2)
boost.pred2c = predict(boost.fit2,newdata = vald,n.trees = bestd5cv2)
boost.predc = predict(boost.fit,newdata = vald,n.trees = bestd5cv)
rfht.predc = predict(rffit.ht, newdata=x2)
hybridtc=rfht.predc+yhat.rlassoc


#Get weights

#First, we form a design matrix where the X variables are the assorted forecasts on the validation set:
fmatu=cbind(lasso.predc,yhat.rlassoc,yhat.rlassopc,ridge.predc,elnet.predc,rft.predc,boost.pred2c,boost.predc,hybridtc)

#GR weights, no constant, all restrictions in place
gru=lsei(fmatu, y2, c=rep(1,9), d=1, e=diag(9), f=rep(0,9))

View(gru) #Examine weights

#Combine the forecasts with nonzero weights:
combpredu=gru[6]*rft.pred+gru[8]*boost.pred+gru[9]*as.matrix(hybridt)

mse_no_constant_all_restrictions = MSE(y3,combpredu) #check MSE

#Redefine the X matrix for forecasts by adding a column of ones - a constant regressor
fmatb=cbind(rep(1,nrow(lasso.predc)),lasso.predc,yhat.rlassoc,yhat.rlassopc,ridge.predc,elnet.predc,rft.predc,boost.pred2c,boost.predc,hybridtc)

#From auxiliary matrix that is identity, but the (1,1) element is 0 (so that the constant regressor is not constrained to sum to one with the other weights)
temp=diag(10)
temp[1,1]=0

#Find the GR weights under constraints, but with constant in the regression:
grb=lsei(fmatb, y2, c=c(0,rep(1,9)), d=1, e=temp, f=rep(0,10))

#From the forecasts using nonzero weights:
combpredb=grb[1]+grb[7]*rft.pred+grb[9]*boost.pred+grb[10]*as.matrix(hybridt)
mse_nonzero_weights_constant = MSE(y3,combpredb) #check MSE


#unrestricted weights: no constraints, no constant
grunr=lsei(fmatu, y2)

#Form combined forecast (almost all weights nonzero, so do vector product to sum):
combpredur=cbind(lasso.pred1,yhat.rlasso,yhat.rlassop,ridge.pred1,elnet.pred1,rft.pred,boost.pred2,boost.pred,hybridt)%*%grunr
mse_no_constraints_no_constant = MSE(y3,combpredur)

#unrestricted weights: no constraints, but include constant
grunrc=lsei(fmatb, y2)
combpredurc=cbind(rep(1,nrow(lasso.pred1)),lasso.pred1,yhat.rlasso,yhat.rlassop,ridge.pred1,elnet.pred1,rft.pred,boost.pred2,boost.pred,hybridt)%*%grunrc
mse_no_contraints_constant = MSE(y3,combpredurc)

#LASSO weights, no constant:
rlassocomb.fit = rlasso(y2~fmatu,  post=FALSE,intercept=FALSE)
#Form combination:
comblasso=cbind(lasso.pred1,yhat.rlasso,yhat.rlassop,ridge.pred1,elnet.pred1,rft.pred,boost.pred2,boost.pred,hybridt)%*%rlassocomb.fit$coefficients
mse_lasso_weights_no_constant = MSE(y3,comblasso)

#LASSO weights, include constant:
rlassocombc.fit = rlasso(y2~fmatu,  post=FALSE,intercept=TRUE)
comblassoc=cbind(rep(1,nrow(lasso.pred1)),lasso.pred1,yhat.rlasso,yhat.rlassop,ridge.pred1,elnet.pred1,rft.pred,boost.pred2,boost.pred,hybridt)%*%rlassocombc.fit$coefficients
mse_lasso_weights_constant = MSE(y3,comblassoc)

mse_results <- data.frame(
  Model = character(),
  RMSE = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

mse_results <- rbind(mse_results, 
                     data.frame(Model = "Hybrid (rLASSO + RF", RMSE = sqrt(mse_lasso_rf[1]), SE = mse_lasso_rf[2]),
                     data.frame(Model = "GR, no constant", RMSE = sqrt(mse_no_constant_all_restrictions[1]), SE = mse_no_constant_all_restrictions[2]),
                     data.frame(Model = "GR, constant", RMSE = sqrt(mse_nonzero_weights_constant[1]), SE = mse_nonzero_weights_constant[2]),
                     data.frame(Model = "GR, unrestricted, no constant", RMSE = sqrt(mse_no_constraints_no_constant[1]),SE = mse_no_constraints_no_constant[2]),
                     data.frame(Model = "GR, unrestricted, constant", RMSE = sqrt(mse_no_contraints_constant[1]), SE = mse_no_contraints_constant[2]),
                     data.frame(Model = "LASSO, no constant", RMSE = sqrt(mse_lasso_weights_no_constant[1]), SE = mse_lasso_weights_no_constant[2]),
                     data.frame(Model = "LASSO, constant", RMSE = sqrt(mse_lasso_weights_constant[1]), SE = mse_lasso_weights_constant[2]))

mse_results
weights = list(gru,grb,grunr,grunrc,rlassocomb.fit$coefficients,rlassocombc.fit$coefficients)
#Save results as a lot of heavy computation was done:
save.image("HybridModels.RData")