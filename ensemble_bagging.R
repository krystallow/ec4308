###################################################
### Ensemble - Boosting for Variable Selection  ###
###################################################
library(gbm)
library(readxl)

set.seed(2457829) #seed of the random number generator for replicability

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]


ntrain <- floor(0.7 * nrow(df))  # take 70% train
tr = sample(1:nrow(df),ntrain)  
train = df[tr,]   
test = df[-tr,]  

boostfit = gbm(INDPRO ~.,data=train,distribution='gaussian',bag.fraction = .5,
               interaction.depth=2,n.trees=10000,shrinkage=.01)

#We can then compute the estimated optimal M (number of iterations) using the gbm.perf() function:
#Note the option method: choices are "OOB" for out-of-bag or "cv" for cross-validation (need to specify cv.folds in the gbm call)
best = gbm.perf(boostfit, method="OOB")

ntreev = c(5,best,10000) # different numbers of iterations considered

nset = length(ntreev) #no of iterations considered

fmat = matrix(0,ntrain,nset) #blank for predictions

#Get the three boosting fits in a loop

for(i in 1:nset) {
  
  boostfit = gbm(INDPRO ~.,data=train,distribution='gaussian',
                 interaction.depth=2,n.trees=ntreev[i],shrinkage=.1)
  
  #NB: here I use a larger  shrinkage parameter (0.1 vs 0.01) in order to have a more illustrative figure
  
  fmat[,i] = predict(boostfit,n.trees=ntreev[i])
}

boosting_pred <- predict(boostfit)
#Plot the fits:
windows()
par(mfrow=c(1,3))
oo = order(train$FEDFUNDS)
for(i in 1:nset) {
  plot(train$FEDFUNDS,train$INDPRO,xlab='FEDFUNDS',ylab='INDPRO')
  lines(train$FEDFUNDS[oo],fmat[oo,i],col=i+1,lwd=3,lty=1)
  title(main=paste('boosting, ntree= ',ntreev[i]))
}

par(mfrow=c(1,3))
oo = order(train$IPMAT)
for(i in 1:nset) {
  plot(train$IPMAT,train$INDPRO,xlab='IPMAT',ylab='INDPRO')
  lines(train$IPMAT[oo],fmat[oo,i],col=i+1,lwd=3,lty=1)
  title(main=paste('boosting, ntree= ',ntreev[i]))
}

boosting_pred <- predict(boostfit)

###################################################
### Ensemble - Boosting for Variable Selection  ### > ON CV
###################################################
boostfit_cv <- gbm(
  INDPRO ~.,
  data = train,
  distribution = 'gaussian',
  bag.fraction = 0.5,
  interaction.depth = 2,
  n.trees = 10000,
  shrinkage = 0.01,
  cv.folds = 5
)


best_cv <- gbm.perf(boostfit_cv, method = "cv")

set.seed(2457829) 

ntree_cv = c(5,best_cv,10000) # different numbers of iterations considered

nset_cv = length(ntree_cv) #no of iterations considered

fmat = matrix(0,ntrain,nset_cv) #blank for predictions

#Get the three boosting fits in a loop

for(i in 1:nset_cv) {
  boostfit_cv <- gbm(
    INDPRO ~.,
    data = train,
    distribution = 'gaussian',
    bag.fraction = 0.5,
    interaction.depth = 2,
    n.trees = 10000,
    shrinkage = 0.01,
    cv.folds = 5
  )
  fmat[,i] = predict(boostfit_cv,n.trees=ntree_cv[i])
}

boostingcv_pred <- predict(boostfit_cv)

#Plot the fits:
windows()
par(mfrow=c(1,3))
oo = order(train$FEDFUNDS)
for(i in 1:nset) {
  plot(train$FEDFUNDS,train$INDPRO,xlab='FEDFUNDS',ylab='INDPRO')
  lines(train$FEDFUNDS[oo],fmat[oo,i],col=i+1,lwd=3,lty=1)
  title(main=paste('boosting, ntree= ',ntreev[i]))
}

par(mfrow=c(1,3))
oo = order(train$IPMAT)
for(i in 1:nset) {
  plot(train$IPMAT,train$INDPRO,xlab='IPMAT',ylab='INDPRO')
  lines(train$IPMAT[oo],fmat[oo,i],col=i+1,lwd=3,lty=1)
  title(main=paste('boosting, ntree= ',ntreev[i]))
}

