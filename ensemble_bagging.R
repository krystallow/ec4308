###################################################
### Ensemble - Boosting for Variable Selection  ###
###################################################

set.seed(2457829) #seed of the random number generator for replicability

#Here we boost the trees and get the best number on OOB (because cv option doesn't work)

#Model and data definition proceed same as in trees, randomForest and lm function
#Options of note here: 
#1) distribution - this will determine what kind of loss minimization problem we solve.
#For us, the relevant choices are 'gaussian' (regression problem - MSE), 'bernoulli' (classification problem).
#There are some more specialized options, e.g. 'adaboost' for classification with an older algorithm, 'laplace' for regression with absolute loss.

#2) bag.fraction: if set to >0, gbm computes the OOB error. Recommended value: 0.5

#3) interaction.depth: depth of each tree - keep it short!

#4) n.trees: number of boosting iterations

#5) shrinkage: the penalty parameter ("learning rate"). Recommended small, typically 0.01 or so. Note reducing this will necessitate increase in n.trees 
#to keep fitting well.
boostfit = gbm(formula,data=train,distribution='gaussian',bag.fraction = .5,
               interaction.depth=2,n.trees=10000,shrinkage=.01)

#We can then compute the estimated optimal M (number of iterations) using the gbm.perf() function:
#Note the option method: choices are "OOB" for out-of-bag or "cv" for cross-validation (need to specify cv.folds in the gbm call)
best = gbm.perf(boostfit, method="OOB")

#NB: Iterations number chosen with OOB tends to be conservative, i.e., a bit less than the number that
#truly minimizes test set error. So expect K-fold CV to produce larger M than the OOB method. 
#Also, using a bit more iterations than the OOB selection would usually improve performance.

#Examine boosting for various number of trees:
set.seed(2457829) #seed



ntreev = c(5,best,10000) # different numbers of iterations considered

nset = length(ntreev) #no of iterations considered

fmat = matrix(0,ntrain,nset) #blank for predictions

#Get the three boosting fits in a loop

for(i in 1:nset) {
  
  boostfit = gbm(formula,data=train,distribution='gaussian',
                 interaction.depth=2,n.trees=ntreev[i],shrinkage=.1)
  
  #NB: here I use a larger  shrinkage parameter (0.1 vs 0.01) in order to have a more illustrative figure
  
  fmat[,i] = predict(boostfit,n.trees=ntreev[i])
}


#Plot the fits:

windows()
par(mfrow=c(1,3))
oo = order(train$SoldAge)
for(i in 1:nset) {
  plot(train$SoldAge,train$Price,xlab='SoldAge',ylab='Price')
  lines(train$SoldAge[oo],fmat[oo,i],col=i+1,lwd=3,lty=1)
  title(main=paste('boosting, ntree= ',ntreev[i]))
}