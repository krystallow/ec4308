load("boostingtree.RData")
load("lasso_rlasso_postlasso_elasticnet.RData")


##################################
### Diebold-Mariano (DM) tests ###
##################################

## AR squared loss for different horizons ##

## Boosting squared loss for different horizons ##
lboosting1c <- (oosy - boost1c$pred[1:230])^2
lboosting3c <- (oosy - boost3c$pred[1:230])^2
lboosting6c <- (oosy - boost6c$pred[1:230])^2
lboosting12c <- (oosy - boost12c$pred[1:230])^2

## LASSO squared loss for different horizons ## (rLASSO)
llasso1c <- (oosy - rlasso1c$pred)^2
llasso3c <- (oosy - rlasso3c$pred)^2
llasso6c <- (oosy - rlasso6c$pred)^2
llasso12c <- (oosy - rlasso12c$pred)^2

## post-LASSO squared loss for different horizons ##
lpostlasso1c <- (oosy - prlasso1c$pred)^2
lpostlasso3c <- (oosy - prlasso3c$pred)^2
lpostlasso6c <- (oosy - prlasso6c$pred)^2
lpostlasso12c <- (oosy - prlasso12c$pred)^2

## Elastic Net squared loss for different horizons ##
lelasticnet1c <- (oosy - elasticnet1c$pred)^2
lelasticnet3c <- (oosy - elasticnet3c$pred)^2
lelasticnet6c <- (oosy - elasticnet6c$pred)^2
lelasticnet12c <- (oosy - elasticnet12c$pred)^2


## Random Forest squared loss for different horizons ##




## (Boosting-postLASSO) >> just testing
## Compute loss differentials (d_t) for different horizons ##
dboost_postlasso1 = lboosting1c - lpostlasso1c
dboost_postlasso3 = lboosting3c - lpostlasso3c
dboost_postlasso6 = lboosting6c - lpostlasso6c
dboost_postlasso12 = lboosting12c - lpostlasso12c

#Create ts object containing loss differentials
dtboost_postlasso.ts=ts(cbind(dboost_postlasso1,dboost_postlasso3,dboost_postlasso6,dboost_postlasso12), start=c(1960,3), end=c(2024,7), freq=12)
#Plot them to examine stationarity:
colnames(dtboost_postlasso.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dtboost_postlasso.ts, main="Loss differential Boosting-PostLASSO",cex.axis=1.8)

#Regress d_t (Boosting-postLASSO) for 1-step forecasts on a constant - get estimate of mean(d_t)
#3-step forecast test
dmboost_postlasso1=lm(dboost_postlasso1~1) #regression
acf(dmboost_postlasso1$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmboost_postlasso1$coefficients/sqrt(NeweyWest(dmboost_postlasso1,lag=6)) #form the DM t-statistic


#3-step forecast test
dmboost_postlasso3=lm(dboost_postlasso3~1)
acf(dmboost_postlasso3$residuals)
dmboost_postlasso3$coefficients/sqrt(NeweyWest(dmboost_postlasso3,lag=6))

#6-step forecast test
dmboost_postlasso6=lm(dboost_postlasso6~1)
acf(dmboost_postlasso6$residuals)
dmboost_postlasso6$coefficients/sqrt(NeweyWest(dmboost_postlasso6,lag=6))

#12-step forecast test
dmboost_postlasso12=lm(dboost_postlasso12~1)
acf(dmboost_postlasso12$residuals)
dmboost_postlasso12$coefficients/sqrt(NeweyWest(dmboost_postlasso12,lag=6))

