load("boostingtree.RData")
load("lasso_rlasso_postlasso_elasticnet.RData")
load("HybridModels.RData")
load("rf.RData")

##################################
### Diebold-Mariano (DM) tests ###
##################################

## AR squared loss for different horizons ##
lbench1c = (oosy - bar1c$pred)^2
lbench3c = (oosy - bar3c$pred)^2
lbench6c = (oosy - bar6c$pred)^2
lbench12c = (oosy - bar12c$pred)^2

## AR squared loss for with test set of 101 for different horizons ##
lbench1c_101 = (y3 - bar1c_101$pred)^2
lbench3c_101 = (y3 - bar3c_101$pred)^2
lbench6c_101 = (y3 - bar6c_101$pred)^2
lbench12c_101 = (y3 - bar12c_101$pred)^2

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
lrf1c=(oosy-rf1c$pred)^2
lrf3c=(oosy-rf3c$pred)^2
lrf6c=(oosy-rf6c$pred)^2
lrf12c=(oosy-rf12c$pred)^2

## Hybrid (rLASSO + RF) squared loss
llrf = (y3-hybridt)^2 #y3 (test set for hybrid and combi models) instead of oosy due to length difference

## GR, unrestricted, no constant squared loss
lgr = (y3-combpredur)^2 #y3 (test set for hybrid and combi models) instead of oosy due to length difference

## Loss differentials for RF and AR
drfar1 = lrf1c - lbench1c
drfar3 = lrf1c - lbench3c
drfar6 = lrf1c - lbench6c
drfar12 = lrf1c - lbench12c

# Regress (RF-AR) for 1-step forecasts 
dmrfar1 = lm(drfar1~1)
DM_RandomForest_Autoregression_1_Step = dmrfar1$residuals
acf(DM_RandomForest_Autoregression_1_Step)
dmrfar1$coefficients/sqrt(NeweyWest(dmrfar1, lag=5))

# Regress (RF-AR) for 3-step forecasts 
dmrfar3 = lm(drfar3~1)
DM_RandomForest_Autoregression_3_Steps = dmrfar3$residuals
acf(DM_RandomForest_Autoregression_3_Steps)
dmrfar3$coefficients/sqrt(NeweyWest(dmrfar3, lag=5))

# Regress (RF-AR) for 6-step forecasts 
dmrfar6 = lm(drfar6~1)
DM_RandomForest_Autoregression_6_Steps = dmrfar6$residuals
acf(DM_RandomForest_Autoregression_6_Steps)
dmrfar6$coefficients/sqrt(NeweyWest(dmrfar6, lag=5))

# Regress (RF-AR) for 12-step forecasts 
dmrfar12 = lm(drfar12~1)
DM_RandomForest_Autoregression_12_Steps = dmrfar12$residuals
acf(DM_RandomForest_Autoregression_12_Steps)
dmrfar12$coefficients/sqrt(NeweyWest(dmrfar12, lag=5))

## Loss differentials for Elastic Net and AR
delar1 = lelasticnet1c - lbench1c
delar3 = lelasticnet3c - lbench3c
delar6 = lelasticnet6c - lbench6c
delar12 = lelasticnet12c -lbench12c

# Regress (EL-AR) for 1-step forecasts 
dmelar1 = lm(delar1~1)
DM_ElasticNet_Autoregression_1_Step = dmelar1$residuals
acf(DM_ElasticNet_Autoregression_1_Step)
dmelar1$coefficients/sqrt(NeweyWest(dmelar1, lag=5))

# Regress (EL-AR) for 3-step forecasts 
dmelar3 = lm(delar3~1)
DM_ElasticNet_Autoregression_3_Steps = dmelar3$residuals
acf(DM_ElasticNet_Autoregression_3_Steps)
dmelar3$coefficients/sqrt(NeweyWest(dmelar3, lag=5))

# Regress (EL-AR) for 6-step forecasts 
dmelar6 = lm(delar6~1)
DM_ElasticNet_Autoregression_6_Steps = dmelar6$residuals
acf(DM_ElasticNet_Autoregression_6_Steps)
dmelar6$coefficients/sqrt(NeweyWest(dmelar6, lag=5))

# Regress (EL-AR) for 12-step forecasts 
dmelar12 = lm(delar12~1)
DM_ElasticNet_Autoregression_12_Steps = dmelar12$residuals
acf(DM_ElasticNet_Autoregression_12_Steps)
dmelar12$coefficients/sqrt(NeweyWest(dmelar12, lag=5))


## Loss differentials for GR (unrestricted, no constant) and AR
dgrar = lgr - lbench1c_101

# Regress (GR-AR)
dmgrar = lm(dgrar~1)
DM_GR_Autoregression = dmgrar$residuals
acf(DM_GR_Autoregression)
dmgrar$coefficients/sqrt(NeweyWest(dmgrar, lag=5))

## Loss differential for Hybrid (rLASSO + RF) and AR
dhar = llrf -lbench1c_101

# Regress (Hybrid-AR)
dmhar = lm(dhar~1)
DM_Hybrid_Autoregression = dmhar$residuals
acf(DM_Hybrid_Autoregression)
dmhar$coefficients/sqrt(NeweyWest(dmhar, lag=5))

## Loss differentials for GR (unrestricted, no constant) and Hybrid (rLASSO + RF)
dhgr = lgr - llrf

# Regress (Hybrid-GR)
dmhgr = lm(dhgr~1)
DM_Hybrid_GR = dmhgr$residuals
acf(DM_Hybrid_GR)
dmhgr$coefficients/sqrt(NeweyWest(dmhgr, lag=5))

save.image(file = "DM_test.RData")


