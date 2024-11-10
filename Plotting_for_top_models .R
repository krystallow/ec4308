load("lasso_rlasso_postlasso_elasticnet.RData")
load("HybridModels.RData")
load("rf.RData")

#Create the time series object collecting 1-step best=performing ML forecasts
ml1.ts=ts(cbind(oosy,bar1c$pred,elasticnet1c$pred,rf1c$pred), start=c(2010,1), end=c(2024,7), freq=12)
plot.ts(ml1.ts[,1], main="1-step ML forecasts", cex.axis=1.5, lwd=2, ylab="Industrial Production")
points(ml1.ts[,2], type="l", col="blue",lwd=2.3)
points(ml1.ts[,3], type="l", col="red",lwd=2.3)
points(ml1.ts[,4], type="l", col="green",lwd=2.3)
legend("bottomleft", c("Industrial Production","AR","Elastic Net","Random Forest"), lty=c(1,1,1,1) ,col=c("black","blue","red","green"))
#Create the time series object collecting 6-step best=performing ML forecasts
ml6.ts=ts(cbind(oosy,bar6c$pred,elasticnet6c$pred,rf6c$pred), start=c(2010,1), end=c(2024,7), freq=12)
plot.ts(ml6.ts[,1], main="6-step ML forecasts", cex.axis=1.5, lwd=2, ylab="Industrial Production")
points(ml6.ts[,2], type="l", col="blue",lwd=2.3)
points(ml6.ts[,3], type="l", col="red",lwd=2.3)
points(ml6.ts[,4], type="l", col="green",lwd=2.3)
legend("bottomleft", c("Industrial Production","AR","Elastic Net","Random Forest"), lty=c(1,1,1,1) ,col=c("black","blue","red","green"))

