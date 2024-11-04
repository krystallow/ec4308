library(glmnet)
library(hdm)
library(pls)
library(readxl)

set.seed(2457829)
library(readxl)
data <- read_excel("Final_Transformed_Dataset.xlsx", range = "B1:DK775")
df <- read_excel("Final_Transformed_Dataset.xlsx", 
                 col_types = c("date", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric"), range = "A1:DK775" )
data = data[-1,]
df = df[-1,]
data_clean <- na.omit(data)
df = na.omit(df)
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


###############
# FORECASTING
###############
library(githubinstall)
githubinstall("HDeconometrics")
library(HDeconometrics)
library(sandwich)
df$sasdate <- as.Date(df$sasdate) 
dum <- ifelse(format(df$sasdate, "%Y") == "2020", 1, 0)
df <- cbind(df, dum = dum)

nprev=230 #number of out-of-sample observations (test window )

oosy=tail(y,nprev) #auxiliary:get the out-of-sample true values (last 180 obs. using tail())

RMSE <- function(pred, truth){
  return(sqrt(mean((truth - pred)^2)))
}


################
# Random Forest
################

library(randomForest)
runrf=function(Y,indice,lag){
  
  dum=Y[,ncol(Y)] # extract dummy from data
  Y=Y[,-ncol(Y)]  #data without the dummy
  comp=princomp(scale(Y,scale=FALSE)) # compute principal components to add as predictors
  Y2=cbind(Y,comp$scores[,1:35]) #augment predictors by the first 35 principal components
  aux=embed(Y2,6+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))]  # lags of Y (predictors) corresponding to forecast horizon 
  
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)]   #retrieve the last  observations if one-step forecast
  }else{
    X.out=aux[,-c(1:(ncol(Y2)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)]  #last observations: y_T,y_t-1...y_t-h
  }
  
  dum=tail(dum,length(y)) #cut the dummy to size to account for lost observations due to lags
  
  model=randomForest(cbind(X,dum),y,importance = TRUE) #fit the random forest on default settings
  pred=predict(model,c(X.out,0)) #generate forecast
  
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}

rf.rolling.window=function(Y,nprev,indice=1,lag=1){
  
  save.importance=list() #blank for saving variable importance
  save.pred=matrix(NA,nprev,1) ##blank for forecasts
  for(i in nprev:1){#NB: backwards FOR loop: going from 180 down to 1
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 491, then 2 to 492 etc.)
    rf=runrf(Y.window,indice,lag)#call the function to fit the Random Forest and generate h-step forecast
    save.pred[(1+nprev-i),]=rf$pred #save the forecast
    save.importance[[i]]=importance(rf$model) #save variable importance
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  #Some helpful stuff:
  real=Y[,indice]#get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  
  return(list("pred"=save.pred,"errors"=errors,"save.importance"=save.importance)) #return forecasts, history of variable importance, and RMSE and MAE for the period.
}


df = df[,-1]
df = as.matrix(df)
rf1c=rf.rolling.window(df,nprev,1,1)
rf3c=rf.rolling.window(df,nprev,1,3)
rf6c=rf.rolling.window(df,nprev,1,6)
rf12c=rf.rolling.window(df,nprev,1,12)
save(rf1c,rf3c,rf6c,rf12c, file = "rf.RData")
par(mfrow = c(2, 2)) # Arrange plots in a 2x2 layout

# Define actual values
real_values <- df[,1] # or whichever variable corresponds to the actual values in df

# Plot for 1-step forecast
plot(real_values, type = "l", main = "RF 1-Step Forecast", ylab = "INDPRO", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf1c$pred)), rf1c$pred), col = "red")

# Plot for 3-step forecast
plot(real_values, type = "l", main = "RF 3-Step Forecast", ylab = "INDPRO", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf3c$pred)), rf3c$pred), col = "red")

# Plot for 6-step forecast
plot(real_values, type = "l", main = "RF 6-Step Forecast", ylab = "INDPRO", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf6c$pred)), rf6c$pred), col = "red")

# Plot for 12-step forecast
plot(real_values, type = "l", main = "RF 12-Step Forecast", ylab = "INDPRO", xlab = "Time")
lines(c(rep(NA, length(real_values) - length(rf12c$pred)), rf12c$pred), col = "red")

# Restore plotting layout
par(mfrow = c(1, 1))


#See the RMSE:
rf.rmse1=rf1c$errors[1]
rf.rmse3=rf3c$errors[1]
rf.rmse6=rf6c$errors[1]
rf.rmse12=rf12c$errors[1]

#Compute squared loss for different horizons (RF)
lrf1c=(oosy-rf1c$pred)^2
lrf3c=(oosy-rf3c$pred)^2
lrf6c=(oosy-rf6c$pred)^2
lrf12c=(oosy-rf12c$pred)^2

