#Here we follow Medeiros et al. (2019) in defining four functions:

#One for forming forecasts using LASSO or Elastic Net model selected on BIC, which will be called
#on each iteration of the rolling window forecasting exercise.

#The second one for producing the series of h-step LASSO forecasts using rolling window.

#The third one will take the resulting coefficients from LASSO, seek out the nonzero ones and rerun OLS with selected predictors to produce post-LASSO forecast.
#The fourth one will call on the third one for producing the series of h-step post-LASSO forecasts using rolling window.

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#3) lag - the forecast horizon

#4) pl - set to TRUE if wnat to do post-LASSO or FALSE for vanilla LASSO.
runrlasso=function(Y,indice,lag,pl){
  
  dum=Y[,ncol(Y)] # extract dummy from data
  Y=Y[,-ncol(Y)] #data without the dummy
  comp=princomp(scale(Y,scale=FALSE)) # compute principal components to add as predictors
  Y2=cbind(Y,comp$scores[,1:6]) #augment predictors by the first 6 principal components
  aux=embed(Y2,1+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))]   # lags of Y (predictors) corresponding to forecast horizon   
  
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)] #retrieve the last  observations if one-step forecast  
  }else{
    X.out=aux[,-c(1:(ncol(Y2)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)] #last observations: y_T,y_t-1...y_t-h
  }
  dum=tail(dum,length(y)) #cut the dummy to size to account for lost observations due to lags
  
  
  #Here we replace the glmnet-based LASSO by rlasso:
  
  #If we want to use predict() with rlasso, we need to stuff things into a dataframe
  #type object first:
  
  #All the predictors:
  tempd=cbind(X,dum)
  #Append the test set at the end and set as dataframe:
  df=as.data.frame(rbind(tempd,c(X.out,0)))
  
  #Run (post-)lasso (default setting, heteroskedasticity-adjusted):
  model=rlasso(df[1:nrow(X),],y,post=pl) 
  
  
  pred=predict(model,newdata=df[nrow(X)+1,]) #generate the forecast (note that c(X.out,0) is the last row of the dataframe)
  
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}

runrlasso_all=function(Y,indice,lag,pl){
  
  dum=Y[,ncol(Y)] # extract dummy from data
  Y=Y[,-ncol(Y)] #data without the dummy
  comp=princomp(scale(Y,scale=FALSE)) # compute principal components to add as predictors
  Y2=cbind(Y,comp$scores[,1:35]) #augment predictors by the first 6 principal components
  aux=embed(Y2,1+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))]   # lags of Y (predictors) corresponding to forecast horizon   
  
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)] #retrieve the last  observations if one-step forecast  
  }else{
    X.out=aux[,-c(1:(ncol(Y2)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)] #last observations: y_T,y_t-1...y_t-h
  }
  dum=tail(dum,length(y)) #cut the dummy to size to account for lost observations due to lags
  
  
  #Here we replace the glmnet-based LASSO by rlasso:
  
  #If we want to use predict() with rlasso, we need to stuff things into a dataframe
  #type object first:
  
  #All the predictors:
  tempd=cbind(X,dum)
  #Append the test set at the end and set as dataframe:
  df=as.data.frame(rbind(tempd,c(X.out,0)))
  
  #Run (post-)lasso (default setting, heteroskedasticity-adjusted):
  model=rlasso(df[1:nrow(X),],y,post=pl) 
  
  
  pred=predict(model,newdata=df[nrow(X)+1,]) #generate the forecast (note that c(X.out,0) is the last row of the dataframe)
  
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}
#This function will repeatedly call the previous function in the rolling window h-step forecasting

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) nprev - number of out-of-sample observations (at the end of the sample)

#3) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#4) lag - the forecast horizon

#5) pl - set to TRUE if wnat to do post-LASSO or FALSE for vanilla LASSO.


rlasso.rolling.window=function(Y,nprev,indice=1,lag=1,pl=FALSE){      
  # Use a flag to initialize save.coef and save.pred inside the loop dynamically
  initialized = FALSE
  
  for(i in nprev:1){ # Backwards FOR loop: going from nprev down to 1     
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] # Define the estimation window
    lasso=runrlasso(Y.window, indice, lag, pl) # Call the function to fit the rLASSO model
    
    # Dynamically initialize save.coef and save.pred after the first iteration
    if (!initialized) {
      save.coef=matrix(NA, nprev, length(lasso$model$coef)) # Initialize save.coef based on number of coefficients
      save.pred=matrix(NA, nprev, length(lasso$pred)) # Initialize save.pred based on the number of predictions
      initialized = TRUE
    }
    
    # Save coefficients and predictions
    save.coef[(1+nprev-i),]=lasso$model$coef
    save.pred[(1+nprev-i),]=lasso$pred
    
    # Debugging output
    cat("Iteration:", (1+nprev-i), " | Length of lasso$pred:", length(lasso$pred), "\n")
  }
  
  # Some helpful stuff:   
  real=Y[,indice] # Get actual values   
  plot(real, type="l")   
  
  # Plot predictions (only the first prediction if multiple values)
  lines(c(rep(NA, length(real)-nprev), save.pred[,1]), col="red") # Padded with NA for blanks, plot predictions vs. actual 
  
  # Compute errors (use first prediction if multiple)
  rmse=sqrt(mean((tail(real,nprev)-save.pred[,1])^2)) # Compute RMSE   
  mae=mean(abs(tail(real,nprev)-save.pred[,1])) # Compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse, "mae"=mae) # Stack errors in a vector
  
  return(list("pred"=save.pred, "coef"=save.coef, "errors"=errors)) # Return forecasts, history of estimated coefficients, and RMSE and MAE for the period.
}


rlasso.rolling.window.all=function(Y,nprev,indice=1,lag=1,pl=FALSE){      
  # Use a flag to initialize save.coef and save.pred inside the loop dynamically
  initialized = FALSE
  
  for(i in nprev:1){ # Backwards FOR loop: going from nprev down to 1     
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] # Define the estimation window
    lasso=runrlasso_all(Y.window, indice, lag, pl) # Call the function to fit the rLASSO model
    
    # Dynamically initialize save.coef and save.pred after the first iteration
    if (!initialized) {
      save.coef=matrix(NA, nprev, length(lasso$model$coef)) # Initialize save.coef based on number of coefficients
      save.pred=matrix(NA, nprev, length(lasso$pred)) # Initialize save.pred based on the number of predictions
      initialized = TRUE
    }
    
    # Save coefficients and predictions
    save.coef[(1+nprev-i),]=lasso$model$coef
    save.pred[(1+nprev-i),]=lasso$pred
    
    # Debugging output
    cat("Iteration:", (1+nprev-i), " | Length of lasso$pred:", length(lasso$pred), "\n")
  }
  
  # Some helpful stuff:   
  real=Y[,indice] # Get actual values   
  plot(real, type="l")   
  
  # Plot predictions (only the first prediction if multiple values)
  lines(c(rep(NA, length(real)-nprev), save.pred[,1]), col="red") # Padded with NA for blanks, plot predictions vs. actual 
  
  # Compute errors (use first prediction if multiple)
  rmse=sqrt(mean((tail(real,nprev)-save.pred[,1])^2)) # Compute RMSE   
  mae=mean(abs(tail(real,nprev)-save.pred[,1])) # Compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse, "mae"=mae) # Stack errors in a vector
  
  return(list("pred"=save.pred, "coef"=save.coef, "errors"=errors)) # Return forecasts, history of estimated coefficients, and RMSE and MAE for the period.
}
