#Here we follow Medeiros et al. (2019) in defining two functions:

#One for forming forecasts using AR(p) model selected on BIC, which will be called
#on each iteration of the rolling window forecasting exercise.

#The other one for producing the series of h-step forecasts using rolling window.

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#3) lag - the forecast horizon

#4) type -  "fixed" will use the AR(4) in all cases; "bic" will select lags using BIC

runAR=function(Y,indice,lag,type="fixed"){
  
  dum=Y[,ncol(Y)] # extract dummy from data
  Y=Y[,-ncol(Y)] #data without the dummy
  
  Y2=cbind(Y[,indice]) #Y variable (CPI if indice=1, or PCE-based inflation indice=2)
  aux=embed(Y2,4+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,1] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y2)*lag))] # lags of Y (predictors) corresponding to forecast horizon   
  
  if(lag==1){ 
    X.out=tail(aux,1)[1:ncol(X)] #retrieve 4 last observations if one-step forecast 
  }else{
    X.out=aux[,-c(1:(ncol(Y2)*(lag-1)))] #delete first (h-1) columns of aux,  
    X.out=tail(X.out,1)[1:ncol(X)] #last observations: y_T,y_t-1...y_t-h
  }
  
  dum=tail(dum,length(y)) #cut the dummy to size to account for lost observations due to lags
  
  if(type=="fixed"){ #if fixed at AR(4)
    model=lm(y~X+dum) #estimate direct h-step AR(4) by OLS with the dummy
    coef=coef(model)[1:(ncol(X)+1)] #extract coefficients
  }
  
  if(type=="bic"){ #if selection on BIC
    bb=Inf #initialize the "best BIC" at a huge number
    for(i in seq(1,ncol(X),1)){ #try for every lag length 1:4
      m=lm(y~X[,1:i]+dum) #estimate AR(i) by OLS with the dummy
      crit=BIC(m) #retrieve BIC
      if(crit<bb){ #if BIC improved
        bb=crit #updated the new best BIC
        model=m #save the model object
        ar.coef=coef(model) #save coefficients
        ar.coef=ar.coef[-length(ar.coef)] #remove the dummy coefficient
      }
    }
    coef=rep(0,ncol(X)+1) #blank vector for coefficients
    coef[1:length(ar.coef)]=ar.coef #fill in the model coefficients (zeros may remain if model less than order 4)
  }
  pred=c(1,X.out)%*%coef #make a forecast using the last few observations: a direct h-step forecast
  
  return(list("model"=model,"pred"=pred,"coef"=coef)) #save estimated AR regression, prediction, and estimated coefficients
}


#This function will repeatedly call the previous function in the rolling window h-step forecasting

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) nprev - number of out-of-sample observations (at the end of the sample)

#3) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#4) lag - the forecast horizon

#5) type -  "fixed" will use the AR(4) in all cases; "bic" will select lags using BIC


ar.rolling.window=function(Y,nprev,indice=1,lag=1,type="fixed"){
  
  save.coef=matrix(NA,nprev,5) #blank matrix for coefficients at each iteration
  save.pred=matrix(NA,nprev,1) #blank for forecasts
  for(i in nprev:1){  #NB: backwards FOR loop: going from 180 down to 1
    
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 491, then 2 to 492 etc.)
    fact=runAR(Y.window,indice,lag) #call the function to fit the AR(p) selected on BIC and generate h-step forecast
    save.coef[(1+nprev-i),]=fact$coef #save estimated coefficients
    save.pred[(1+nprev-i),]=fact$pred #save the forecast
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  
  #Some helpful stuff:
  real=Y[,indice] #get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  
  return(list("pred"=save.pred,"coef"=save.coef,"errors"=errors)) #return forecasts, history of estimated coefficients, and RMSE and MAE for the period.
}

