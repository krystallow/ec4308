############################
### Ensemble - Boosting  ###
############################
library(gbm)
library(readxl)
library(Metrics)

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


multi_step_forecast <- function(model, train_data, test_data, horizon) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  # Initialise a matrix to store predictions for each horizon
  predictions <- numeric(n_test)
  
  combined_data <- train_data  # Initialise with training data
  
  for (i in 1:(n_test - horizon + 1)) {
    # Combine train data with the test data up to the i-th point
    combined_data <- rbind(combined_data, test_data[1:(i - 1), , drop = FALSE])
    
    for (h in 1:horizon) {
      # Add the current test point to combined data for step-by-step forecasting
      current_test_point <- test_data[i + h - 1, , drop = FALSE]
      combined_data <- rbind(combined_data, current_test_point)
      
      # Predict using the model for the current step
      predictions[i] <- predict(model, newdata = combined_data[nrow(combined_data), , drop = FALSE], n.trees = best)
      
      # Use the predicted value to update INDPRO in combined_data
      combined_data[nrow(combined_data), "INDPRO"] <- predictions[i]
    }
  }
  
  return(predictions)
}


h_steps <- c(1, 3, 6, 12)


forecast_results <- list()
rmse_results <- list()

for (h in h_steps) {
  # Generate forecasts for the entire test set for h-step horizon
  forecast_results[[as.character(h)]] <- multi_step_forecast(boostfit, train, test, horizon = h)
  
  # Get the actual values from the test set
  actual_values <- test$INDPRO
  
  # Calculate RMSE for each horizon
  rmse_results[[as.character(h)]] <- rmse(actual_values, forecast_results[[as.character(h)]])
}


print(forecast_results)
print(rmse_results)
