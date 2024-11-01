###################################################
### Ensemble - Boosting for Variable Selection  ### > ON CV
###################################################
library(glmnet)
library(readxl)
library(gbm)
library(Metrics)
library(ggplot2)
library(dplyr)

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, -1]

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
                               "numeric"))

df <- df[-1, ]   

## Create dummy variable for 2020 Covid
df$sasdate <- as.Date(df$sasdate)
dum <- ifelse(format(df$sasdate, "%Y") == "2020", 1, 0)

Y <- cbind(df, dum = dum)
Y <- Y[,-1]
Y <- Y[apply(Y, 1, function(row) all(is.finite(row) & !is.na(row))), ]
df <- Y


ntrain <- floor(0.7 * nrow(df))  # take 70% train
tr = sample(1:nrow(df),ntrain)  
train = df[tr,]   
test = df[-tr,]  

boostfit_cv <- gbm(
  INDPRO ~.,
  data = train,
  distribution = 'gaussian',
  bag.fraction = 0.5,
  interaction.depth = 5,
  n.trees = 10000,
  shrinkage = 0.01,
  cv.folds = 10
)


best_cv <- gbm.perf(boostfit_cv, method = "cv")

set.seed(2457829) 

ntree = c(5,best_cv,10000) # different numbers of iterations considered

nset = length(ntree) #no of iterations considered

fmat = matrix(0,ntrain,nset) #blank for predictions

#Get the three boosting fits in a loop

for(i in 1:nset) {
  boostfit_cv <- gbm(
    INDPRO ~.,
    data = train,
    distribution = 'gaussian',
    bag.fraction = 0.5,
    interaction.depth = 5,
    n.trees = 10000,
    shrinkage = 0.01,
    cv.folds = 10
  )
  fmat[,i] = predict(boostfit_cv,n.trees=ntree[i])
}

boostingcv_pred <- predict(boostfit_cv)

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
      predictions[i] <- predict(model, newdata = combined_data[nrow(combined_data), , drop = FALSE], n.trees = best_cv)
      
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
  forecast_results[[as.character(h)]] <- multi_step_forecast(boostfit_cv, train, test, horizon = h)
  
  # Get the actual values from the test set
  actual_values <- test$INDPRO
  
  # Calculate RMSE for each horizon
  rmse_results[[as.character(h)]] <- rmse(actual_values, forecast_results[[as.character(h)]])
}


print(forecast_results)
print(rmse_results)


###############################################
### Selected variables on importance score  ###
###############################################
importance_scores <- summary(boostfit_cv, plotit = FALSE)
selected_vars <- importance_scores$var[importance_scores$rel.inf > 1]  # Keep variables with > 1% importance


#######################
### Forecast Plots  ###
#######################

full_indices <- seq_len(nrow(df))         
test_indices <- (ntrain + 1):nrow(df)   

actual_values <- data.frame(Index = full_indices, INDPRO = df[, "INDPRO"])

# 1-step forecast plot
forecasted_values_1 <- data.frame(Index = test_indices, INDPRO = forecast_results[["1"]])
p1 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_1, aes(x = Index, y = INDPRO), color = "red", size = 1, linetype = "dashed") +
  labs(title = "Forecast Horizon: 1 Step", x = "Observation Index", y = "INDPRO") +
  theme_minimal()
print(p1)


# 3-step forecast plot
forecasted_values_3 <- data.frame(Index = test_indices, INDPRO = forecast_results[["3"]])
p3 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_3, aes(x = Index, y = INDPRO), color = "red", size = 1, linetype = "dashed") +
  labs(title = "Forecast Horizon: 3 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()
print(p3)

# 6-step forecast plot
forecasted_values_6 <- data.frame(Index = test_indices, INDPRO = forecast_results[["6"]])
p6 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_6, aes(x = Index, y = INDPRO), color = "red", size = 1, linetype = "dashed") +
  labs(title = "Forecast Horizon: 6 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()
print(p6)

# 12-step forecast plot
forecasted_values_12 <- data.frame(Index = test_indices, INDPRO = forecast_results[["12"]])
p12 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_12, aes(x = Index, y = INDPRO), color = "red", size = 1, linetype = "dashed") +
  labs(title = "Forecast Horizon: 12 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()
print(p12)
