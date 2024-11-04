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

boost1c <- list()
boost3c <- list()
boost6c <- list()
boost12 <- list()

# Function for multi-step forecasting with additional outputs
multi_step_forecast <- function(model, train_data, test_data, horizon) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  # Initialise a vector to store predictions for each horizon
  predictions <- numeric(n_test)
  
  # Combine the training and test data for sequential prediction
  combined_data <- train_data
  
  for (i in 1:(n_test - horizon + 1)) {
    combined_data <- rbind(combined_data, test_data[1:(i - 1), , drop = FALSE])
    
    for (h in 1:horizon) {
      current_test_point <- test_data[i + h - 1, , drop = FALSE]
      combined_data <- rbind(combined_data, current_test_point)
      
      # Predict the current step
      predictions[i] <- predict(model, newdata = combined_data[nrow(combined_data), , drop = FALSE], n.trees = best_cv)
      
      # Update INDPRO for step-by-step forecasting
      combined_data[nrow(combined_data), "INDPRO"] <- predictions[i]
    }
  }
  
  # Calculate model coefficients as feature importances
  coefficients <- summary(model, n.trees = best_cv, plotit = FALSE)
  
  # Calculate RMSE for predictions
  actual_values <- test_data$INDPRO
  rmse_value <- sqrt(mean((actual_values - predictions)^2))
  
  return(list(pred = predictions, coef = coefficients, errors = c(rmse = rmse_value)))
}

boost1c <- multi_step_forecast(boostfit_cv, train, test, horizon = 1)
boost3c <- multi_step_forecast(boostfit_cv, train, test, horizon = 3)
boost6c <- multi_step_forecast(boostfit_cv, train, test, horizon = 6)
boost12c <- multi_step_forecast(boostfit_cv, train, test, horizon = 12)


###############################################
### Selected variables on importance score  ###
###############################################
importance_scores <- summary(boostfit_cv, plotit = FALSE)
selected_vars <- importance_scores$var[importance_scores$rel.inf > 1]  # Keep variables with > 1% importance

#######################
### Forecast Plots  ###
#######################
library(gridExtra)

full_indices <- seq_len(nrow(df))       
test_indices <- (ntrain + 1):nrow(df)    

actual_values <- data.frame(Index = test_indices, INDPRO = df[test_indices, "INDPRO"])

# 1-step forecast plot
forecasted_values_1 <- data.frame(Index = test_indices, INDPRO = boost1c$pred)
p1 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_1, aes(x = Index, y = INDPRO), color = "red", size = 1) +
  labs(title = "Forecast Horizon: 1 Step", x = "Observation Index", y = "INDPRO") +
  theme_minimal()

# 3-step forecast plot
forecasted_values_3 <- data.frame(Index = test_indices, INDPRO = boost3c$pred)
p3 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_3, aes(x = Index, y = INDPRO), color = "red", size = 1) +
  labs(title = "Forecast Horizon: 3 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()

# 6-step forecast plot
forecasted_values_6 <- data.frame(Index = test_indices, INDPRO = boost6c$pred)
p6 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_6, aes(x = Index, y = INDPRO), color = "red", size = 1) +
  labs(title = "Forecast Horizon: 6 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()

# 12-step forecast plot
forecasted_values_12 <- data.frame(Index = test_indices, INDPRO = boost12c$pred)
p12 <- ggplot() +
  geom_line(data = actual_values, aes(x = Index, y = INDPRO), color = "black", size = 1) +
  geom_line(data = forecasted_values_12, aes(x = Index, y = INDPRO), color = "red", size = 1) +
  labs(title = "Forecast Horizon: 12 Steps", x = "Observation Index", y = "INDPRO") +
  theme_minimal()

grid.arrange(p1, p3, p6, p12, ncol = 2)



