################
### ADL Model###
################


############################################
###### Optimal Lags = 1 for all ############ >> No dimension warning >> Overfit
############################################
library(readxl)
library(dynlm)
library(car)

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]   
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]

# Check for multicollinearity and remove problematic variables
vif_values <- vif(lm(INDPRO ~ ., data = df))
print(vif_values)
# Remove variables with high VIF values
df <- df[, !(names(df) %in% c("IPFPNSS", "IPMAT", "IPMANSICS"))]

# Recheck VIF values after removal
vif_values_after <- vif(lm(INDPRO ~ ., data = df))
print(vif_values_after)

y <- df$INDPRO
x <- df[, setdiff(colnames(df), "INDPRO")]


# Define cross-validation function with rolling window
cv_rmse <- function(p, q, y, x_var, data) {
  n <- nrow(data)
  rmse <- numeric() 
  
  for (train_end in (p + max(q) + 1):(n - 1)) {
    train_indices <- 1:train_end
    test_index <- train_end + 1
    
    # Construct the formula as a string
    formula_str <- paste(y, "~ L(", y, ", ", p, ") + L(", x_var, ", ", q, ")", sep = "")
    
    # Fit the model using dynlm
    model <- dynlm(as.formula(formula_str), data = data[train_indices, ])
    
    # Prepare new data for prediction
    newdata <- data[test_index, , drop = FALSE]
    
    # Ensure the newdata has the correct columns
    if (!all(c(y, x_var) %in% colnames(newdata))) {
      stop("newdata does not contain the required columns.")
    }
    
    # Make predictions
    pred <- predict(model, newdata = newdata)
    actual <- data[[y]][test_index]
    
    # Calculate RMSE
    rmse <- c(rmse, sqrt(mean((pred - actual)^2, na.rm = TRUE)))
  }
  
  return(mean(rmse, na.rm = TRUE))
}

# Initialize to store best lag (q) for each independent variable
best_results <- data.frame(variable = character(), best_q = numeric(), best_rmse = numeric(), stringsAsFactors = FALSE)

p_range <- 1:6  # Dependent variable lag
q_range <- 1:4   # Independent variable lags

p_results <- data.frame(p = p_range, rmse = NA)

for (p in p_range) { 
  rmse_value <- cv_rmse(p, 0, "INDPRO", "INDPRO", df)
  p_results$rmse[p_results$p == p] <- rmse_value
}


best_p_index <- which.min(p_results$rmse)
best_p <- p_results$p[best_p_index]
cat("Best lag p for y:", best_p, "\n")

# New function to find RMSE for a single predictor while holding others constant
cv_rmse_single <- function(p, q, y, x_var, data) {
  n <- nrow(data)
  rmse <- numeric() 
  
  for (train_end in (p + q + 1):(n - 1)) {
    train_indices <- 1:train_end
    test_index <- train_end + 1
    
    # Construct the formula for one predictor
    formula_str <- paste(y, "~ L(", y, ", ", p, ") + L(", x_var, ", ", q, ")", sep = "")
    
    # Fit the model
    model <- dynlm(as.formula(formula_str), data = data[train_indices, ])
    
    # Prepare new data for prediction
    newdata <- data[test_index, , drop = FALSE]
    
    # Ensure the newdata has the correct columns
    if (!all(c(y, x_var) %in% colnames(newdata))) {
      stop("newdata does not contain the required columns.")
    }
    
    # Make predictions
    pred <- predict(model, newdata = newdata)
    actual <- data[[y]][test_index]
    
    # Calculate RMSE
    rmse <- c(rmse, sqrt(mean((pred - actual)^2, na.rm = TRUE)))
  }
  
  return(mean(rmse, na.rm = TRUE))
}

# Loop through each independent variable to determine optimal q lag
for (var in colnames(x)) {
  q_results <- data.frame(q = q_range, rmse = NA) 
  
  for (q in q_range) {
    p <- best_p  # evaluate q for the best p
    rmse_value <- cv_rmse_single(p, q, "INDPRO", var, df)  
    q_results$rmse[q_results$q == q] <- rmse_value
  }
  
  best_q_index <- which.min(q_results$rmse)
  best_q <- q_results$q[best_q_index]
  best_rmse <- q_results$rmse[best_q_index]
  
  best_results <- rbind(best_results, data.frame(variable = var, best_q = best_q, best_rmse = best_rmse))
}

print(best_results)

formula_str <- paste("INDPRO ~ L(INDPRO,", best_p, ") + ",  
                     paste(mapply(function(var, q) paste0("L(", var, ",", q, ")"),  
                                  best_results$variable, best_results$best_q), collapse = " + "))

final_model <- dynlm(as.formula(formula_str), data = df)
summary(final_model)




#######################
###### LOOCV - MSE ####
#######################
library(readxl)
library(dynlm)
library(car)

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]   
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]

# Check for multicollinearity and remove problematic variables
#vif_values <- vif(lm(INDPRO ~ ., data = df))
#print(vif_values)
#df <- df[, !(names(df) %in% c("IPFPNSS", "IPMAT", "IPMANSICS"))]
#vif_values_after <- vif(lm(INDPRO ~ ., data = df))
#print(vif_values_after)

y <- "INDPRO"
x <- df[, setdiff(colnames(df), y)]

# Define LOOCV function to calculate MSE
loocv_mse <- function(p, q, y, x_var, data) {
  n <- nrow(data)
  mse_values <- numeric(n)
  
  for (i in 1:n) {
    # Training data excluding the ith observation
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]
    
    # formula for ADL model
    formula_str <- paste(y, "~ L(", y, ", ", p, ") + L(", x_var, ", ", q, ")", sep = "")
    
    model <- tryCatch(dynlm(as.formula(formula_str), data = train_data), error = function(e) NULL)
    
    if (!is.null(model)) {
      # Predict for the left-out observation
      pred <- predict(model, newdata = test_data)
      actual <- test_data[[y]]
      
      # Calculate 
      mse_values[i] <- (pred - actual)^2
    }
  }
  
  return(mean(mse_values, na.rm = TRUE))
}


best_results <- data.frame(variable = character(), best_q = numeric(), best_mse = numeric(), stringsAsFactors = FALSE)

p_range <- 1:12  # Lags for the dependent variable
q_range <- 1:12  # Lags for the independent variables

# optimal lag for dependent variable
p_results <- data.frame(p = p_range, mse = NA)

for (p in p_range) {  
  mse_value <- loocv_mse(p, 0, y, y, df)
  p_results$mse[p_results$p == p] <- mse_value
}

best_p_index <- which.min(p_results$mse)
best_p <- p_results$p[best_p_index]
cat("Best lag p for y:", best_p, "\n")

# to evaluate a single independent variable with a specific lag q
loocv_mse_single <- function(p, q, y, x_var, data) {
  n <- nrow(data)
  mse_values <- numeric(n)
  
  for (i in 1:n) {
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]
    
    formula_str <- paste(y, "~ L(", y, ", ", p, ") + L(", x_var, ", ", q, ")", sep = "")
    
    model <- tryCatch(dynlm(as.formula(formula_str), data = train_data), error = function(e) NULL)
    
    if (!is.null(model)) {
      pred <- predict(model, newdata = test_data)
      actual <- test_data[[y]]
      mse_values[i] <- (pred - actual)^2
    }
  }
  
  return(mean(mse_values, na.rm = TRUE))
}

# optimal lag for each independent variable
for (var in colnames(x)) {
  q_results <- data.frame(q = q_range, mse = NA) 
  
  for (q in q_range) {
    p <- best_p  # using the best lag for the dependent variable
    mse_value <- loocv_mse_single(p, q, y, var, df)
    q_results$mse[q_results$q == q] <- mse_value
  }
  
  best_q_index <- which.min(q_results$mse)
  best_q <- q_results$q[best_q_index]
  best_mse <- q_results$mse[best_q_index]
  
  best_results <- rbind(best_results, data.frame(variable = var, best_q = best_q, best_mse = best_mse))
}

print(best_results)

# formula for the final ADL model using optimal lags
formula_str <- paste(y, "~ L(", y, ",", best_p, ") + ",   
                     paste(mapply(function(var, q) paste0("L(", var, ",", q, ")"),   
                                  best_results$variable, best_results$best_q), collapse = " + "))

final_model <- dynlm(as.formula(formula_str), data = df)
summary(final_model)
