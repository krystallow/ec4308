################
### ADL Model###
################

################################
###### Optimal p = 9, q = 3 ############ >> With dimension warning >> Overfit
################################

library(readxl)
library(dynlm)
library(car)

# Load and clean data
df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]  
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]
df <- df[, !(names(df) %in% c("IPFPNSS", "IPMAT", "IPMANSICS"))]

y <- df$INDPRO
x <- df[, setdiff(colnames(df), "INDPRO")]

## Multicollinearity Check##
correlation_matrix <- cor(x)
print(correlation_matrix)

# Check VIF values for multicollinearity
vif_values <- vif(lm(y ~ ., data = df))
print(vif_values)


# Define cross-validation function with rolling window
cv_rmse <- function(p, q, y, x_var, data) {
  n <- nrow(data)
  rmse <- numeric()  #
  
  # Rolling window cross-validation
  for (train_end in (p + max(q) + 1):(n - 1)) {
    # Define training and test sets
    train_indices <- 1:train_end
    test_index <- train_end + 1
    
    # Construct formula dynamically with lagged variables
    formula_str <- paste0("y ~ L(y, ", p, ") + L(x_var, ", q, ")")
    
    model <- dynlm(as.formula(formula_str), data = data[train_indices, ])
    
    # Predict on test set
    test_data <- data[test_index, , drop = FALSE]
    pred <- as.numeric(predict(model, newdata = test_data))
    
    actual <- data$INDPRO[test_index]
    rmse <- c(rmse, sqrt((pred - actual)^2))
  }
  
  return(mean(rmse, na.rm = TRUE))
}


# Initialise to store best lag (q) for each independent variable
best_results <- data.frame(variable = character(), best_q = numeric(), best_rmse = numeric(), stringsAsFactors = FALSE)

p_range <- 1:12  # Dependent variable lag
q_range <- 1:6  # Independent variable lags

# Find the optimal lag p for the dependent variable (y)
p_results <- data.frame(p = p_range, rmse = NA)  # Store RMSE for each p

for (p in p_range) {
  # Cross-validate with only lagged y
  rmse_value <- cv_rmse(p, 0, df$INDPRO, df$INDPRO, df)
  p_results$rmse[p_results$p == p] <- rmse_value
}

best_p_index <- which.min(p_results$rmse)
best_p <- p_results$p[best_p_index]
cat("Best lag p for y:", best_p, "\n")

# Loop through each independent variable to determine optimal q lag
for (var in colnames(x)) {
  q_results <- data.frame(q = q_range, rmse = NA)  # Store RMSE for each q of the variable
  
  for (q in q_range) {
    # Cross-validate with best_p for y and optimizing q for current variable
    rmse_value <- cv_rmse(best_p, q, df$INDPRO, df[[var]], df)
    q_results$rmse[q_results$q == q] <- rmse_value
  }
  
  # Find best q for the current variable
  best_q_index <- which.min(q_results$rmse)
  best_q <- q_results$q[best_q_index]
  best_rmse <- q_results$rmse[best_q_index]
  
  # Append to results
  best_results <- rbind(best_results, data.frame(variable = var, best_q = best_q, best_rmse = best_rmse))
}

print(best_results)

formula_str <- paste("y ~ L(y,", best_p, ") + ", 
                     paste(mapply(function(var, q) paste0("L(", var, ",", q, ")"), 
                                  best_results$variable, best_results$best_q), collapse = " + "))


final_model <- dynlm(as.formula(formula_str), data = df)
summary(final_model)




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



############################################
###### Optimal Lags = 1 for all ############ >>> AIC
############################################
library(readxl)
library(dynlm)
library(car)

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]   
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]

remove_multicollinear <- function(data, vif_threshold = 5) {
  repeat {
    vif_values <- vif(lm(INDPRO ~ ., data = data))
    max_vif <- max(vif_values)
    if (max_vif < vif_threshold) break
    var_to_remove <- names(vif_values)[which.max(vif_values)]
    data <- data[, !colnames(data) %in% var_to_remove]
  }
  return(data)
}

# Apply multicollinearity check
df <- remove_multicollinear(df)
print("VIF check complete; remaining variables:")
print(names(df))

# Define the target (y) and predictor (x) variables
y <- "INDPRO"
x_vars <- setdiff(colnames(df), y)

# Function to find optimal lag p for the dependent variable using AIC
aic_single <- function(p, y, data) {
  formula_str <- paste(y, "~ L(", y, ", ", p, ")", sep = "")
  model <- try(dynlm(as.formula(formula_str), data = data), silent = TRUE)
  if (inherits(model, "try-error")) return(Inf)  # Return a large value if the model fails
  aic_value <- AIC(model)
  cat("p =", p, "AIC =", aic_value, "\n")  # Debugging output
  return(aic_value)
}

p_range <- 1:6
p_aic_values <- sapply(p_range, function(p) aic_single(p, y, df))
best_p <- p_range[which.min(p_aic_values)]
cat("Best lag p for y:", best_p, "\n")

# Function to find optimal q lag for each predictor while holding best_p constant using AIC
aic_single_q <- function(p, q, y, x_var, data) {
  formula_str <- paste(y, "~ L(", y, ", ", p, ") + L(", x_var, ", ", q, ")", sep = "")
  model <- try(dynlm(as.formula(formula_str), data = data), silent = TRUE)
  if (inherits(model, "try-error")) return(Inf)  # Return a large value if the model fails
  aic_value <- AIC(model)
  cat("p =", p, "q =", q, "AIC =", aic_value, "\n")  # Debugging output
  return(aic_value)
}

# Loop to find best q lag for each independent variable
best_results <- data.frame(variable = character(), best_q = numeric(), best_aic = numeric(), stringsAsFactors = FALSE)

for (var in x_vars) {
  q_aic_values <- sapply(1:4, function(q) aic_single_q(best_p, q, y, var, df))
  best_q <- 1:4[which.min(q_aic_values)]
  best_aic <- min(q_aic_values, na.rm = TRUE)
  best_results <- rbind(best_results, data.frame(variable = var, best_q = best_q, best_aic = best_aic))
}

print("Best lag q for each variable:")
print(best_results)

formula_str <- paste(y, "~ L(", y, ",", best_p, ")", sep = "")
for (i in seq_along(best_results$variable)) {
  formula_str <- paste(formula_str, "+ L(", best_results$variable[i], ",", best_results$best_q[i], ")", sep = "")
}

final_model <- try(dynlm(as.formula(formula_str), data = df), silent = TRUE)
summary(final_model)



#################
################# >> Using autoARIMA()
#################
library(readxl)
library(forecast)

df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]  
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]


# optimal lag for INDPRO
best_model <- auto.arima(df$INDPRO)
summary(best_model)

optimal_lags <- list()
for (var in names(df)[names(df) != "INDPRO"]) {
  model <- auto.arima(df[[var]], xreg = df$INDPRO)
  p <- model$arma[1]  # AR order
  q <- model$arma[2]  # MA order
  optimal_lags[[var]] <- list(p = p, q = q)
}

adl_formula <- "INDPRO ~"
#p_indpro <- best_model$arma[1]
#for (lag in 1:p_indpro) {
  #adl_formula <- paste(adl_formula, paste0("lag(INDPRO, ", lag, ")"), sep = " + ")
#}

# Loop through each independent variable to construct the formula
for (var in names(optimal_lags)) {
  p <- optimal_lags[[var]]$p
  q <- optimal_lags[[var]]$q
  
  # Add terms for the independent variable with the specified lags
  for (lag in 1:p) {
    adl_formula <- paste(adl_formula, paste0("lag(", var, ", ", lag, ")"), sep = " + ")
  }
}

adl_formula <- as.formula(adl_formula)

adl_model <- lm(adl_formula, data = df)
summary(adl_model)



## Model validity 
set.seed(2457829)

train_index <- sample(1:nrow(df), 0.7 * nrow(df))
train_data <- df[train_index, ]
test_data <- df[-train_index, ]

# Fit the model on training data
adl_model_train <- lm(adl_formula, data = train_data)

# Predict on the test data
predictions <- predict(adl_model_train, newdata = test_data)

# Calculate RMSE (Root Mean Squared Error) for evaluation
rmse <- sqrt(mean((predictions - test_data$INDPRO)^2))
cat("RMSE on test data:", rmse)


# Residual diagnostics
par(mfrow = c(2, 2))  # Set up the plotting area for multiple plots

# Plot residuals vs fitted values
plot(adl_model_train$fitted.values, resid(adl_model_train), 
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# Q-Q plot for normality of residuals
qqnorm(resid(adl_model_train))
qqline(resid(adl_model_train), col = "red")

# Scale-Location plot
plot(adl_model_train$fitted.values, sqrt(abs(resid(adl_model_train))), 
     xlab = "Fitted Values", ylab = "Sqrt of |Residuals|", 
     main = "Scale-Location")
abline(h = 0, col = "red")

# Residuals vs Leverage plot
plot(hatvalues(adl_model_train), resid(adl_model_train), 
     xlab = "Leverage", ylab = "Residuals", 
     main = "Residuals vs Leverage")
abline(h = 0, col = "red")
