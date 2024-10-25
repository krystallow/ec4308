################
### ADL Model###
################

# Load required libraries
library(readxl)
library(dynlm)

# Load and clean data
df <- na.omit(read_excel("Final_Transformed_Dataset.xlsx"))
df <- df[-1, ]  
selected_variables <- read.csv("variables_selected.csv")
selected_cols <- c("INDPRO", colnames(selected_variables))
df <- df[, selected_cols]

y <- df$INDPRO
x <- df[, setdiff(colnames(df), "INDPRO")]
x <- as.matrix(x)

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

# Identify best p for y
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

# Construct final formula based on best p for y and best q for each predictor
formula_str <- paste("y ~ L(y,", best_p, ") + ", 
                     paste(mapply(function(var, q) paste0("L(", var, ",", q, ")"), 
                                  best_results$variable, best_results$best_q), collapse = " + "))


final_model <- dynlm(as.formula(formula_str), data = df)
summary(final_model)
