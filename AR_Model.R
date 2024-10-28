# Install and load necessary libraries
required_packages <- c("forecast", "dplyr", "dynlm", "ggplot2", "plotly")

# Check for missing packages and install them
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

# Load the libraries
library(forecast)
library(dplyr)
library(dynlm)
library(ggplot2)
library(plotly)

# Set parameters
csv_file_path <- "Final_Transformed_Dataset.csv"
target_var <- "INDPRO"
# x_vars <- c("IPFPNSS", "IPMAT", "IPDMAT", "IPMANSICS", "UNRATE",
#             "CES1021000001", "AWOTMAN", "CES0600000008",
#             "CES3000000008", "CIVPART", "LNS11300036", "FEDFUNDS", "COVID_Dummy")
x_vars <- "COVID_Dummy"
max_lags <- 12
forecast_steps <- 12

# Load the CSV file
df <- read.csv(csv_file_path, stringsAsFactors = FALSE)

# Convert 'sasdate' column to Date format
df$sasdate <- as.Date(df$sasdate, format = "%m/%d/%Y")

# Fill null values in numeric columns with the mean of those columns
for (col in names(df)) {
  if (is.numeric(df[[col]])) {
    mean_value <- mean(df[[col]], na.rm = TRUE)
    df[[col]][is.na(df[[col]])] <- mean_value
  }
}

# Check for null values after filling
null_values_after <- sapply(df, function(x) sum(is.na(x)))

# Create a dummy variable for the COVID period
df$COVID_Dummy <- ifelse(df$sasdate >= as.Date("2020-03-01") & df$sasdate <= as.Date("2020-08-01"), 1, 0)

# Function to perform LOOCV and calculate MSE for an AR model
calculate_loocv_mse_ar <- function(df, target_var, max_lags) {
  mse_results <- numeric(max_lags)
  
  for (p in 1:max_lags) {
    mse_fold <- numeric(nrow(df))
    
    # Create lagged variables
    for (lag in 1:p) {
      df[paste0(target_var, "_lag", lag)] <- dplyr::lag(df[[target_var]], lag)
    }
    
    for (i in 1:nrow(df)) {
      # Exclude the i-th observation for LOOCV
      train_data <- df[-i, ]
      test_data <- df[i, , drop = FALSE]
      
      # Formula for AR model with lags of the target variable only
      model_formula <- as.formula(paste(target_var, "~", 
                                        paste(paste0(target_var, "_lag", 1:p), collapse = " + ")))
      
      # Fit the model
      model <- dynlm(model_formula, data = train_data)
      
      # Predict the target variable for the left-out observation
      prediction <- predict(model, newdata = test_data)
      
      # Calculate squared error for this fold
      mse_fold[i] <- (test_data[[target_var]] - prediction)^2
    }
    
    # Store the mean squared error for the current lag
    mse_results[p] <- mean(mse_fold, na.rm = TRUE)
    
    # Remove lagged variables to avoid contamination in the next loop
    for (lag in 1:p) {
      df[[paste0(target_var, "_lag", lag)]] <- NULL
    }
  }
  
  return(mse_results)
}

# Calculate MSE for different lag orders
mse_results <- calculate_loocv_mse_ar(df, target_var, max_lags)

# Find the optimal number of lags (the one with the lowest MSE)
optimal_lag <- which.min(mse_results)
optimal_mse <- min(mse_results)

# Output results
cat("Optimal number of lags:", optimal_lag, "\n")
cat("Minimum MSE:", optimal_mse, "\n")

# Plot MSE results for visualization
plot(1:max_lags, mse_results, type = "b", xlab = "Number of Lags", ylab = "MSE",
     main = "MSE for Different Lag Orders")

# Create lagged variables for AR model
for (i in 1:max_lags) {
  df <- df %>%
    mutate(!!paste0("INDPRO_lag", i) := lag(INDPRO, i))
}

# Drop any rows with NA values introduced by lags
df <- na.omit(df)

# Build the AR model
ar_formula <- as.formula(paste(target_var, "~", 
                               paste(paste0("INDPRO_lag", 1:optimal_lag), collapse = " + "), 
                               "+", paste(x_vars, collapse = " + ")))
ar_model <- dynlm(ar_formula, data = df)

# AR Model Summary
cat("AR Model Summary:\n")
print(summary(ar_model))

# Generate optimal number of lagged variables for AR model
for (i in 1:optimal_lag) {
  df[[paste0("INDPRO_lag", i)]] <- lag(df$INDPRO, i)
}

# Drop rows with NA values introduced by the optimal lags
df <- na.omit(df)

# Fit the AR(optimal_lag) model
ar_formula <- as.formula(
  paste("INDPRO ~ COVID_Dummy +", paste(paste0("INDPRO_lag", 1:optimal_lag), collapse = " + "))
)
ar_model <- dynlm(ar_formula, data = df)

# Prepare for 12-step forecast
last_observed <- tail(df, 1)
forecast_df <- data.frame(matrix(ncol = 1, nrow = forecast_steps))
colnames(forecast_df) <- target_var

# Forecast loop for AR(optimal_lag) model
for (i in 1:forecast_steps) {
  new_data <- data.frame(
    COVID_Dummy = ifelse(as.Date(last_observed$sasdate) + i > as.Date("2020-03-01") &
                           as.Date(last_observed$sasdate) + i <= as.Date("2020-08-01"), 1, 0)
  )
  
  # Add the 12 lags for each forecast step
  for (lag in 1:12) {
    lag_value <- if (i <= lag) {
      last_observed[[paste0("INDPRO_lag", lag - i + 1)]]
    } else {
      forecast_df[i - lag, target_var]
    }
    new_data[[paste0("INDPRO_lag", lag)]] <- lag_value
  }
  
  # Predict using the model
  forecast_df[i, target_var] <- predict(ar_model, newdata = new_data)
}

# Fan chart plot function using Plotly
plot_fan_chart_plotly <- function(actual_data, forecast_data, title) {
  # Prepare data for Plotly
  forecast_data$Time <- seq(nrow(actual_data) + 1, nrow(actual_data) + forecast_steps)
  combined_df <- rbind(
    data.frame(Time = seq_len(nrow(actual_data)), Forecast = actual_data[[target_var]], Type = "Actual"),
    data.frame(Time = forecast_data$Time, Forecast = forecast_data[[target_var]], Type = "Forecast")
  )
  
  # Calculate ymin and ymax for prediction intervals
  forecast_data$lower_80 <- forecast_data[[target_var]] * 0.9
  forecast_data$upper_80 <- forecast_data[[target_var]] * 1.1
  forecast_data$lower_90 <- forecast_data[[target_var]] * 0.8
  forecast_data$upper_90 <- forecast_data[[target_var]] * 1.2
  
  # Create fan chart with Plotly
  p <- plot_ly() %>%
    add_lines(data = combined_df[combined_df$Type == "Forecast", ], x = ~Time, y = ~Forecast, name = "Forecast", line = list(color = "red")) %>%
    add_ribbons(data = forecast_data, 
                x = ~Time, 
                ymin = ~lower_80, ymax = ~upper_80,
                fillcolor = "rgba(0, 0, 255, 0.3)",
                line = list(color = "transparent"),
                name = "80% Prediction Interval") %>%
    add_ribbons(data = forecast_data, 
                x = ~Time, 
                ymin = ~lower_90, ymax = ~upper_90,
                fillcolor = "rgba(173, 216, 230, 0.2)",
                line = list(color = "transparent"),
                name = "90% Prediction Interval") %>%
    layout(title = title, xaxis = list(title = "Time"), yaxis = list(title = target_var))
  
  return(p)
}

# Plot the fan chart
fan_chart_plot <- plot_fan_chart_plotly(df, forecast_df, "12-Step Ahead Forecast with Fan Chart (AR(1) Model)")
fan_chart_plot
