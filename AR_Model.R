library(forecast)
library(tseries)
library(plotly)
library(readxl)
Final_Transformed_Dataset <- read_excel("Final_Transformed_Dataset.xlsx")[-1,c(1,2)]
Final_Transformed_Dataset$sasdate <- as.Date(as.numeric(Final_Transformed_Dataset$sasdate), origin = "1899-12-30")
Final_Transformed_Dataset$sasdate <- as.Date(Final_Transformed_Dataset$sasdate, format = "%Y-%d-%m")
Final_Transformed_Dataset$sasdate <- gsub("^(\\d{4})-(\\d{2})-(\\d{2})$", "\\1-\\3-\\2", Final_Transformed_Dataset$sasdate)
Final_Transformed_Dataset$sasdate <- as.Date(Final_Transformed_Dataset$sasdate, format = "%Y-%m-%d")
df <- Final_Transformed_Dataset

# find optimal number of lags
ts_data <- ts(df$INDPRO, start = c(1960, 3), frequency = 12)

# creating a dummy variable for the COVID period
df$COVID_Dummy <- ifelse(df$Date >= as.Date("2020-01-01") & df$Date <= as.Date("2020-07-01"), 1, 0)
# Convert the dummy variable to a time series
covid_dummy_ts <- ts(df$COVID_Dummy, start = c(1960, 3), frequency = 12)

# function to calculate LOOCV for AR models including exogenous regressors (COVID dummy)
loocv_mse_with_dummy <- function(ts_data, covid_dummy_ts, max_lag) {
  mse_per_lag <- numeric(max_lag)
  
  for (lag in 1:max_lag) {
    errors <- numeric(length(ts_data) - lag)
    
    for (i in (lag + 1):length(ts_data)) {
      train_set <- ts_data[1:(i - 1)]
      train_dummy <- covid_dummy_ts[1:(i - 1)]  # Corresponding dummy variables
      
      # Fit AR model with the specified lag and COVID dummy as an exogenous regressor
      fit <- try(Arima(train_set, order = c(lag, 0, 0), xreg = train_dummy), silent = TRUE)
      
      if (!inherits(fit, "try-error")) {
        forecasted_value <- forecast(fit, h = 1, xreg = covid_dummy_ts[i])$mean
        actual_value <- ts_data[i]
        errors[i - lag] <- (forecasted_value - actual_value)^2
      } else {
        errors[i - lag] <- NA
      }
    }
    
    mse_per_lag[lag] <- mean(errors, na.rm = TRUE)
  }
  
  optimal_lag <- which.min(mse_per_lag)
  return(optimal_lag)
}

# Find the optimal lag including the COVID dummy variable (testing up to lag 8)
optimal_lag_with_dummy <- loocv_mse_with_dummy(ts_data, covid_dummy_ts, max_lag = 8)
print(paste("Optimal Lag with Dummy:", optimal_lag_with_dummy))

optimal_lag <- optimal_lag_with_dummy
fit <- Arima(ts_data, order = c(optimal_lag, 0, 0), xreg = covid_dummy_ts)
future_covid_dummy <- rep(0, 12)  
forecasted_values <- forecast(fit, h = 12, xreg = future_covid_dummy)
print(forecasted_values)
plot(forecasted_values, main = "12-Step Forecast with COVID Dummy", 
     xlab = "Year", ylab = "INDPRO", xlim = c(2019, 2026))