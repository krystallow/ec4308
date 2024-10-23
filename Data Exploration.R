library(readxl) 
library(dplyr)
library(ggplot2)
library(corrplot)
library(tseries)
library(writexl)
library(tidyr)

# Loading datasets
Pre_Dataset <- read_excel("Pre_Dataset.xlsx")[-1,-1]
Pre_Dataset <- Pre_Dataset %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
Pre_Dataset <- Pre_Dataset %>%
  mutate(across(everything(), as.numeric))
Transformed_Dataset <- read_excel("Final_Transformed_Dataset.xlsx")[-1,-1]
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), as.numeric))

####################################################
####### 1a. Auto Correlation Plot for INDPRO #######
####################################################
#Pre:
# Calculate autocorrelation for the 'INDPRO' column without plotting 
acf_results_pre <- lapply(Pre_Dataset['INDPRO'], acf, plot = FALSE)
acf_results_pre 
# Plotting ACF
par(mfrow = c(1, 1)) 
lapply(names(Pre_Dataset['INDPRO']), function(col) {
  acf(Pre_Dataset[[col]], main = paste("ACF of", col), lag.max = 12, col = "blue", lwd = 2)
})

#Transformed:
acf_results_transformed <- lapply(Pre_Dataset['INDPRO'], acf, plot = FALSE)
acf_results_transformed 
par(mfrow = c(1, 1)) 
lapply(names(Transformed_Dataset['INDPRO']), function(col) {
  acf(Transformed_Dataset[[col]], main = paste("ACF of", col), lag.max = 12, col = "blue", lwd = 2)
})

####################################################
###### 1b. Augmented Dickey-Fuller (ADF) test ######
####################################################
# Pre:
adf_result <- adf.test(Pre_Dataset$INDPRO)
adf_result
# The Dickey-Fuller statistic is positive in this case, which is a key indicator that the time series might be non-stationary. 
# p-value = 0.99 -> not statistically significant & cannot reject the numm hypothesis that the series has a unit root
# The lag order refers to the number of lagged differences included in the test. 
# The ADF test is "augmented" with lagged values of the differenced series to account for higher-order autocorrelation that might be present in the data. 
# lag order 9 -> 9 lagged differences of the INDPRO series were included in the test to account for autocorrelation.

# Transformed: 
adf_result <- adf.test(Transformed_Dataset$INDPRO)
adf_result
# A more negative Dickey-Fuller statistic indicates stronger evidence against the presence of a unit root (non-stationarity)
# -7.3901 is quite negative, suggesting strong evidence that the series is stationary and does not have a unit root.
# p-value = 0.01 (significant) -> only a 1% chance that the INDPRO time series is non-stationary. 
# reject null hypothesis, conclude that it is likely stationary

####################################################
##### 2. Correlation matrix between variables #####
####################################################
subset_data <- Pre_Dataset[, c("INDPRO", "DPCERA3M086SBEA", "RPI", "IPFINAL", "CMRMTSPLx", "PAYEMS", "UNRATE", "CUMFNS", "BAA", "AAA", "RETAILx", "MANEMP")]  
correlation_matrix_subset <- cor(subset_data, use = "pairwise.complete.obs")
correlation_matrix_subset
write_xlsx(as.data.frame(correlation_matrix_subset), "correlation_matrix.xlsx")

# Heatmap of the correlation matrix
corrplot(correlation_matrix_subset, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black")  

####################################################
############ 3. Descriptive Statistics ############# 
####################################################
subset_data <- Pre_Dataset[, c("INDPRO", "DPCERA3M086SBEA", "RPI", "IPFINAL", "CMRMTSPLx", "PAYEMS", "UNRATE", "CUMFNS", "BAA", "AAA", "RETAILx", "MANEMP")]  

descriptive_stats <- function(x) {
  stats <- c(
    Count = length(x),
    Mean = mean(x, na.rm = TRUE),
    `Standard error` = sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))),
    Minimum = min(x, na.rm = TRUE),
    `25%` = quantile(x, 0.25, na.rm = TRUE),
    `50%` = median(x, na.rm = TRUE),
    `75%` = quantile(x, 0.75, na.rm = TRUE),
    Maximum = max(x, na.rm = TRUE)
  )
  return(stats)
}

descriptive_table <- sapply(subset_data, descriptive_stats)
descriptive_table <- as.data.frame(descriptive_table)
descriptive_table <- round(descriptive_table, 2)
descriptive_table <- cbind(Statistic = c("Count", "Mean", "Standard error", "Minimum", "25%", "50%", "75%", "Maximum"), descriptive_table)
descriptive_table
write_xlsx(as.data.frame(descriptive_table), "descriptive_table.xlsx")

####################################################
######### 4. Trends of INDPRO/ X-variables ######### 
####################################################
pre_with_columns <- read_excel("Pre_Dataset.xlsx")[-1,]
subset_data <- pre_with_columns[, c("sasdate", "INDPRO", "DPCERA3M086SBEA", "RPI", "IPFINAL", "CMRMTSPLx", "PAYEMS", "UNRATE", "CUMFNS", "BAA", "AAA", "RETAILx", "MANEMP")]  
filtered_data <-subset_data %>%
  pivot_longer(cols = -sasdate, names_to = "Indicator", values_to = "Value")

# Create the stacked plots using ggplot2
ggplot(filtered_data, aes(x = Date, y = Value, color = dates)) +
  geom_line(size = 1) +  # Create a line plot with different colors for each indicator
  facet_wrap(~ dates, scales = "free_y", ncol = 1) +  # Stack the plots vertically
  theme_minimal() +  # Apply a minimal theme for cleaner visuals
  labs(title = "Trends of Indicators Over Time", x = "Date", y = "Value") +
  theme(legend.position = "none")  # Remove the legend as each plot is labeled individually
