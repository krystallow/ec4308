library(readxl) 
library(dplyr)
library(ggplot2)
library(corrplot)
library(tseries)
library(writexl)
library(tidyr)

# Loading datasets
Pre <- read_excel("Pre_Dataset.xlsx")[-1,]
Pre <- rename(Pre, Date = sasdate,GDP = INDPRO,Consumption = DPCERA3M086SBEA,Income = RPI,Production = IPFINAL,Manufacturing_Trade = CMRMTSPLx,Employment_Levels = PAYEMS,Unemployment_Rate = UNRATE,Capacity_Utilisation = CUMFNS,BAA = BAA,AAA = AAA,Retail = RETAILx,Employment_Industrial = MANEMP)
Pre_Dataset <- read_excel("Pre_Dataset.xlsx")[-1,-1]
Pre_Dataset <- Pre_Dataset %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
Pre_Dataset <- Pre_Dataset %>%
  mutate(across(everything(), as.numeric))
Pre_Dataset <- rename(Pre_Dataset,GDP = INDPRO,Consumption = DPCERA3M086SBEA,Income = RPI,Production = IPFINAL,Manufacturing_Trade = CMRMTSPLx,Employment_Levels = PAYEMS,Unemployment_Rate = UNRATE,Capacity_Utilisation = CUMFNS,BAA = BAA,AAA = AAA,Retail = RETAILx,Employment_Industrial = MANEMP)
Transformed_Dataset <- read_excel("Final_Transformed_Dataset.xlsx")[-1,-1]
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), as.numeric))
Transformed_Dataset <- rename(Transformed_Dataset,GDP = INDPRO,Consumption = DPCERA3M086SBEA,Income = RPI,Production = IPFINAL,Manufacturing_Trade = CMRMTSPLx,Employment_Levels = PAYEMS,Unemployment_Rate = UNRATE,Capacity_Utilisation = CUMFNS,BAA = BAA,AAA = AAA,Retail = RETAILx,Employment_Industrial = MANEMP)

####################################################
####### 1a. Auto Correlation Plot for INDPRO #######
####################################################
#Pre:
# Calculate autocorrelation for the 'INDPRO' column without plotting 
acf_results_pre <- lapply(Pre_Dataset['GDP'], acf, plot = FALSE)
acf_results_pre 
# Plotting ACF
par(mfrow = c(1, 1)) 
lapply(names(Pre_Dataset['GDP']), function(col) {
  acf(Pre_Dataset[[col]], main = paste("ACF of", col), lag.max = 12, col = "blue", lwd = 2)
})

#Transformed:
acf_results_transformed <- lapply(Transformed_Dataset['GDP'], acf, plot = FALSE)
acf_results_transformed 
par(mfrow = c(1, 1)) 
lapply(names(Transformed_Dataset['GDP']), function(col) {
  acf(Transformed_Dataset[[col]], main = paste("ACF of", col), lag.max = 12, col = "blue", lwd = 2)
})

####################################################
###### 1b. Augmented Dickey-Fuller (ADF) test ######
####################################################
# Pre:
adf_result <- adf.test(Pre_Dataset$GDP)
adf_result
# The Dickey-Fuller statistic is positive in this case, which is a key indicator that the time series might be non-stationary. 
# p-value = 0.99 -> not statistically significant & cannot reject the numm hypothesis that the series has a unit root
# The lag order refers to the number of lagged differences included in the test. 
# The ADF test is "augmented" with lagged values of the differenced series to account for higher-order autocorrelation that might be present in the data. 
# lag order 9 -> 9 lagged differences of the INDPRO series were included in the test to account for autocorrelation.

# Transformed: 
adf_result <- adf.test(Transformed_Dataset$GDP)
adf_result
# A more negative Dickey-Fuller statistic indicates stronger evidence against the presence of a unit root (non-stationarity)
# -7.3901 is quite negative, suggesting strong evidence that the series is stationary and does not have a unit root.
# p-value = 0.01 (significant) -> only a 1% chance that the INDPRO time series is non-stationary. 
# reject null hypothesis, conclude that it is likely stationary

####################################################
##### 2. Correlation matrix between variables #####
####################################################
subset_data <- Pre_Dataset[, c("GDP", "Consumption", "Income", "Production", "Manufacturing_Trade", "Employment_Levels", "Unemployment_Rate", "Capacity_Utilisation", "BAA", "AAA", "Retail", "Employment_Industrial")]  
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
subset_data <- Pre_Dataset[, c("GDP", "Consumption", "Income", "Production", "Manufacturing_Trade", "Employment_Levels", "Unemployment_Rate", "Capacity_Utilisation", "BAA", "AAA", "Retail", "Employment_Industrial")]  

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
Pre$Date <- as.Date(as.numeric(Pre$Date), origin = "1899-12-30")
Pre$Date <- as.Date(Pre$Date, format = "%Y-%d-%m")
Pre$Date <- gsub("^(\\d{4})-(\\d{2})-(\\d{2})$", "\\1-\\3-\\2", Pre$Date)
Pre$Date <- as.Date(Pre$Date, format = "%Y-%m-%d")
subset_data <- Pre[, c("Date","GDP", "Consumption", "Income", "Production", "Manufacturing_Trade", "Employment_Levels", "Unemployment_Rate", "Capacity_Utilisation", "BAA", "AAA", "Retail", "Employment_Industrial")]  

# Plot INDPRO  
ggplot(subset_data, aes(x = Date, y = GDP)) +
  geom_line(color = "blue") +
  labs(title = "GDP over Time", x = "Year", y = "GDP") +
  theme_minimal()

# Get column names of the indicators you want to plot excluding date
indicators <- names(subset_data)[!names(subset_data) %in% "Date"]
# Plotting the other indicators

library(ggplot2)
library(gridExtra)

# Create a list to store the plots
plots <- list()

# Loop through each indicator and generate the plot
for (indicator in indicators) {
  p <- ggplot(subset_data, aes(x = Date, y = .data[[indicator]])) +
    geom_line() +
    labs(title = paste(indicator, "over Time"), x = "Year", y = indicator) +
    theme_minimal()
  
  plots[[indicator]] <- p
}

# Arrange the plots in a grid with 3 rows and 4 columns
grid.arrange(grobs = plots, nrow = 3, ncol = 4)


