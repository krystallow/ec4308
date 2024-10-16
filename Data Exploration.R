library(readxl) 
library(dplyr)
library(ggplot2)
library(corrplot)

Transformed_Dataset <- read_excel("Transformed_Dataset.xlsx")[-1,-1]
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
Transformed_Dataset <- Transformed_Dataset %>%
  mutate(across(everything(), as.numeric))

# Auto Correlation
acf_results <- lapply(Transformed_Dataset, acf, plot = FALSE)  # Calculate autocorrelation for each variable
acf_results

par(mfrow = c(1, 1))  # Adjust this based on how many plots you want to display at once

lapply(names(Transformed_Dataset), function(col) { # Plot the autocorrelation for each variable
  acf(Transformed_Dataset[[col]], main = paste("ACF of", col))
})
  
# Correlation between variables
correlation_matrix <- cor(Transformed_Dataset_numeric, use = "pairwise.complete.obs") # for all variables

subset_data <- Transformed_Dataset[, 1:5]  # Adjust the index as necessary (this is for the first 5 variables)
correlation_matrix_subset <- cor(subset_data, use = "pairwise.complete.obs")
correlation_matrix_subset

# Heatmap of the correlation matrix
corrplot(correlation_matrix_subset, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black")  
