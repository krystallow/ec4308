load("HybridModels.RData")
library(ggplot2)

selected_rlasso_variables <- names(which(coef(rlasso.fit) != 0))

importance_rf <- importance(rffit.ht)
# Convert importance values to a data frame
importance_df <- as.data.frame(importance_rf)
# If IncNodePurity is in the first column, you can name it explicitly
colnames(importance_df) <- c("IncNodePurity")
threshold <- quantile(importance_df$IncNodePurity, 0.75) #75-percentile threshold 
# Convert importance scores to a data frame and filter based on a threshold
importance_df <- as.data.frame(importance_rf)
significant_rf_variables <- rownames(importance_df[importance_df$MeanDecreaseGini > threshold, ])

# Combine variables from rLASSO and significant RF variables
hybrid_selected_variables <- unique(c(selected_rlasso_variables, significant_rf_variables))
# Filter importance values for selected variables
hybrid_importance_df <- importance_df[rownames(importance_df) %in% hybrid_selected_variables, , drop = FALSE]

# Add variable names as a column for easier plotting
custom_labels <- c(
  "CES2000000008" = "Construction Wage",
  "CES3000000008" = "Manufacturing Wage",
  "IPMAT" = "Materials Production",
  "UNRATE" = "Unemployment Rate",
  "CES1021000001" = "Employees: Mining and Logging",
  "IPNMAT" = "Nondurable Goods Production",
  "UEMP15T26" = "Unemployed for 15-26 Weeks",
  "USGOOD" = "Employees: Goods Producing",
  "DPCERA3M086SBEA" = "Real Personal Consumption Expenditure",
  "CE16OV" = "Employment Level",
  "IPMANSICS" = "Manufacturing Production",
  "AWHMAN" = "Average Weekly Hours of Production",
  "IPDCONGD" = "Durable Goods Production",
  "IPBUSEQ" = "Business Equipment Production",
  "IPFPNSS" = "Final and Nonindustrial production",
  "PERMITMW" = "New Prvately-Owned Housing (Midwest Region)",
  "IPFUELS" = "Fuel Production",
  "GS5" = "Yield on Treasury Securities (5-Year Constant Maturity)",
  "EXSZUSx" = "Switzerland/US Exchange Rate",
  "RETAILx" = "Retail and Food Service Sales",
  "GS10" = "Yield on Treasury Securities (10-Year Constant Maturity)",
  "(Intercept)" = "Intercept"
)
# Replace rownames (variable names) in the filtered data frame with custom labels
hybrid_importance_df$Variable <- rownames(hybrid_importance_df)
hybrid_importance_df$Variable <- custom_labels[hybrid_importance_df$Variable]

# Plotting the relative importance
ggplot(hybrid_importance_df, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Relative Importance of Selected Variables (Hybrid Model)",
    x = "Variable",
    y = "Node Purity"
  ) +
  theme_minimal()