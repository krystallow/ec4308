library(readxl)  
library(writexl)
library(dplyr) 

# Read the dataset
Pre_Dataset <- read_excel("Desktop/Pre_Dataset.xlsx")

# Define transformation functions
apply_transform <- function(x, transform_type) {
  if (transform_type == 1) {
    return(x)  # No transformation
  } else if (transform_type == 2) {
    return(c(NA, diff(x, differences = 1)))  # First difference ∆xt, with NA padding
  } else if (transform_type == 3) {
    return(c(NA, NA, diff(x, differences = 2)))  # Second difference ∆²xt, with NA padding
  } else if (transform_type == 4) {
    return(log(x))  # Log transformation
  } else if (transform_type == 5) {
    return(c(NA, diff(log(x), differences = 1)))  # First difference of log ∆log(xt), with NA padding
  } else if (transform_type == 6) {
    return(c(NA, NA, diff(log(x), differences = 2)))  # Second difference of log ∆²log(xt), with NA padding
  } else if (transform_type == 7) {
    return(c(NA, diff(x / lag(x) - 1)))  # Relative change ∆(xt/xt−1 − 1.0), with NA padding
  } 
}

# Function to apply the transformations and handle row alignment
transformation <- function(Pre_data, output_file) {
  
  # Initialize an empty dataframe with the same structure as the original
  Transformed_data <- data.frame(matrix(ncol = ncol(Pre_data), nrow = nrow(Pre_data)))
  colnames(Transformed_data) <- colnames(Pre_data)
  
  # Copy the first column (sasdate) and first row (transformation row) to the transformed dataset
  Transformed_data[,1] <- Pre_data[,1]  # sasdate
  Transformed_data[1,] <- Pre_data[1,]  # Transformation row
  
  # Extract the transformation row
  transform_row <- as.numeric(unlist(Pre_data[1, -1]))
  
  # Exclude the first row and first column (date column) from Pre_data for transformation
  Pre_data <- as.data.frame(Pre_data[-1,-1])
  
  # Apply transformations column by column
  for (i in 1:ncol(Pre_data)) {  # For each data column
    col_data <- as.numeric(Pre_data[, i])  # Column data excluding the transformation row
    transform_type <- transform_row[i]  # Get the transformation type for this column
    
    # Apply the appropriate transformation
    transformed_col <- apply_transform(col_data, transform_type)
    
    # Store the transformed data, with NA alignment to match the original dataset
    Transformed_data[-1, i + 1] <- transformed_col  # Assigning the transformed column
  }

  # Write the transformed data to an Excel file
  write_xlsx(Transformed_data, output_file)
  #Transformed_data[is.na(Transformed_data)] <- 0
  return(Transformed_data)
}

# Usage
Pre_data <- read_excel("Desktop/Pre_Dataset.xlsx")
transformation(Pre_data, "Desktop/Transformed_Dataset.xlsx")
