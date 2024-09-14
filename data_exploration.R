library(dplyr)
library(readr)

data_path <- "../../Data/FARS Data/"

year_folders <- list.dirs(data_path, full.names = TRUE, recursive = FALSE)

accident_data_list <- list()

for (folder in year_folders) {
  accident_file <- file.path(folder, "ACCIDENT.CSV")
  
  if (file.exists(accident_file)) {
    accident_data <- read_csv(accident_file, col_types = cols(.default = "c")) 
    accident_data_list[[length(accident_data_list) + 1]] <- accident_data
  } else {
    warning(paste("File not found:", accident_file))
  }
}


normalize_data <- function(df) { ## LATTIDENAME saved as char in some, double in some ##
  if ("LATITUDENAME" %in% names(df)) {
    df$LATITUDENAME <- as.character(df$LATITUDENAME)
  }
  return(df)
}

accident_data_list <- lapply(accident_data_list, normalize_data)

all_columns <- unique(unlist(lapply(accident_data_list, names)))

add_missing_columns <- function(df, all_columns) {
  missing_columns <- setdiff(all_columns, names(df))
  df[missing_columns] <- NA
  return(df)
}

accident_data_list <- lapply(accident_data_list, add_missing_columns, all_columns = all_columns)
all_accident_data <- bind_rows(accident_data_list)

print(dim(all_accident_data))  # Print the dimensions of the final dataset



## Check that no rows were dropped ##
rows_per_file <- list()

for (folder in year_folders) {
  accident_file <- file.path(folder, "ACCIDENT.CSV")
  
  if (file.exists(accident_file)) {
    accident_data <- read_csv(accident_file, col_types = cols(.default = "c")) 
    rows_per_file[[basename(folder)]] <- nrow(accident_data)  
  } else {
    warning(paste("File not found:", accident_file))
    rows_per_file[[basename(folder)]] <- 0  
  }
}

total_rows_from_files <- sum(unlist(rows_per_file))
total_rows_combined <- nrow(all_accident_data)

cat("Total rows from individual files:", total_rows_from_files, "\n")
cat("Total rows in combined data frame:", total_rows_combined, "\n")

