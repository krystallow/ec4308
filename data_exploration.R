library(dplyr)
library(readr)


### Read in accident data ###
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

print(dim(all_accident_data))  


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


## Check if its correct number of FIPS Code >> https://transition.fcc.gov/oet/info/maps/census/fips/fips.txt
all_accident_data$STATE <- as.factor(all_accident_data$STATE)
all_accident_data$COUNTY <- as.factor(all_accident_data$COUNTY)


unique_counties_per_state <- all_accident_data %>%
  group_by(STATE) %>%
  summarise(unique_counties = n_distinct(COUNTY))

print(unique_counties_per_state) ## 51 unique state codes >> correct

# Groupby STATE
accidents_by_state <- all_accident_data %>%
  group_by(STATE) %>%
  summarise(total_accidents = n()) %>%
  arrange(desc(total_accidents))


## Groupby STATE, YEAR, MONTH
all_accident_data$YEAR <- as.factor(all_accident_data$YEAR)
all_accident_data$MONTH <- as.factor(all_accident_data$MONTH)

accidents_by_state_year_month <- all_accident_data %>%
  group_by(STATE, YEAR, MONTH) %>%
  summarise(total_accidents = n()) %>%
  arrange(STATE, YEAR, MONTH)



