library(dplyr)

setwd("C:/Users/marce/Documents/Repos/coding_task/data")

getDataFromTable <- function(data_path, vars_toDrop_beg = 0, vars_toDrop_end = 0){
  dataTable <- readr::read_csv(data_path, col_names = FALSE, )
  dataTable <- t(dataTable)
  colnames(dataTable) <- dataTable[1,]
  # Dropping rows that are not used, including the column where colnames where
  dataTable <- dataTable[(2+vars_toDrop_beg):(nrow(dataTable)-vars_toDrop_end),]
  dataTable <- data.frame(dataTable, row.names = NULL)
  # Return df with correct rownames and removing column where rownames where
  return(tibble::column_to_rownames(dataTable, colnames(dataTable)[1]))
}

# Reading and setting up metadata df 
metadata <- getDataFromTable("Metadata.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

metadata <-transform(metadata, ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth))

metadata_surface <- metadata %>%
  transform(ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth)) %>%
  filter(ATTRIBUTE_Depth_Range == "0-30") %>%
  arrange(ATTRIBUTE_Location, ATTRIBUTE_Depth)

# Extracting list of samples in order of cycle and day.
samples_in_order <- rownames(metadata_surface)
samples_in_order


