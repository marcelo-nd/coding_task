library(readr)
library(dplyr)

metadata <- read_csv("C:/Users/marce/Desktop/Metadata.csv", col_names = FALSE)

metadata <- t(metadata)

colnames(metadata) <- metadata[1,]

metadata <- metadata[2:(nrow(metadata)-10),]

rownames(metadata) <- metadata[,1]

metadata <- data.frame(metadata[,-1])


getDataFromTable <- function(data_path, vars_toDrop_beg = 0, vars_toDrop_end = 0){
  dataTable <- readr::read_csv(data_path, col_names = FALSE)
  dataTable <- t(dataTable)
  colnames(dataTable) <- dataTable[1,]
  dataTable <- dataTable[(2+vars_toDrop_beg):(nrow(dataTable)-vars_toDrop_end),]
  rownames(dataTable) <- dataTable[,1]
  return(data.frame(dataTable[,-1]))
}

meta <- getDataFromTable("C:/Users/marce/Desktop/Metadata.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

metadata_surface <- metadata %>% 
  filter(ATTRIBUTE_Depth_Range == "0-30") %>%
  arrange(ATTRIBUTE_Location)
