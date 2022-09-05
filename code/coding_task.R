library(dplyr)
library(vegan)
library(ggplot2)

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
surface_samples_in_order <- rownames(metadata_surface)
surface_samples_in_order

# Reading all data
asv_16s <- getDataFromTable("ASV_16S.csv", vars_toDrop_beg = 0, vars_toDrop_end = 8)

#Converting char variables into numeric without loosing DF structure
for (colX in colnames(asv_16s)) {
  print(colX)
  asv_16s[colX] <- as.numeric(asv_16s[colX][,1])
}

# Extracting surface samples in cycle-day order
surface_asv16s <- asv_16s[surface_samples_in_order,]


data_richness <- estimateR(surface_asv16s)

t(data_richness)

data.frame(t(data_richness))$S.chao1

test2 <- data.frame(cbind(row.names(surface_asv16s), as.numeric(data.frame(t(data_richness))$S.chao1)))

test2[,2] <- as.numeric(test2[,2])

ggplot(test2, aes(x=X1, y=X2))+
  geom_point()







asv_18s <- getDataFromTable("ASV_18SV9.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

metabolites <- getDataFromTable("Metabolites_Feature Table.csv", vars_toDrop_beg = 2, vars_toDrop_end = 0)
