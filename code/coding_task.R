library(dplyr)
library(vegan)
library(ggplot2)

# Set working directory
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

##### Read and set up metadata DataFrame
metadata <- getDataFromTable("Metadata.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

metadata <-transform(metadata, ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth))

metadata <- mutate(metadata, sample_name = paste(ATTRIBUTE_Location, ATTRIBUTE_Depth, sep = "_"))

metadata_surface <- metadata %>%
  transform(ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth)) %>%
  filter(ATTRIBUTE_Depth_Range == "0-30") %>%
  arrange(ATTRIBUTE_Location, ATTRIBUTE_Depth)

##### Read ASVs data
# Load 16s AVSs data
asv_16s <- getDataFromTable("ASV_16S.csv", vars_toDrop_beg = 0, vars_toDrop_end = 8)

# Convert char variables into numeric without loosing DF structure
for (colX in colnames(asv_16s)) {
  asv_16s[colX] <- as.numeric(asv_16s[colX][,1])
}

# Extract surface samples in cycle-day order
surface_asv16s <- asv_16s[row.names(metadata_surface),]

# Calculate alpha diversity indices for surface samples
alpha_diversity_indices <- vegan::estimateR(surface_asv16s)

# I will only use chao1 for this exercise
surface_asv16s_chao1 <- as.numeric(data.frame(t(alpha_diversity_indices))$S.chao1)

# Create DF with necessary data to plot
surface_alpha_diversity <- bind_cols(select(metadata_surface, ATTRIBUTE_Location, sample_name), surface_asv16s_chao1)
colnames(surface_alpha_diversity) <- c("ATTRIBUTE_Location", "sample_name", "Chao1")

# Plot all samples separately
ggplot(surface_alpha_diversity, aes(x=sample_name, y=Chao1))+
  geom_point() +
  labs(x = "Samples",
       title = "Alpha diversity for 16S AVSs of Surface Samples") +
  theme(title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot average alpha diversity per cycle-day
ggplot(surface_alpha_diversity, aes(x=ATTRIBUTE_Location, y=Chao1))+
  geom_boxplot() +
  labs(x = "Samples",
       title = "Average Alpha diversity per cycle-day for 16S AVSs of Surface Samples") +
  theme(title = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




asv_18s <- getDataFromTable("ASV_18SV9.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

metabolites <- getDataFromTable("Metabolites_Feature Table.csv", vars_toDrop_beg = 2, vars_toDrop_end = 0)
