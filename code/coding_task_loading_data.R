library(dplyr)
library(vegan)
library(ggplot2)

# Set working directory
setwd("C:/Users/marce/Documents/Repos/coding_task/data")

source("C:/Users/marce/Documents/Repos/coding_task/code/helper_functions.R")

##### Read and set up metadata DataFrame #####
##########################################################################################
metadata <- getDataFromTable("Metadata.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

# Convert ATTRIBUTE_Depth to numeric
metadata <-transform(metadata, ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth))

# Create new variable (sample_name) that includes cycle + day + depth information
metadata <- mutate(metadata, sample_name = paste(ATTRIBUTE_Location, ATTRIBUTE_Depth, sep = "_"))

# Select only surface samples (Depth <= 30)
metadata_surface <- metadata %>%
  transform(ATTRIBUTE_Depth = as.numeric(ATTRIBUTE_Depth)) %>%
  filter(ATTRIBUTE_Depth_Range == "0-30") %>%
  arrange(ATTRIBUTE_Location, ATTRIBUTE_Depth)

##########################################################################################

##### Read ASVs data #####
##########################################################################################
# 16s AVSs =================================

# Load data
asv_16s <- getDataFromTable("ASV_16S.csv", vars_toDrop_beg = 0, vars_toDrop_end = 8)

# Convert char variables into numeric without loosing DF structure
asv_16s <- convert_all_cols_to_numeric(asv_16s)

# Extract surface samples in cycle-day-depth order
surface_asv16s <- asv_16s[row.names(metadata_surface),]

# Check samples are in same order
row.names(surface_asv16s) == row.names(metadata_surface)

# Get dataframe with metadata and chao1 index information to plot
surface_alpha_diversity_asv16s <- get_alpha_diversity(dataframe = surface_asv16s, metadata_df = metadata_surface, useful_metadata_cols = c("ATTRIBUTE_Location", "sample_name", "ATTRIBUTE_Filament_Possition"))

# 18s AVSs ================================= 
# Load data
asv_18s <- getDataFromTable("ASV_18SV9.csv", vars_toDrop_beg = 0, vars_toDrop_end = 10)

# Convert char variables into numeric without loosing DF structure
asv_18s <- convert_all_cols_to_numeric(asv_18s)

# Extract surface samples in cycle-day-depth order
surface_asv18s <- asv_18s[row.names(metadata_surface),]

# Check samples are in same order
row.names(surface_asv18s) == row.names(metadata_surface)

# Get dataframe with metadata and chao1 index information to plot
surface_alpha_diversity_asv18s <- get_alpha_diversity(dataframe = surface_asv18s, metadata_df = metadata_surface, useful_metadata_cols = c("ATTRIBUTE_Location", "sample_name", "ATTRIBUTE_Filament_Possition"))


##########################################################################################

##### Read metabolites data #####
##########################################################################################
# Load data
metabolites <- getDataFromTable("Metabolites_Feature Table.csv", vars_toDrop_beg = 2, vars_toDrop_end = 0)

# Convert char variables into numeric without loosing DF structure
metabolites <- convert_all_cols_to_numeric(metabolites)

metabolites <- metabolites * 100000

metabolites <- mutate_all(metabolites, round)

metabolites <- mutate_if(metabolites, is.numeric, list(~ifelse(is_less_than_value(., 10), 0, .)))

#hist(metabolites$X4969, breaks=500)

# Extract surface samples in cycle-day-depth order
surface_metabolites <- metabolites[row.names(metadata_surface),]

# Check samples are in same order
row.names(surface_metabolites) == row.names(metadata_surface)

# Get dataframe with metadata and chao1 index information to plot
surface_alpha_diversity_metabolites <- get_alpha_diversity(dataframe = surface_metabolites, metadata_df = metadata_surface, useful_metadata_cols = c("ATTRIBUTE_Location", "sample_name", "ATTRIBUTE_Filament_Possition"))

##########################################################################################