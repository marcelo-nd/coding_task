
# This function reads a table from "data_path" and transposes it.
# Assumes first column in original is important such as IDs
# Drops "vars_toDrop_beg" number of columns (counting from the 2nd col) from the original table
# Also, drops "vars_toDrop_end" number of columns at the end of the original table if they are not needed.
getDataFromTable <- function(data_path, vars_toDrop_beg = 0, vars_toDrop_end = 0){
  dataTable <- readr::read_csv(data_path, col_names = FALSE, )
  # transpose table to achieve the format: rows = samples and cols = variables
  dataTable <- t(dataTable)
  # setting first row as colnames
  colnames(dataTable) <- dataTable[1,]
  # Dropping rows that are not used, including the column where colnames where
  dataTable <- dataTable[(2+vars_toDrop_beg):(nrow(dataTable)-vars_toDrop_end),]
  # erasing rownames because, ,if present, cannot be assigned with column_to_rownames
  dataTable <- data.frame(dataTable, row.names = NULL)
  # Return df with correct rownames and removing column where rownames where
  return(tibble::column_to_rownames(dataTable, colnames(dataTable)[1]))
}

# This function converts all columns to numeric by iterating. This prevents structure changes induced by other functions.
convert_all_cols_to_numeric <- function(dataframe){
  # Iterate over all columns
  for (colX in colnames(dataframe)) {
    #print(colX)
    # If colX char variable is not null
    if (!is.na(colX)) {
      # reassing the column after transforming it to numeric
      dataframe[colX] <- as.numeric(dataframe[colX][,1])
    }
  }
  
  return(dataframe)
}


#' Calculates Shannon for a microbiome or metabolite table 
#' 
#' @param dataframe Dataframe containing microbiome or metabolite data in rows = samples cols= variables format.
#' @param metadata_df Dataframe containg metadata including sample names in the same order than "dataframe".
#' @param useful_metadata_cols A list of strings that are the names of the variables to be included in the resulting dataframe.
#' @return A dataframe containgn Shanon index and metadata of samples.
get_alpha_diversity <- function(dataframe, metadata_df, useful_metadata_cols){
  # Calculate shannon index
  shannon <- vegan::diversity(dataframe)
  
  # Check samples are in same order. Here I should return error and stopping
  #row.names(dataframe) == colnames(shannon)
  
  # Create DF with necessary data to plot
  alpha_diversity <- bind_cols(select(metadata_surface, useful_metadata_cols), shannon)
  colnames(alpha_diversity) <- c(useful_metadata_cols, "Shannon")
  
  return(alpha_diversity)
}

get_alpha_diversity_chao <- function(dataframe, metadata_df, useful_metadata_cols){
  alpha_diversity_indices <- vegan::estimateR(dataframe)
  
  # Check samples are in same order
  row.names(dataframe) == colnames(alpha_diversity_indices)
  
  # I will only use chao1 for this exercise
  chao1 <- as.numeric(data.frame(t(alpha_diversity_indices))$S.chao1)
  
  # Create DF with necessary data to plot
  alpha_diversity <- bind_cols(select(metadata_surface, useful_metadata_cols), chao1)
  colnames(alpha_diversity) <- c(useful_metadata_cols, "Chao1")
  
  return(alpha_diversity)
}

# Simple function to determine if x is less than value.
is_less_than_value <- function(x, value){
  return(x<value)
}
