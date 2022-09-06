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

convert_all_cols_to_numeric <- function(dataframe){
  for (colX in colnames(dataframe)) {
    dataframe[colX] <- as.numeric(dataframe[colX][,1])
  }
  
  return(dataframe)
}

get_alpha_diversity <- function(dataframe, metadata_df, useful_metadata_cols){
  shannon <- vegan::diversity(dataframe)
  
  # Check samples are in same order
  row.names(dataframe) == colnames(shannon)
  
  # I will only use chao1 for this exercise
  #chao1 <- as.numeric(data.frame(t(alpha_diversity_indices))$S.chao1)
  
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

is_less_than_value <- function(x, value){
  return(x<value)
}
