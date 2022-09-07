library(dplyr)
library(vegan)
library(ggplot2)

source("C:/Users/marce/Documents/Repos/coding_task/code/helper_functions.R")

##### Plot alpha diversity #####
##########################################################################################
# 16s AVSs =================================
# Plot all samples separately
ggplot(surface_alpha_diversity_asv16s, aes(x=sample_name, y=Shannon))+
  geom_point() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Alpha diversity (Shannon) for 16S AVSs of Surface Samples") +
  theme(title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot average alpha diversity per cycle-day
ggplot(surface_alpha_diversity_asv16s, aes(x=ATTRIBUTE_Location, y=Shannon))+
  geom_boxplot() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Average Alpha diversity (Shannon) per cycle-day for 16S AVSs of Surface Samples") +
  theme(title = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# 18s AVSs =================================
# Plot all samples separately
ggplot(surface_alpha_diversity_asv18s, aes(x=sample_name, y=Shannon))+
  geom_point() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Alpha diversity (Shannon) for 18S AVSs of Surface Samples") +
  theme(title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot average alpha diversity per cycle-day
ggplot(surface_alpha_diversity_asv18s, aes(x=ATTRIBUTE_Location, y=Shannon))+
  geom_boxplot() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Average Alpha diversity (Shannon) per cycle-day for 18S AVSs of Surface Samples") +
  theme(title = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Metabolites =================================
# Plot all samples separately
ggplot(surface_alpha_diversity_metabolites, aes(x=sample_name, y=Shannon))+
  geom_point() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Alpha diversity (Shannon) for Metabolites of Surface Samples") +
  theme(title = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot average alpha diversity per cycle-day
ggplot(surface_alpha_diversity_metabolites, aes(x=ATTRIBUTE_Location, y=Shannon))+
  geom_boxplot() +
  labs(x = "Samples",
       y = "Shannon Index",
       title = "Average Alpha diversity (Shannon) per cycle-day for Metabolites of Surface Samples") +
  theme(title = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##########################################################################################

##### MDS/PCoA of surface samples using BrayCurtis distance #####
##########################################################################################
# 16s AVSs =================================
# Calculate bray-curtis distance matrix
surface_bray_asv16s <- vegan::vegdist(surface_asv16s, method = "bray")
surface_bray_asv16s

# Calculate PCoA eigen values and vectors
surface_bray_pcoa_asv16s <- ecodist::pco(surface_bray_asv16s)

# Extract data about vectors for the first two components. Add "cycle" variable to color samples in graph.
surface_bray_pcoa_df_asv16s <- data.frame(pcoa1 = surface_bray_pcoa_asv16s$vectors[,1], 
                                  pcoa2 = surface_bray_pcoa_asv16s$vectors[,2],
                                  cycle = metadata_surface$ATTRIBUTE_Filament_Possition)

# Graph PCoA
bray_pcoa_asv16s <- ggplot(data = surface_bray_pcoa_df_asv16s, aes(x=pcoa1, y=pcoa2, label = surface_alpha_diversity_asv16s$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of surface samples using BrayCurtis distance on 16S AVSs") +
  theme(title = element_text(size = 8))

bray_pcoa_asv16s

# 18s AVSs =================================
# Calculate bray-curtis distance matrix
surface_bray_asv18s <- vegan::vegdist(surface_asv18s, method = "bray")
surface_bray_asv18s

# Calculate PCoA eigen values and vectors
surface_bray_pcoa_asv18s <- ecodist::pco(surface_bray_asv18s)

# Extract data about vectors for the first two components. Add cycle variable to color samples in graph.
surface_bray_pcoa_df_asv18s <- data.frame(pcoa1 = surface_bray_pcoa_asv18s$vectors[,1], 
                                   pcoa2 = surface_bray_pcoa_asv18s$vectors[,2],
                                   cycle = metadata_surface$ATTRIBUTE_Filament_Possition)

# Graph PCoA
bray_pcoa_asv18s <- ggplot(data = surface_bray_pcoa_df_asv18s, aes(x=pcoa1, y=pcoa2, label = surface_alpha_diversity_asv18s$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of surface samples using BrayCurtis distance on 18S AVSs") +
  theme(title = element_text(size = 8))

bray_pcoa_asv18s

# Metabolites =================================
# Calculate bray-curtis distance matrix
surface_bray_metabolites <- vegan::vegdist(surface_metabolites, method = "bray")
surface_bray_metabolites

# Calculate PCoA eigen values and vectors
surface_bray_pcoa_metabolites <- ecodist::pco(surface_bray_metabolites)

# Extract data about vectors for the first two components. Add cycle variable to color samples in graph.
surface_bray_pcoa_df_metabolites <- data.frame(pcoa1 = surface_bray_pcoa_metabolites$vectors[,1], 
                                          pcoa2 = surface_bray_pcoa_metabolites$vectors[,2],
                                          cycle = metadata_surface$ATTRIBUTE_Filament_Possition)

# Graph PCoA
bray_pcoa_metabolites <- ggplot(data = surface_bray_pcoa_df_metabolites, aes(x=pcoa1, y=pcoa2, label = metadata_surface$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of surface samples using BrayCurtis distance on metabolites") +
  theme(title = element_text(size = 8))

bray_pcoa_metabolites


##########################################################################################

##### MDS/PCoA of all samples using BrayCurtis distance #####
##########################################################################################
# 16s AVSs =================================
# Calculate bray-curtis distance matrix
all_bray_16s <- vegan::vegdist(asv_16s, method = "bray")
all_bray_16s

# Calculate PCoA eigen values and vectors
all_bray_pcoa_16s <- ecodist::pco(all_bray_16s)

# Extract data about vectors for the first two components. Add cycle variable to color samples in graph.
all_bray_pcoa_df_16s <- data.frame(pcoa1 = all_bray_pcoa_16s$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_16s$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

# Graph PCoA
all_bray_plot_16s <- ggplot(data = all_bray_pcoa_df_16s, aes(x=pcoa1, y=pcoa2, label = metadata$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of all samples using BrayCurtis distance on 16S AVSs") +
  theme(title = element_text(size = 8))

all_bray_plot_16s

# 18s AVSs =================================
# Calculate bray-curtis distance matrix
all_bray_18s <- vegan::vegdist(asv_18s, method = "bray")
all_bray_18s

# Calculate PCoA eigen values and vectors
all_bray_pcoa_18s <- ecodist::pco(all_bray_18s)

# Extract data about vectors for the first two components. Add cycle variable to color samples in graph.
all_bray_pcoa_df_18s <- data.frame(pcoa1 = all_bray_pcoa_18s$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_18s$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

# Graph PCoA
all_bray_plot_18s <- ggplot(data = all_bray_pcoa_df_18s, aes(x=pcoa1, y=pcoa2, label = metadata$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of all samples using BrayCurtis distance on 18S AVSs") +
  theme(title = element_text(size = 8))

all_bray_plot_18s

# Metabolites =================================
# Extract only samples present in ASVs tables
metabolites <- metabolites[row.names(asv_16s),]

# Calculate bray-curtis distance matrix
all_bray_metabolites <- vegan::vegdist(metabolites, method = "bray")
all_bray_metabolites

# Calculate PCoA eigen values and vectors
all_bray_pcoa_metabolites <- ecodist::pco(all_bray_metabolites)

# Extract data about vectors for the first two components. Add cycle variable to color samples in graph.
all_bray_pcoa_df_metabolites <- data.frame(pcoa1 = all_bray_pcoa_metabolites$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_metabolites$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

# Graph PCoA
all_bray_plot_metabolites <- ggplot(data = all_bray_pcoa_df_metabolites, aes(x=pcoa1, y=pcoa2, label = metadata$sample_name, colour = cycle)) +
  geom_point() +
  geom_text(size = 2, nudge_x = 0.005, nudge_y = 0.005, angle = 45) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = 0) +
  labs(x = "PCoA1",
       y = "PCoA2", 
       title = "PCoA of all samples using BrayCurtis distance on metabolites") +
  theme(title = element_text(size = 8))

all_bray_plot_metabolites

##########################################################################################

##### Correlation between 16s AVSs and metabolites #####
##########################################################################################
# First lets collapse our ASVs table using only Genus level counts

# Read our ASVs data
ASV_16S <- read_csv("ASV_16S.csv")

# Remove IDs column. This cannot be collapsed because it is type char.
ASV_16S_genus <- ASV_16S[, 2:(ncol(ASV_16S)-8)]

# Lets put genus as first column that serve as IDs.
# At this point several ASVs will have repeated Genus IDs
ASV_16S_genus <- cbind(ASV_16S["Genus16S"], ASV_16S_genus)

# Now lets use dplyr to collapse by (group_by) Genus IDs.
# ASVs with repeated Genus IDs will have their counts summed up.
ASV_16S_genus <- ASV_16S_genus %>%
  group_by(Genus16S) %>%
  summarise_all(funs(sum)) #"summarise_all" will do "sum" to "grouped_by" rows along all columns

# Let's transpose our df to use with the "cor" function. Samples will be rows and ASVs will be cols
ASV_16S_genus <- data.frame(t(ASV_16S_genus[1:157,]))

# Lets set the first row as colnames and remove it
colnames(ASV_16S_genus) <- ASV_16S_genus[1,]
ASV_16S_genus <- ASV_16S_genus[2:nrow(ASV_16S_genus),]

# Convert all columns to numeric, because previous transposing converts them to char
ASV_16S_genus <- convert_all_cols_to_numeric(ASV_16S_genus)

# Extract samples present in ASVs tables in the metabolites DF.
metabolites <- metabolites[row.names(asv_16s),]

# Check that row names in both DFs are the same to prevent correlation errors.
row.names(ASV_16S_genus) == row.names(metabolites)

# Calculate pearson correlations between ASVs16s and metabolites.
cor_table <- cor(ASV_16S_genus, metabolites, method = "pearson")

# Let's just use the first 100 correlations
cor_table_small <- cor_table[0:100, 0:100]

# Plot heatmap
heatmap(cor_table_small)

# Convert correlation matrix into a edges list DF
edges_genus <- data.frame(row=rownames(cor_table_small)[row(cor_table_small)], col=colnames(cor_table_small)[col(cor_table_small)], corr=c(cor_table_small))

# Save edges list
write.csv(edges_genus, "edges_genus.csv", row.names = FALSE)
