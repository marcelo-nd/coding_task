library(dplyr)
library(vegan)
library(ggplot2)

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
surface_bray_asv16s <- vegan::vegdist(surface_asv16s, method = "bray")

surface_bray_asv16s

surface_bray_pcoa_asv16s <- ecodist::pco(surface_bray_asv16s)

surface_bray_pcoa_df_asv16s <- data.frame(pcoa1 = surface_bray_pcoa_asv16s$vectors[,1], 
                                  pcoa2 = surface_bray_pcoa_asv16s$vectors[,2],
                                  cycle = metadata_surface$ATTRIBUTE_Filament_Possition)
  
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
surface_bray_asv18s <- vegan::vegdist(surface_asv18s, method = "bray")

surface_bray_asv18s

surface_bray_pcoa_asv18s <- ecodist::pco(surface_bray_asv18s)

surface_bray_pcoa_df_asv18s <- data.frame(pcoa1 = surface_bray_pcoa_asv18s$vectors[,1], 
                                   pcoa2 = surface_bray_pcoa_asv18s$vectors[,2],
                                   cycle = metadata_surface$ATTRIBUTE_Filament_Possition)

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
surface_bray_metabolites <- vegan::vegdist(surface_metabolites, method = "bray")

surface_bray_metabolites

surface_bray_pcoa_metabolites <- ecodist::pco(surface_bray_metabolites)

surface_bray_pcoa_df_metabolites <- data.frame(pcoa1 = surface_bray_pcoa_metabolites$vectors[,1], 
                                          pcoa2 = surface_bray_pcoa_metabolites$vectors[,2],
                                          cycle = metadata_surface$ATTRIBUTE_Filament_Possition)

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
all_bray_16s <- vegan::vegdist(asv_16s, method = "bray")

all_bray_16s

all_bray_pcoa_16s <- ecodist::pco(all_bray_16s)

all_bray_pcoa_df_16s <- data.frame(pcoa1 = all_bray_pcoa_16s$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_16s$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

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
all_bray_18s <- vegan::vegdist(asv_18s, method = "bray")

all_bray_18s

all_bray_pcoa_18s <- ecodist::pco(all_bray_18s)

all_bray_pcoa_df_18s <- data.frame(pcoa1 = all_bray_pcoa_18s$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_18s$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

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
# Extract samples present in ASVs tables
metabolites <- metabolites[row.names(asv_16s),]

all_bray_metabolites <- vegan::vegdist(metabolites, method = "bray")

all_bray_metabolites

all_bray_pcoa_metabolites <- ecodist::pco(all_bray_metabolites)

all_bray_pcoa_df_metabolites <- data.frame(pcoa1 = all_bray_pcoa_metabolites$vectors[,1], 
                                   pcoa2 = all_bray_pcoa_metabolites$vectors[,2],
                                   cycle = metadata$ATTRIBUTE_Filament_Possition)

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
