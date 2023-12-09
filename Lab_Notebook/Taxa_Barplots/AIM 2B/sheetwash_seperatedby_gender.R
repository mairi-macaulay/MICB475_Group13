# Loading libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)


# Load Rdata (phyloseq object)
load("../../Phyloseq/dorms_final_sheetwashfreq.RData")



# Extracting OTU data from the phyloseq object
otu_table <- data.frame(t(otu_table(dorms_final)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata from the phyloseq object
metadata <- data.frame(sample_data(dorms_final))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file 
tax <- data.frame(tax_table(dorms_final))

# Formatting the taxa dataframe and clean names
tax_mat <- tax[,-1]
tax_mat <- data.frame(tax_mat)
tax_mat$Phylum <- gsub("^...","",tax_mat$Phylum)
tax_mat$Class <- gsub("^...","",tax_mat$Class)
tax_mat$Order <- gsub("^...","",tax_mat$Order)
tax_mat$Family <- gsub("^...","",tax_mat$Family)
tax_mat$Genus <- gsub("^...","",tax_mat$Genus)
tax_mat$Species <- gsub("^...","",tax_mat$Species)
tax_mat$ASV <- rownames(tax_mat)

# Join OTU, metadata, and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

# Transform the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

# Join the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

# Generate a column that combined two variables together - sheet washing freqeuncy and taxa, sex and taxa
grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned,grouped_taxa$sex) #Can add gender for when gender becomes a consideration

# Determine which groups to loop through based on the line before ^
levels <- unique(grouped_taxa$legend)

#Generate an empty dataframe
data_rel = data.frame()
#Generate the relative abundance at the phylum level for sheetwashing freq split by sex
for (i in levels){
  
  df = grouped_taxa %>%
    filter(legend == i, sheetwashfreq_binned != "medium")
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

# Combine the timepoints together from the multiple readings from the same individual.
data_rel_sum = data_rel %>%
  group_by(sex,sheetwashfreq_binned,Phylum) %>%
  summarise(mean_rel_abs = sum(rel_abs))

# Filter for phyla that represent relative abundance greater than 1% of the group.
data_rel_sum_filtered = data_rel_sum %>%
  filter(mean_rel_abs> 1)

# Use ggplot 2 to create taxa bar plots representing the average relative abundance for each phylum across the different sheetwash frequency levels.
data_rel_sum_filtered$sheetwashfreq_binned = factor(data_rel_sum_filtered$sheetwashfreq_binned, levels = c("low", "high")) #create the order for low, and high in the plot
ggplot(data =data_rel_sum_filtered, aes(sex,mean_rel_abs, fill = Phylum))+ 
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90),
        axis.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")
        ,legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"))+
  facet_grid(cols = vars(sheetwashfreq_binned), scales = "free_x", space = "free_x")+
  labs(x = "", y = "Relative abundance (%)")


