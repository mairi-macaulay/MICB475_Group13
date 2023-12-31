### Aim 1A - removing medium from taxa bar plots for sheet washing frequency ###

# Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)


# Load Rdata (phyloseq object)
load("../../Phyloseq/dorms_final_sheetwashfreq.RData")



# Extract OTU data from the phyloseq object
otu_table <- data.frame(t(otu_table(dorms_final)))
otu_table$ID <- rownames(otu_table)

# Extract metadata from the phyloseq object
metadata <- data.frame(sample_data(dorms_final))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_final))

# Format the taxa dataframe and clean names
tax_mat <- tax[,-1]
tax_mat <- data.frame(tax_mat)
tax_mat$Phylum <- gsub("^...","",tax_mat$Phylum)
tax_mat$Class <- gsub("^...","",tax_mat$Class)
tax_mat$Order <- gsub("^...","",tax_mat$Order)
tax_mat$Family <- gsub("^...","",tax_mat$Family)
tax_mat$Genus <- gsub("^...","",tax_mat$Genus)
tax_mat$Species <- gsub("^...","",tax_mat$Species)
tax_mat$ASV <- rownames(tax_mat)

# Join OTU and metadata and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

# Transform the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

# Join the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

# Generate a new column containing taxa and sheet washing frequency data 
grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) 

#Generating the relative abundance for all individuals within low/high sheetwashing freq. Filter out "medium" sheet washing frequency.
levels <- unique(grouped_taxa$sheetwashfreq_binned)
data_rel = data.frame()
for (i in levels){
  
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned == i, sheetwashfreq_binned != "medium")
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

data_rel_sum = data_rel %>%
  group_by(legend,Phylum) %>%
  summarise(mean_rel_abs = sum(rel_abs))

#Filter for phyla that represent relative abundance greater than 1% of the group.

data_rel_sum_filtered = data_rel_sum %>%
  filter(mean_rel_abs> 1)

# Generate taxa bar plots representing the average relative abundance for each phylum across the different sheetwash frequency levels.
data_rel_sum_filtered$legend = factor(data_rel_sum_filtered$legend, levels = c("low","high")) # Create the order for low, high in the plot
ggplot(data =data_rel_sum_filtered, aes(legend,mean_rel_abs, fill = Phylum))+ # Generating the plot with X axis equal to sheetwash_freq_binned
  labs(y="Relative Abundance", x = "Sheet Washing Freqeuncy") +
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))
