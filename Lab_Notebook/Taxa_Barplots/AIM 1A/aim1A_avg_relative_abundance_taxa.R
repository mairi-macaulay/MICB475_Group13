### Taxa bar plots for sheet washing frequency - Relative abundance of phyla ###

# Load libraries
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

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

# Generate a new column containing sheet washing and taxa information
grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) 




#Generate the relative abundance for all individuals within low/medium/high sheetwashing freq.
levels <- unique(grouped_taxa$sheetwashfreq_binned)
data_rel = data.frame()
for (i in levels){
  
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

# Generate taxa bar plots at the phylum level using ggplot2 with an x axis containing sheet washing frequency data. 
data_rel$sheetwashfreq_binned = factor(data_rel$sheetwashfreq_binned, levels = c("low","medium","high")) # Create the order for low, medium, high in the plot
ggplot(data =data_rel, aes(sheetwashfreq_binned,rel_abs, fill = Phylum))+
  geom_col()+
  theme_bw()+
  labs(y="Relative Abundance", x = "Sheet Wash Frequency")+
  theme(axis.text.x = element_text(angle = -90),
        axis.title = element_text(face = "bold"))

