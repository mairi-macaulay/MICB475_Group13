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

# Formatting the taxa dataframe and cleaning names
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

# Transforming the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

#Join the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

# Create a new column that contains taxa and sex information
grouped_taxa$legend = paste(grouped_taxa$sex) 

# Calculate the relative abundance for each indivudal
ppl <- unique(grouped_taxa$ID) 
data_rel = data.frame() #Create empty dataframe
for (i in ppl){
  df = grouped_taxa %>%
    filter(ID == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex,sheetwashfreq_binned,Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

# Generate taxa bar plots representing data at the phylum level
data_rel$sheetwashfreq_binned = factor(data_rel$sheetwashfreq_binned, levels = c("low","medium","high")) # Create the order for low, medium, high in the plot
ggplot(data =data_rel, aes(ID,rel_abs, fill = Phylum))+ # Generating the plot with X axis equal to individual
  geom_col()+
  labs(y="Relative Abundance", x = "Individuals (ID)")+
  theme(axis.text.x = element_text(angle = -90))+
  facet_grid(cols = vars(sex), scales = "free_x", space = "free_x")




