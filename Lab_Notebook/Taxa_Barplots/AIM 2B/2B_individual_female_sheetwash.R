library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

# Load Rdata
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_female.RData")

# Extracting OTU data
otu_table <- data.frame(t(otu_table(dorms_rare_sheetwashfreq_female)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata
metadata <- data.frame(sample_data(dorms_rare_sheetwashfreq_female))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_rare_sheetwashfreq_female))

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


# Joining OTU and metadata and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

# Transforming the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we wnat to exlcude.
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

# Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) 

#Calculating the relative abundance for each indivudal
ppl <- unique(grouped_taxa$ID) 
data_rel = data.frame() #Create empty dataframe
for (i in ppl){
  df = grouped_taxa %>%
    filter(ID == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

# Create taxa bar plots using ggplot2
data_rel$sheetwashfreq_binned = factor(data_rel$sheetwashfreq_binned, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel, aes(ID,rel_abs, fill = Phylum))+ #Generating the plot with X axis equal to individual
  geom_col()+
  theme(axis.text.x = element_text(angle = -90))+
  labs(y="Relative Abundance", x = "Individuals (ID)")+
  ggtitle("FEMALE")+
  theme(axis.text.x = element_text(size = 4),
        axis.title = element_text(size = 15,face = "bold")) +
  facet_grid(cols = vars(sheetwashfreq_binned), scales = "free_x", space = "free_x")
