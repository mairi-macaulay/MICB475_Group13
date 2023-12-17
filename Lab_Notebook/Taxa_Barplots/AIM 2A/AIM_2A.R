### Taxa bar plots of relative abundace of phyla for males and females ###

# Load the libraries required
library(phyloseq) 
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

# Load the Rdata (phyloseq object)
load("../../Phyloseq/dorms_final_sheetwashfreq.RData")


# Extract the OTU data from phyloseq object
otu_table = data.frame(t(otu_table(dorms_final)))
otu_table$ID = rownames(otu_table)

# Extract metadata from phyloseq object
metadata = data.frame(sample_data(dorms_final))
metadata$ID = rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_final))

# Formatting the taxa dataframe and clean names
tax_mat <- tax[,-1]
tax_mat = data.frame(tax_mat)
tax_mat$Phylum = gsub("^...","",tax_mat$Phylum)
tax_mat$Class = gsub("^...","",tax_mat$Class)
tax_mat$Order = gsub("^...","",tax_mat$Order)
tax_mat$Family = gsub("^...","",tax_mat$Family)
tax_mat$Genus = gsub("^...","",tax_mat$Genus)
tax_mat$Species = gsub("^...","",tax_mat$Species)
tax_mat$ASV = rownames(tax_mat)

# Joining OTU and metadata
otu_meta = inner_join(metadata,otu_table , by = "ID")

# Transform the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we want to exclude.
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27) )

# Joining the taxa information to this transformed dataframe
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


# Collect the list of unique sexes for the loop (male, female)
vars = unique(grouped_taxa$sex)

# Creating an empty dataframe that the loop will fill. Re run this line before re running the loop.
data_rel = data.frame()

# Loop though each sex to create a relative abundance measure.
for (i in vars){
  df = grouped_taxa %>%
    filter(sex == i)
  
  df_sum = df %>%
    group_by(sex, Phylum) %>%
    summarize(rel_abs = sum(abundance))
  
  df_sum$rel_abs = as.factor(df_sum$rel_abs)
  row_sub = apply(df_sum, 1, function(row) all(row !=0 ))
  
  df_sum = df_sum[row_sub,]
  df_sum$rel_abs = as.numeric(as.character(df_sum$rel_abs))
  df_sum = na.omit(df_sum)
  count = sum(df_sum$rel_abs)
  df_sum$rel_abs = df_sum$rel_abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  data_rel = na.omit(data_rel)
}

# Using ggplot2, generate taxa bar plots at the phylum level with sex on the x axis.
ggplot(data = data_rel, aes(sex,rel_abs, fill = Phylum))+
  geom_col()+
  theme_bw()+
  labs(y = "Relative abundance", x = "Gender")+
  theme(axis.text.x = element_text(angle = -90),
        axis.title = element_text(face = "bold"))


