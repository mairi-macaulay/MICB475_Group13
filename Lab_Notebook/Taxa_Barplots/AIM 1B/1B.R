### Taxa bar plots for shower recency ###

# Load libraries
library(phyloseq) 
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

# Load Rdata
load("../../Phyloseq/dorms_final_showerrecency.RData")


# Extract OTU data
otu_table = data.frame(t(otu_table(dorms_final)))
otu_table$ID = rownames(otu_table)

# Extract metadata
metadata = data.frame(sample_data(dorms_final))
metadata$ID = rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_final))

# Format the taxa dataframe and clean names
tax_mat <- tax[,-1]
tax_mat = data.frame(tax_mat)
tax_mat$Phylum = gsub("^...","",tax_mat$Phylum)
tax_mat$Class = gsub("^...","",tax_mat$Class)
tax_mat$Order = gsub("^...","",tax_mat$Order)
tax_mat$Family = gsub("^...","",tax_mat$Family)
tax_mat$Genus = gsub("^...","",tax_mat$Genus)
tax_mat$Species = gsub("^...","",tax_mat$Species)
tax_mat$ASV = rownames(tax_mat)

# Join OTU and metadata
otu_meta = inner_join(metadata,otu_table , by = "ID")

# Transform the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we want to exclude.
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27) )

# Join the taxa information to this transformed dataframe
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


# Collecting the list of unique names for the loop (recent or not_recent)
vars = unique(grouped_taxa$last_shower_binned)

#Create an empty dataframe that the loop will fill. Re run this line before re running the loop.
data_rel = data.frame()

# Loop though recent and not_recent shower recency to create a relative abundance measure for each.
for (i in vars){
  df = grouped_taxa %>%
    filter(last_shower_binned == i)
  
  df_sum = df %>%
    group_by(last_shower_binned, Phylum) %>%
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

# Customize the labels for the x axis
custom <- c("not recent", "recent")

# Plot the results at the phylum level using ggplot2.
ggplot(data = data_rel, aes(last_shower_binned,rel_abs, fill = Phylum))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90))+
  scale_x_discrete(labels = custom)+
  labs(y = "Relative abundance", x = "Shower Recency")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -90),
        axis.title = element_text(face = "bold"))
