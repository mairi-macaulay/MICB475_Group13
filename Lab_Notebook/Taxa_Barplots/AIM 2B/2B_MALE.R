library(phyloseq) 
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

#load Rdata
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_male.RData")


#Extracting OTU data
otu_table = data.frame(t(otu_table(dorms_rare_sheetwashfreq_male)))
otu_table$ID = rownames(otu_table)

#Extracting metadata
metadata = data.frame(sample_data(dorms_rare_sheetwashfreq_male))
metadata$ID = rownames(metadata)

#load the raw taxonomy file
tax <- data.frame(tax_table(dorms_rare_sheetwashfreq_male))

#Formatting the taxa dataframe and cleaning names
tax_mat <- tax[,-1]
tax_mat = data.frame(tax_mat)
tax_mat$Phylum = gsub("^...","",tax_mat$Phylum)
tax_mat$Class = gsub("^...","",tax_mat$Class)
tax_mat$Order = gsub("^...","",tax_mat$Order)
tax_mat$Family = gsub("^...","",tax_mat$Family)
tax_mat$Genus = gsub("^...","",tax_mat$Genus)
tax_mat$Species = gsub("^...","",tax_mat$Species)
tax_mat$ASV = rownames(tax_mat)

#joining OTU and metadata
otu_meta = inner_join(metadata,otu_table , by = "ID")

#Transforming the OTU matrix to a single column called abundance. 53 represents the number of metadata columns we want to exclude.
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:53) )

#joining the taxa information to this transformed dataframe
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


#collecting the list of unique severity names for the loop (high, medium, low)
vars = unique(grouped_taxa$sheetwashfreq_binned)

#Creating an empty dataframe that the loop will fill. Re run this line before re running the loop.
data_rel = data.frame()

#looping though each severity index to create a relative abundance measure.
for (i in vars){
  #vars = "normal BM.Soup.Broth"
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum = df %>%
    group_by(sheetwashfreq_binned, Phylum) %>%
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

#plotting the results at the phylum level
ggplot(data = data_rel, aes(sheetwashfreq_binned,rel_abs, fill = Phylum))+
  geom_col(color = "black")+
  theme(axis.text.x = element_text(angle = -90))+
  labs(y = "Relative abundance", x = "Sheet Wash Frequency")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 15,face = "bold"))

