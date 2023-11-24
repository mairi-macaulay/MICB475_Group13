library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)


# Load Rdata
load("../../Phyloseq/dorms_final_sheetwashfreq.RData")



# Extracting OTU data
otu_table <- data.frame(t(otu_table(dorms_final)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata
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

# Joining OTU and metadata and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

#transforming the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) #Can add gender for when gender becomes a consideration
#Generating relative abundance at the genus level

#Define the low, medium, and high levels
levels <- unique(grouped_taxa$sheetwashfreq_binned)
#Create a new empty dataframe to put firmicute data in.
data_rel_genus = data.frame()

#Run a loop to generate the relative abundance of each genus for low, medium, and high groups.
for (i in levels){
  
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned != "medium")
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned,Phylum,Order,Family,Class, Genus) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel_genus = rbind(data_rel_genus, df_sum)
  
}

#Filter for only Proteobacteria - repeat this same step for the other phylum but replace Proteobacteria with Firmicutues, Actinobacteriota, Fusobacteriota, and Bacteroidota
data_rel_proteobacteria = data_rel_genus %>%
  filter(Phylum == "Proteobacteria") %>% #Keep only bacteria that are in the Proteobacteria phylum. Change this line for different phyla.
  group_by(legend,Phylum, Order, Class, Family,Genus) %>%
  summarise(mean_rel_abs = sum(rel_abs))%>% #Add relative abundances together for multiple species that have the same genus
  filter(mean_rel_abs>= 1) #Remove Genus that have a relative abundance less than 1%

data_rel_proteobacteria$legend = factor(data_rel_proteobacteria$legend, levels = c("low","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel_proteobacteria, aes(legend,mean_rel_abs, fill = Family))+#Generating the plot with X axis equal to sheetwash_freq_binned
  labs(x = "Relative Abundance", y = "Sheet Washing Frequency")+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))
