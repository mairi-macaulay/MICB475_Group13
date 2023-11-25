library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)


# Load Rdata
load("../Phyloseq/dorms_final_sheetwashfreq.RData")



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

#transforming the OTU matrix to a single column called abundance. 53 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) #Can add gender for when gender becomes a consideration

#Generating the relative abundance for all individuals within low/medium/high sheetwashing freq.
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

data_rel_sum = data_rel %>%
  group_by(legend,Phylum) %>%
  summarise(mean_rel_abs = sum(rel_abs))

#Filtering for phyla that represent relative abundance greater than 1% of the group.
data_rel_sum_filtered = data_rel_sum %>%
  filter(mean_rel_abs> 1)

#This plot represents the average relative abundance for each phylum across the different sheetwash frequency levels.
data_rel_sum_filtered$legend = factor(data_rel_sum_filtered$legend, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel_sum_filtered, aes(legend,mean_rel_abs, fill = Phylum))+#Generating the plot with X axis equal to sheetwash_freq_binned
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))

#####################################################################
#NEW CODE FOR PLOTS ABOUT SPECIFIC PHYLUM
####################################################################

#Generating relative abundance at the genus level

#Define the low, medium, and high levels
levels <- unique(grouped_taxa$sheetwashfreq_binned)
#Create a new empty dataframe to put firmicute data in.
data_rel_genus = data.frame()

#Run a loop to generate the relative abundance of each genus for low, medium, and high groups.
for (i in levels){
  
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned,Phylum,Order,Family,Class, Genus) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel_genus = rbind(data_rel_genus, df_sum)
  
}

#Filter for only Proteobacteria 
data_rel_proteobacteria = data_rel_genus %>%
  filter(Phylum == "Bacteroidota") %>% #Keep only bacteria that are in the Proteobacteria phylum. Change this line for different phyla.
  group_by(legend,Phylum, Order, Class, Family,Genus) %>%
  summarise(mean_rel_abs = sum(rel_abs))%>% #Add relative abundances together for multiple species that have the same genus
  filter(mean_rel_abs>= 1) #Remove Genus that have a relative abundance less than 1%

data_rel_proteobacteria$legend = factor(data_rel_proteobacteria$legend, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel_proteobacteria, aes(legend,mean_rel_abs, fill = Family))+#Generating the plot with X axis equal to sheetwash_freq_binned
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))

#This plot shows that Proteobacteria abundance increases with sheet wash frequency. You can also see the trends for different families.

