####################################################
# TAXABAR PLOTS SPLITTING SHEETWASH FREQ. BY SEX
###################################################


#Loading libraries
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

#transforming the OTU matrix to a single column called abundance. 53 represents the number of metadata columns we wnat to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

#Generating a column that combined two variables together
grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned,grouped_taxa$sex) #Can add gender for when gender becomes a consideration

#Determining which groups to loop though based on the line before ^
levels <- unique(grouped_taxa$legend)

#Generating an empty dataframe
data_rel = data.frame()
#Generating the relative abundance at the phylum level for sheetwashing freq split by sex
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

#There are multiple readings from the same individual at different timepoints. This bit combines the timepoints together. I'm realizing we should have thought about this before. Will probably come back into discussion.
data_rel_sum = data_rel %>%
  group_by(sex,sheetwashfreq_binned,Phylum) %>%
  summarise(mean_rel_abs = sum(rel_abs))

#Filtering for phyla that represent relative abundance greater than 1% of the group.
data_rel_sum_filtered = data_rel_sum %>%
  filter(mean_rel_abs> 1)

#This plot represents the average relative abundance for each phylum across the different sheetwash frequency levels.
data_rel_sum_filtered$sheetwashfreq_binned = factor(data_rel_sum_filtered$sheetwashfreq_binned, levels = c("low","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel_sum_filtered, aes(sex,mean_rel_abs, fill = Phylum))+#Generating the plot with X axis equal to sheetwash_freq_binned
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


