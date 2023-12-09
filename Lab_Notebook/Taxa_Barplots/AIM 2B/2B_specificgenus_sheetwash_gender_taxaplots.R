# Loading libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)


# Load the Rdata (Phyloseq object)
load("../../Phyloseq/dorms_final_sheetwashfreq.RData")



# Extracting OTU data from the phyloseq object
otu_table <- data.frame(t(otu_table(dorms_final)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata from the phyloseq object
metadata <- data.frame(sample_data(dorms_final))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(dorms_final))

# Format the taxa dataframe and clean the names
tax_mat <- tax[,-1]
tax_mat <- data.frame(tax_mat)
tax_mat$Phylum <- gsub("^...","",tax_mat$Phylum)
tax_mat$Class <- gsub("^...","",tax_mat$Class)
tax_mat$Order <- gsub("^...","",tax_mat$Order)
tax_mat$Family <- gsub("^...","",tax_mat$Family)
tax_mat$Genus <- gsub("^...","",tax_mat$Genus)
tax_mat$Species <- gsub("^...","",tax_mat$Species)
tax_mat$ASV <- rownames(tax_mat)

# Join the OTU, metadata, and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

# Transform the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we want to exlcude
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

# Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")

#Generate a column that combines two variables together - sheet washing frequency and taxa
grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned,grouped_taxa$sex) #Can add gender for when gender becomes a consideration


#Generating relative abundance at the genus level

#Define the low, medium, and high levels
levels <- unique(grouped_taxa$legend)
#Create a new empty dataframe to put the phylum data in.
data_rel_genus = data.frame()

#Run a loop to generate the relative abundance of each genus for low, medium, and high groups.
for (i in levels){
  
  df = grouped_taxa %>%
    filter(legend == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned,Phylum,Order,Family,Class, Genus) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel_genus = rbind(data_rel_genus, df_sum)
  
}

# Filter for only Actinobacteriota, repeat this step for Firmicutes, Bacteroidota, Fusobacteriota, Proteobacteria by replacing "Actinobacteriota" in the code with the various phyla
data_rel_Actinobacteriota = data_rel_genus %>%
  filter(Phylum == "Actinobacteriota", sheetwashfreq_binned != "medium") %>% #Keep only bacteria that are in the Actinobacteriota phylum and remove medium sheetwashing freq. Change this line for different phyla.
  group_by(legend,Phylum, Order, Class, Family,Genus,sex,sheetwashfreq_binned) %>%
  summarise(mean_rel_abs = sum(rel_abs))%>% # Add relative abundances together for multiple species that have the same genus
  filter(mean_rel_abs>= 1) #Remove Genus that have a relative abundance less than 1%

# Use ggplot2 to generate taxa bar plots of the various phyla of interest - in this case it is for Actinobacteriota
data_rel_Actinobacteriota$sheetwashfreq_binned = factor(data_rel_Actinobacteriota$sheetwashfreq_binned, levels = c("low","high")) # Create the order for low and high in the plot
ggplot(data =data_rel_Actinobacteriota, aes(sex,mean_rel_abs, fill = Genus))+ 
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

