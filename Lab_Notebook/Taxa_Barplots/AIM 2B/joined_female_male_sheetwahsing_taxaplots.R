library(phyloseq) 
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)
library(gridExtra)



#load Rdata
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_male.RData")


#Extracting OTU data
otu_table_male = data.frame(t(otu_table(dorms_rare_sheetwashfreq_male)))
otu_table_male$ID = rownames(otu_table_male)

#Extracting metadata
metadata_male = data.frame(sample_data(dorms_rare_sheetwashfreq_male))
metadata_male$ID = rownames(metadata_male)

#load the raw taxonomy file
tax_male <- data.frame(tax_table(dorms_rare_sheetwashfreq_male))

#Formatting the taxa dataframe and cleaning names
tax_mat_male <- tax[,-1]
tax_mat_male = data.frame(tax_mat_male)
tax_mat_male$Phylum = gsub("^...","",tax_mat_male$Phylum)
tax_mat_male$Class = gsub("^...","",tax_mat_male$Class)
tax_mat_male$Order = gsub("^...","",tax_mat_male$Order)
tax_mat_male$Family = gsub("^...","",tax_mat_male$Family)
tax_mat_male$Genus = gsub("^...","",tax_mat_male$Genus)
tax_mat_male$Species = gsub("^...","",tax_mat_male$Species)
tax_mat_male$ASV = rownames(tax_mat_male)

#joining OTU and metadata
otu_meta_male = inner_join(metadata_male,otu_table_male , by = "ID")

#Transforming the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we want to exclude.
grouped_male = gather(otu_meta_male, key = "ASV", value = "abundance", -(1:27) )

#joining the taxa information to this transformed dataframe
grouped_taxa_male = inner_join(tax_mat_male, grouped_male, by = "ASV", multiple = "all")


#collecting the list of unique severity names for the loop (high, medium, low)
vars_male = unique(grouped_taxa_male$sheetwashfreq_binned)

#Creating an empty dataframe that the loop will fill. Re run this line before re running the loop.
data_rel_male = data.frame()

#looping though each severity index to create a relative abundance measure.
for (i in vars_male){
  #vars = "normal BM.Soup.Broth"
  df_male = grouped_taxa_male %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum_male = df_male %>%
    group_by(sheetwashfreq_binned, Phylum) %>%
    summarize(rel_abs_male = sum(abundance))
  
  df_sum_male$rel_abs_male = as.factor(df_sum_male$rel_abs_male)
  row_sub_male = apply(df_sum_male, 1, function(row) all(row !=0 ))
  
  df_sum_male = df_sum_male[row_sub_male,]
  df_sum_male$rel_abs_male = as.numeric(as.character(df_sum_male$rel_abs_male))
  df_sum_male = na.omit(df_sum_male)
  count_male = sum(df_sum_male$rel_abs_male)
  df_sum_male$rel_abs_male = df_sum_male$rel_abs_male/count*100
  
  data_rel_male = rbind(data_rel_male, df_sum_male)
  data_rel_male = na.omit(data_rel_male)
}

#Store male data in a seperate dataframe with gender information
data_rel_male$Gender <- "Male"


#plotting the results at the phylum level
#data_rel_male$sheetwashfreq_binned = factor(data_rel_male$sheetwashfreq_binned, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
#plot_male <- ggplot(data = data_rel_male, aes(sheetwashfreq_binned,rel_abs_male, fill = Phylum))+
#  theme(axis.text.x = element_text(angle = -90))+
#  labs(y = "Relative abundance", x = "Sheet Wash Frequency")+
#  geom_col()+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle = -90),
#        axis.title = element_text(size = 12,face = "bold"),
#        title = element_text(size = 14),
#        legend.position = "none")














#load Rdata
load("../../Phyloseq/AIM_2B_phyloseq/dorms_rare_sheetwashfreq_female.RData")


#Extracting OTU data
otu_table_female = data.frame(t(otu_table(dorms_rare_sheetwashfreq_female)))
otu_table_female$ID = rownames(otu_table_female)

#Extracting metadata
metadata_female = data.frame(sample_data(dorms_rare_sheetwashfreq_female))
metadata_female$ID = rownames(metadata_female)

#load the raw taxonomy file
tax_female <- data.frame(tax_table(dorms_rare_sheetwashfreq_female))

#Formatting the taxa dataframe and cleaning names
tax_mat_female <- tax_female[,-1]
tax_mat_female = data.frame(tax_mat_female)
tax_mat_female$Phylum = gsub("^...","",tax_mat_female$Phylum)
tax_mat_female$Class = gsub("^...","",tax_mat_female$Class)
tax_mat_female$Order = gsub("^...","",tax_mat_female$Order)
tax_mat_female$Family = gsub("^...","",tax_mat_female$Family)
tax_mat_female$Genus = gsub("^...","",tax_mat_female$Genus)
tax_mat_female$Species = gsub("^...","",tax_mat_female$Species)
tax_mat_female$ASV = rownames(tax_mat_female)

#joining OTU and metadata
otu_meta_female = inner_join(metadata_female,otu_table_female , by = "ID")

#Transforming the OTU matrix to a single column called abundance. 27 represents the number of metadata columns we want to exclude.
grouped_female = gather(otu_meta_female, key = "ASV", value = "abundance", -(1:27) )

#joining the taxa information to this transformed dataframe
grouped_taxa_female = inner_join(tax_mat_female, grouped_female, by = "ASV", multiple = "all")


#collecting the list of unique severity names for the loop (high, medium, low)
vars_female = unique(grouped_taxa_female$sheetwashfreq_binned)

#Creating an empty dataframe that the loop will fill. Re run this line before re running the loop.
data_rel_female = data.frame()

#looping though each severity index to create a relative abundance measure.
for (i in vars_female){
  #vars = "normal BM.Soup.Broth"
  df_female = grouped_taxa_female %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum_female = df_female %>%
    group_by(sheetwashfreq_binned, Phylum) %>%
    summarize(rel_abs_female = sum(abundance))
  
  df_sum_female$rel_abs_female = as.factor(df_sum_female$rel_abs_female)
  row_sub_female = apply(df_sum_female, 1, function(row) all(row !=0 ))
  
  df_sum_female = df_sum_female[row_sub_female,]
  df_sum_female$rel_abs_female = as.numeric(as.character(df_sum_female$rel_abs_female))
  df_sum_female = na.omit(df_sum_female)
  count = sum(df_sum_female$rel_abs_female)
  df_sum_female$rel_abs_female = df_sum_female$rel_abs_female/count*100
  
  data_rel_female = rbind(data_rel_female, df_sum_female)
  data_rel_female = na.omit(data_rel_female)
}

#Store male data in a seperate dataframe with gender information
data_rel_female$Gender <- "Female"



#combine the dataframes
combined_data_rel <- rbind(data_rel_female, data_rel_male)


#plot the taxa bar graph
ggplot(data = combined_data_rel, aes(x = sheetwashfreq_binned, y = rel_abs_male+rel_abs_female, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative abundance", x = "Sheet Wash Frequency") +
  scale_fill_discrete(name = "Phylum") +
  facet_grid(. ~ Gender, scales = "free_x", space = "free_x") +  # Faceting based on Gender
  theme(axis.text.x = element_text(angle = -90),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.placement = "outside")



