# Core Microbiome Aim 1B (Sheet Washing Frequency and Sex, just male)

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(RColorBrewer)
library(reshape2)

# Load data 
load("dorms_final_sheetwashfreq_deseq_male.RData")

# Access the metadata of the phyloseq object
metadata <- sample_data(dorms_final_sheetwashfreq_deseq_male)
sheetwashfreq <- metadata$sheetwashfreq_week

# Making new column "new_sheetwash" that contains only "low" and "high" group
metadata$new_sheetwash <- ifelse(sheetwashfreq > 2, "low", "high")

# Update the sample data in the phyloseq object
sample_data(dorms_final_sheetwashfreq_deseq_male) <- metadata

set.seed(1)

# Convert to relative abundance
dorms_male_RA <- transform_sample_counts(dorms_final_sheetwashfreq_deseq_male, fun=function(x) x/sum(x))

# Filter dataset by sheet washing freq
dorms_male_low <- subset_samples(dorms_male_RA, `new_sheetwash`=="low")
dorms_male_high <- subset_samples(dorms_male_RA, `new_sheetwash`=="high")

# HEATMAP
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

c = plot_core(dorms_male_low, 
              plot.type = "heatmap", 
              colours = rev(brewer.pal(5, "RdBu")),
              prevalences = prevalences, 
              detections = det,
              min.prevalence = 0.75) +
  xlab("Detection Threshold (Relative Abundance (%))")

c

#Sort taxa names from the AVS in the plot
df <- c$data
list <- df$Taxa

#check the ASV ids
print(list)

#get the taxonomy data from the phyloseq object that was used for the plot. This is for low freq
tax <- tax_table(dorms_male_low)
tax <- as.data.frame(tax)

#add the ASVs to last column
tax$OTU <- rownames(tax)

#select taxonomy of only those ASVs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list)

#Make a column with the information from Family and Genus together. This can be changed to look at other taxa levels like "Phylum"
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Family","Genus"), sep = "_;", remove = TRUE)
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

#Now that we have the real names for the ASVs, put this information back into the plot 
df$Taxa <- tax.unit$Taxa_level

#Rerun the plot
#This plot might look different to the original one because it is resolved to the genus level. This means that any ASV's that had an "NA" for the species name is removed. This is fine because we would remove NA's anyways if we were to look at the species level.
dorms_low_male_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_low_male_heat
ggsave("heat_male_low.png", dorms_low_male_heat)

# HEATMAP FOR HIGH_MALE
c = plot_core(dorms_male_high, 
              plot.type = "heatmap", 
              colours = rev(brewer.pal(5, "RdBu")),
              prevalences = prevalences, 
              detections = det,
              min.prevalence = 0.75) +
  xlab("Detection Threshold (Relative Abundance (%))")

c

#Sort taxa names from the AVS in the plot
df <- c$data
list <- df$Taxa

#check the ASV ids
print(list)

#get the taxonomy data from the phyloseq object that was used for the plot. This is for low freq
tax <- tax_table(dorms_male_high)
tax <- as.data.frame(tax)

#add the ASVs to last column
tax$OTU <- rownames(tax)

#select taxonomy of only those ASVs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list)

#Make a column with the information from Family and Genus together. This can be changed to look at other taxa levels like "Phylum"
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Family","Genus"), sep = "_;", remove = TRUE)
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

#Now that we have the real names for the ASVs, put this information back into the plot 
df$Taxa <- tax.unit$Taxa_level

#Rerun the plot
#This plot might look different to the original one because it is resolved to the genus level. This means that any ASV's that had an "NA" for the species name is removed. This is fine because we would remove NA's anyways if we were to look at the species level.
dorms_high_male_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_high_male_heat
ggsave("heat_male_high.png", dorms_high_male_heat)


# Venn Diagram
# Setting detection and prevalence thresholds
low_ASVs <- core_members(dorms_male_low, detection=0.001, prevalence = 0.75)
high_ASVs <- core_members(dorms_male_high, detection=0.001, prevalence = 0.75)

dorms_list_low_high <- list(Low_Frequency = low_ASVs, High_Frequency = high_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to 2 groups
dorms_low_high <- ggVennDiagram(x = dorms_list_low_high)

dorms_low_high

ggsave("venn_dorms_low_high_male.png", dorms_low_high)

#What are these ASVs?
low_male_tax <- tax_table(prune_taxa(low_ASVs,dorms_final_sheetwashfreq_deseq_male))
high_male_tax <- tax_table(prune_taxa(high_ASVs,dorms_final_sheetwashfreq_deseq_male))

#Unique Genera/species for low vs high
tax_male_diff <- setdiff(low_male_tax, high_male_tax)
print(tax_male_diff)

tax_male_diff2 <- setdiff(high_male_tax, low_male_tax)
print(tax_male_diff2)