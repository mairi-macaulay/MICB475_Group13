# Core Microbiome Aim 1B (Sheet Washing Frequency and Sex, NO MEDIUM GROUP, just female)

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(RColorBrewer)
library(reshape2)

# Load data 
load("dorms_final_sheetwashfreq_deseq_female.RData")

set.seed(1)

# Convert to relative abundance
dorms_female_RA <- transform_sample_counts(dorms_final_sheetwashfreq_deseq_female, fun=function(x) x/sum(x))

# Filter dataset by sheet washing freq
dorms_female_low <- subset_samples(dorms_female_RA, `sheetwashfreq_binned`=="low")
dorms_female_high <- subset_samples(dorms_female_RA, `sheetwashfreq_binned`=="high")

metadata_of_female_high <- sample_data(dorms_female_high)
metadata_of_female_low <- sample_data(dorms_female_low)

# HEATMAP
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

c = plot_core(dorms_female_low, 
              plot.type = "heatmap", 
              colours = rev(brewer.pal(5, "RdBu")),
              prevalences = prevalences, 
              detections = det,
              min.prevalence = 0.5) +
  xlab("Detection Threshold (Relative Abundance (%))")

c

#Sort taxa names from the AVS in the plot
df <- c$data
list <- df$Taxa

#check the ASV ids
print(list)

#get the taxonomy data from the phyloseq object that was used for the plot. This is for low freq
tax <- tax_table(dorms_female_low)
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
dorms_low_female_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_low_female_heat
ggsave("heat_female_NM_low.png", dorms_low_female_heat)

# HEATMAP FOR HIGH_FEMALE
c = plot_core(dorms_female_high, 
              plot.type = "heatmap", 
              colours = rev(brewer.pal(5, "RdBu")),
              prevalences = prevalences, 
              detections = det,
              min.prevalence = 0.5) +
  xlab("Detection Threshold (Relative Abundance (%))")

c

#Sort taxa names from the AVS in the plot
df <- c$data
list <- df$Taxa

#check the ASV ids
print(list)

#get the taxonomy data from the phyloseq object that was used for the plot. This is for low freq
tax <- tax_table(dorms_female_high)
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
dorms_high_female_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_high_female_heat
ggsave("heat_female_NM_high.png", dorms_high_female_heat)


# Venn Diagram --> not needed because doing 4 way 
# Setting detection and prevalence thresholds
low_female_ASVs <- core_members(dorms_female_low, detection=0.001, prevalence = 0.5)
high_female_ASVs <- core_members(dorms_female_high, detection=0.001, prevalence = 0.5)

print(low_ASVs)
print(high_ASVs)
dorms_list_low_high <- list(Low_Frequency = low_female_ASVs, High_Frequency = high_female_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to 2 groups
dorms_low_high <- ggVennDiagram(x = dorms_list_low_high)

dorms_low_high

ggsave("venn_dorms_low_high_NM_female.png", dorms_low_high)



#What are these ASVs?
low_female_tax <- tax_table(prune_taxa(low_ASVs,dorms_final_sheetwashfreq_deseq_female))
high_female_tax <- tax_table(prune_taxa(high_ASVs,dorms_final_sheetwashfreq_deseq_female))

#Unique Genera/species for low vs high
tax_female_diff <- setdiff(low_female_tax, high_female_tax)
print(tax_female_diff)

tax_female_diff2 <- setdiff(high_female_tax, low_female_tax)
print(tax_female_diff2)
