# Core Microbiome #

#### Preparations ####

# Load Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(RColorBrewer)
library(reshape2)

# Load data (NOT rare file)
load("dorms_final_sheetwashfreq.RData")

# Convert to relative abundance
dorms_sheetwash_RA <- transform_sample_counts(dorms_final, fun=function(x) x/sum(x))

# Filter dataset by sheet washing freq
dorms_sheetwash_low <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="low")
dorms_sheetwash_med <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="medium")
dorms_sheetwash_high <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="high")

# View filtered tables
view(sample_data(dorms_sheetwash_low))
view(sample_data(dorms_sheetwash_med))
view(sample_data(dorms_sheetwash_low))

# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)


c = plot_core(dorms_sheetwash_low, 
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
tax <- tax_table(dorms_sheetwash_low)
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

df

#Rerun the plot
#This plot might look different to the original one because it is resolved to the genus level. This means that any ASV's that had an "NA" for the species name is removed. This is fine because we would remove NA's anyways if we were to look at the species level.
dorms_low_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_low_heat
ggsave("heat_low.png", dorms_low_heat)

#0.75 prevalence
#0.001 detection

#Plot (medium)
c = plot_core(dorms_sheetwash_med, 
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

#get the taxonomy data from the phyloseq object that was used for the plot. This is for med freq
tax <- tax_table(dorms_sheetwash_med)
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
dorms_med_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_med_heat
ggsave("heat_med.png", dorms_med_heat)

#Plot (high)
c = plot_core(dorms_sheetwash_high, 
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

#get the taxonomy data from the phyloseq object that was used for the plot. This is for med freq
tax <- tax_table(dorms_sheetwash_high)
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
dorms_high_heat <- ggplot(df, aes(x = DetectionThreshold, reorder(Taxa, Prevalence), fill=Prevalence))+
  geom_tile()+
  xlab("Relative Abundance")+
  scale_fill_gradientn(colours = rev(brewer.pal(5, "RdBu")))+
  ylab("Genus")

dorms_high_heat

ggsave("heat_high.png", dorms_high_heat)

