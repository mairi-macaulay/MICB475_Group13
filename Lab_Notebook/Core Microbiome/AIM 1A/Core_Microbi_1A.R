# Core Microbiome Aim 1A (Sheet Washing Frequency)

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load data 
load("dorms_final_sheetwashfreq.RData")

set.seed(1)

# "core" microbiome

# Convert to relative abundance
dorms_sheetwash_RA <- transform_sample_counts(dorms_final, fun=function(x) x/sum(x))

# Filter dataset by sheet washing freq
dorms_sheetwash_low <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="low")
dorms_sheetwash_med <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="medium")
dorms_sheetwash_high <- subset_samples(dorms_sheetwash_RA, `sheetwashfreq_binned`=="high")

# Setting detection and prevalence thresholds
low_ASVs <- core_members(dorms_sheetwash_low, detection=0.001, prevalence = 0.75)
med_ASVs <- core_members(dorms_sheetwash_med, detection=0.001, prevalence = 0.75)
high_ASVs <- core_members(dorms_sheetwash_high, detection=0.001, prevalence = 0.75)

dorms_list_low_med <- list(Low_Frequency = low_ASVs, Medium_Frequency = med_ASVs)
dorms_list_low_high <- list(Low_Frequency = low_ASVs, High_Frequency = high_ASVs)
dorms_list_med_high <- list(Medium_Frequency = med_ASVs, High_Frequency = high_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to 2 groups
dorms_low_med <- ggVennDiagram(x = dorms_list_low_med)
dorms_low_high <- ggVennDiagram(x = dorms_list_low_high)
dorms_med_high <- ggVennDiagram(x = dorms_list_med_high)

dorms_low_med
dorms_low_high
dorms_med_high

ggsave("venn_dorms_low_med.png", dorms_low_med)
ggsave("venn_dorms_low_high.png", dorms_low_high)
ggsave("venn_dorms_med_high.png", dorms_med_high)



#What are these ASVs?
low_tax <- tax_table(prune_taxa(low_ASVs,dorms_final))
med_tax <- tax_table(prune_taxa(med_ASVs,dorms_final))
high_tax <- tax_table(prune_taxa(high_ASVs,dorms_final))

#Unique Genera/species for low vs med
tax_diff <- setdiff(low_tax, med_tax)
print(tax_diff)

tax_diff2 <- setdiff(med_tax, low_tax)
print(tax_diff2)

shared_taxa <- union(low_tax, med_tax)
print(shared_taxa)

#Unique Genera/species for low vs high
tax_diff3 <- setdiff(low_tax, high_tax)
print(tax_diff3)

tax_diff4 <- setdiff(high_tax, low_tax)
print(tax_diff4)

#Unique Genera/species for med vs high
tax_diff5 <- setdiff(med_tax, high_tax)
print(tax_diff5)

tax_diff6 <- setdiff(high_tax, med_tax)
print(tax_diff6)

#Plot those ASVs' relative abundance
prune_taxa(low_ASVs,dorms_sheetwash_RA) %>%
  plot_bar(fill="Genus")+
  facet_wrap(.~sheetwashfreq_binned, scales="free")




