# Core Microbiome Aim 1B (Sheet Washing Frequency and Sex)
# Comparing male high sheetwashing freq vs female high sheetwashing freq

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load data 
load("dorms_final_sheetwashfreq.RData")

# Access the metadata of the phyloseq object
metadata <- sample_data(dorms_final)
sheetwashfreq <- metadata$sheetwashfreq_week

# Making new column "new_sheetwash" that contains only "low" and "high" group
metadata$new_sheetwash <- ifelse(sheetwashfreq > 2, "low", "high")

# Update the sample data in the phyloseq object
sample_data(dorms_final) <- metadata

set.seed(1)

# Convert to relative abundance
dorms_gender_RA <- transform_sample_counts(dorms_final, fun=function(x) x/sum(x))

# Filter dataset with just high group
dorms_gender_high <- subset_samples(dorms_gender_RA, `new_sheetwash`=="high")

# Filter one more time with female and male separately
dorms_high_female <- subset_samples(dorms_gender_high, `sex`=="female")
dorms_high_male <- subset_samples(dorms_gender_high, `sex`=="male")

metadata_of_female_high <- sample_data(dorms_high_female)
# only 8 samples
# but different when it comes to ASVs?
metadata_of_male_high <- sample_data(dorms_high_male)
# 20 samples

# Venn Diagram
# Setting detection and prevalence thresholds
female_high_ASVs <- core_members(dorms_high_female, detection=0.001, prevalence = 0.75)
male_high_ASVs <- core_members(dorms_high_male, detection=0.001, prevalence = 0.75)

dorms_list_high_gender <- list(Female_High_Frequency = female_high_ASVs, Male_High_Frequency = male_high_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to 2 groups
dorms_high_gender <- ggVennDiagram(x = dorms_list_high_gender)

dorms_high_gender

ggsave("venn_dorms_high_bothgender.png", dorms_high_gender)


#What are these ASVs?
high_female_tax <- tax_table(prune_taxa(female_high_ASVs,dorms_final))
high_male_tax <- tax_table(prune_taxa(male_high_ASVs,dorms_final))

#Unique Genera/species for low vs high
tax_high_diff <- setdiff(high_female_tax, high_male_tax)
print(tax_high_diff)

tax_high_diff2 <- setdiff(high_male_tax, high_female_tax)
print(tax_high_diff2)