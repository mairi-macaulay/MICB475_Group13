# Core Microbiome Aim 1B (Sheet Washing Frequency and Sex)
# Comparing male low sheetwashing freq vs female low sheetwashing freq

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

# Filter dataset with only low frequency
dorms_gender_low <- subset_samples(dorms_gender_RA, `new_sheetwash`=="low")

# Filter one more time with female and male separately
dorms_low_female <- subset_samples(dorms_gender_low, `sex`=="female")
dorms_low_male <- subset_samples(dorms_gender_low, `sex`=="male")

metadata_of_female_low <- sample_data(dorms_low_female)
metadata_of_male_low <- sample_data(dorms_low_male)


# Venn Diagram
# Setting detection and prevalence thresholds
female_low_ASVs <- core_members(dorms_low_female, detection=0.001, prevalence = 0.5)
male_low_ASVs <- core_members(dorms_low_male, detection=0.001, prevalence = 0.5)

dorms_list_low_gender <- list(Female_Low_Frequency = female_low_ASVs, Male_Low_Frequency = male_low_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to 2 groups
dorms_low_gender <- ggVennDiagram(x = dorms_list_low_gender)

dorms_low_gender

ggsave("venn_dorms_low_bothgender.png", dorms_low_gender)


#What are these ASVs?
low_female_tax <- tax_table(prune_taxa(female_low_ASVs,dorms_final))
low_male_tax <- tax_table(prune_taxa(male_low_ASVs,dorms_final))

#Unique Genera/species for low vs high
tax_low_diff <- setdiff(low_female_tax, low_male_tax)
print(tax_low_diff)

tax_low_diff2 <- setdiff(low_male_tax, low_female_tax)
print(tax_low_diff2)
