# 4 way venn diagram
#(Sheet Washing Frequency and Sex, NO MEDIUM GROUP)

library(ggVennDiagram)

# Load data 
load("dorms_final_sheetwashfreq.RData")

set.seed(1)

#Convert to relative abundance
dorms_RA <- transform_sample_counts(dorms_final, fun=function(x) x/sum(x))

# Filter dataset with just low group
dorms_low <- subset_samples(dorms_RA, `sheetwashfreq_binned`=="low")
dorms_high <- subset_samples(dorms_RA, `sheetwashfreq_binned`=="high")

# Filter one more time with female and male separately for all four
dorms_low_female <- subset_samples(dorms_low, `sex`=="female")
dorms_high_female <- subset_samples(dorms_high, `sex`=="female")
dorms_low_male <- subset_samples(dorms_low, `sex`=="male")
dorms_high_male <- subset_samples(dorms_high, `sex`=="male")

metadata_of_female_low <- sample_data(dorms_low_female)
metadata_of_female_high <- sample_data(dorms_high_female)

metadata_of_male_low <- sample_data(dorms_low_male)
metadata_of_male_high <- sample_data(dorms_high_male)

# Setting detection and prevalence thresholds
low_female_ASVs <- core_members(dorms_low_female, detection=0.001, prevalence = 0.5)
high_female_ASVs <- core_members(dorms_high_female, detection=0.001, prevalence = 0.5)

# Setting detection and prevalence thresholds
low_male_ASVs <- core_members(dorms_low_male, detection=0.001, prevalence = 0.5)
high_male_ASVs <- core_members(dorms_high_male, detection=0.001, prevalence = 0.5)


dorms_list <- list(Female_Low_Frequency = low_female_ASVs, 
          Female_High_Frequency = high_female_ASVs, 
          Male_Low_Frequency = low_male_ASVs, 
          Male_High_Frequency = high_male_ASVs)

dormsvenn <- ggVennDiagram(x = dorms_list, set_size = 4)
newdormsvenn <- dormsvenn + scale_x_continuous(expand = expansion(mult = .2))

finalvenn <- newdormsvenn + scale_fill_distiller(palette = "RdBu")
finalvenn
ggsave("venn_dorms_4venn.png", finalvenn)

