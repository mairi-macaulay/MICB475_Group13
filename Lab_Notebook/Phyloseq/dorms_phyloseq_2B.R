library(phyloseq)

#Load phyloseq objects
load("dorms_rare_sheetwashfreq.RData")
load("dorms_final_sheetwashfreq_deseq.RData")

####Create phyloseq for 2B alpha diversity and taxonomic bar plots- Females#####
# Filter out males
dorms_rare_sheetwashfreq_female <- subset_samples(dorms_rare, sex != "male")

# Saving
save(dorms_rare_sheetwashfreq_female, file="dorms_rare_sheetwashfreq_female.RData")

####Create phyloseq for 2B alpha diversity and taxonomic bar plots - Males#####
# Filter out females
dorms_rare_sheetwashfreq_male <- subset_samples(dorms_rare, sex != "female")

# Saving
save(dorms_rare_sheetwashfreq_male, file="dorms_rare_sheetwashfreq_male.RData")

####Create phyloseq for 2B DESeq - Females#####
# Filter out males
dorms_final_sheetwashfreq_deseq_female <- subset_samples(dorms_final, sex != "male")

# Saving
save(dorms_final_sheetwashfreq_deseq_female, file="dorms_final_sheetwashfreq_deseq_female.RData")

####Create phyloseq for 2B DESeq - Males#####
# Filter out females
dorms_final_sheetwashfreq_deseq_male <- subset_samples(dorms_final, sex != "female")

# Saving
save(dorms_final_sheetwashfreq_deseq_male, file="dorms_final_sheetwashfreq_deseq_male.RData")
