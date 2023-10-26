library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#Load unfiltered phyloseq object
load("dorms_unfiltered.RData")

# Remove non-bacterial sequences, if any
dorms_filt <- subset_taxa(dorms,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total (DO FOR DESEQ)
dorms_filt_nolow <- filter_taxa(dorms_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
dorms_filt_nolow_samps <- prune_samples(sample_sums(dorms_filt_nolow)>100, dorms_filt_nolow)
# Remove samples where sheet wash frequency is na
dorms_final <- subset_samples(dorms_filt_nolow_samps, !is.na(sheetwashfreq_binned))


# Rarefy samples
rarecurve(t(as.data.frame(otu_table(dorms_final))), cex=0.1)
dorms_rare <- rarefy_even_depth(dorms_final, rngseed = 1, sample.size = 6223)

## Saving
save(dorms_final, file="dorms_final_sheetwashfreq_deseq.RData")
save(dorms_rare, file="dorms_rare_sheetwashfreq_deseq.RData")
