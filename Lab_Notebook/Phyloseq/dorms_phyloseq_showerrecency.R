library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#Load unfiltered phyloseq object
load("Phyloseq/dorms_unfiltered.RData")

# Remove non-bacterial sequences, if any
dorms_filt <- subset_taxa(dorms,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove samples with less than 100 reads
dorms_filt_samps <- prune_samples(sample_sums(dorms_filt)>100, dorms_filt)
# Remove samples where shower recency is na
dorms_final <- subset_samples(dorms_filt_samps, !is.na(last_shower_binned))

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(dorms_final))), cex=0.1)
dorms_rare <- rarefy_even_depth(dorms_final, rngseed = 1, sample.size = 6223)

## Saving
save(dorms_final, file="dorms_final_showerrecency.RData")
save(dorms_rare, file="dorms_rare_showerrecency.RData")
