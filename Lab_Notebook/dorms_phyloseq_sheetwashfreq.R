####Create Phyloseq Object####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

## Load data
metafp <- "dorms_export/dorms_metadata_updated.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "dorms_export/dorms_feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "dorms_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "dorms_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

## Format OTU table
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)

## Format sample metadata
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'#SampleID'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

## Formatting taxonomy
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

# Merge all into a phyloseq object
dorms <- phyloseq(OTU, SAMP, TAX, phylotree)

##Do we need these steps??
# Remove non-bacterial sequences, if any
dorms_filt <- subset_taxa(dorms,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
dorms_filt_nolow <- filter_taxa(dorms_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
dorms_filt_nolow_samps <- prune_samples(sample_sums(dorms_filt_nolow)>100, dorms_filt_nolow)
# Remove samples where month is na
dorms_final <- subset_samples(dorms_filt_nolow_samps, !is.na(sheetwashfreq_binned))


# Rarefy samples
rarecurve(t(as.data.frame(otu_table(dorms_final))), cex=0.1)
dorms_rare <- rarefy_even_depth(dorms_final, rngseed = 1, sample.size = 6223)

## Saving
save(dorms_final, file="dorms_final.RData")
save(dorms_rare, file="dorms_rare.RData")
