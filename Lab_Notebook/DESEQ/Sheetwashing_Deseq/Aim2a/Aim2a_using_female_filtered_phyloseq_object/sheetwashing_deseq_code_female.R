# if you didn't install the DESeq2 package, run the following
BiocManager::install("DESeq2")

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ape)
library(vegan)
library(FSA)

#setting the seed
set.seed(1)

#### Loading data ####
#Load dorms_final (filtered data)
load("Lab_Notebook/Phyloseq/AIM_2A_phyloseq/dorms_final_sheetwashfreq_deseq_female.RData")


#### DESeq Object Creation ####
#adding +1 to all counts in the OTU table to correct for zero's that DESeq cant handle
phyloseq_object_plus1 <- transform_sample_counts(dorms_final_sheetwashfreq_deseq_female, function(x) x+1)
#turning phloseq object to deseq object
sheetwash_deseq <- phyloseq_to_deseq2(phyloseq_object_plus1, ~`sheetwashfreq_binned`)
#running DESeq
DESEQ_sheetwash <- DESeq(sheetwash_deseq)





###viewing DESeq results- comparison group 1 
#high group is the comparison group and low group is reference
res <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sheetwashfreq_binned","high","low"))
View(res)

### Creating the Volcano plot: effect size VS significance ###
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#saving file
ggsave(path = Lab_Notebook/DESEQ/Sheetwashing_deseq/Aim2a/Aim2a- using dorms_final_sheetwashfreq_deseq_female phyloseq object, filename="volcano_plot_high_low_female.png",volcano_plot)

### Getting a table of Results ###
sigASVs <- as.data.frame(res) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
#Significant ASVs
sigASVs_vec <- sigASVs %>%
  pull(ASV)
#There are 45 significant ASV's
view(sigASVs_vec)


### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned <- prune_taxa(sigASVs_vec,dorms_final_sheetwashfreq_deseq_female)

# Phlyum level comparison
phylum_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

barplot_phyla_high_low = ggplot(phylum_sheetwash_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_phyla_high_low_female.png",barplot_phyla_high_low)


# Genus level comparison
genus_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_genus_high_low = ggplot(genus_sheetwash_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_genus_high_low.png",barplot_genus_high_low)


# Species level comparison
species_sheetwash_sigASVs  <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

barplot_species_high_low = ggplot(species_sheetwash_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_species_high_low.png", barplot_species_high_low)












###viewing DESeq results- comparison group 2
#high group is the comparison group and low group is reference
res_med_low <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sheetwashfreq_binned","medium","low"))
View(res_med_low)

### Creating the Volcano plot: effect size VS significance ###
ggplot(res_med_low) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot_med_low =  res_med_low %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#saving file
ggsave(filename="volcano_plot_med_low.png",volcano_plot_med_low)

### Getting a table of Results ###
sigASVs_med_low <- as.data.frame(res) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_med_low)
#Significant ASVs
sigASVs_vec_med_low <- sigASVs_med_low %>%
  pull(ASV)
#There are 45 significant ASV's
view(sigASVs_vec_med_low)


### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned_med_low <- prune_taxa(sigASVs_vec_med_low,dorms_final)

# Phlyum level comparison
phylum_sheetwash_sigASVs_med_low <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

barplot_phyla_med_low = ggplot(phylum_sheetwash_sigASVs_med_low) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_phyla_med_low.png",barplot_phyla_med_low)


# Genus level comparison
genus_sheetwash_sigASVs_med_low <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_genus_med_low = ggplot(genus_sheetwash_sigASVs_med_low) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_genus_med_low.png",barplot_genus_med_low)


# Species level comparison
species_sheetwash_sigASVs_med_low  <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_med_low) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

barplot_species_med_low = ggplot(species_sheetwash_sigASVs_med_low) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + 
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_species_med_low.png", barplot_species_med_low)














###viewing DESeq results- comparison group 3
#high group is the comparison group and low group is reference
res_high_med <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sheetwashfreq_binned","high","medium"))
View(res_high_med)

### Creating the Volcano plot: effect size VS significance ###
ggplot(res_high_med) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot_high_med =  res_high_med %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#saving file
ggsave(filename="volcano_plot_high_med.png",volcano_plot_high_med)

### Getting a table of Results ###
sigASVs_high_med <- as.data.frame(res) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_high_med)
#Significant ASVs
sigASVs_vec_high_med <- sigASVs_high_med %>%
  pull(ASV)
#There are 45 significant ASV's
view(sigASVs_vec_high_med)


### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned_high_med <- prune_taxa(sigASVs_vec_high_med,dorms_final)

# Phlyum level comparison
phylum_sheetwash_sigASVs_high_med <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

barplot_phyla_high_med = ggplot(phylum_sheetwash_sigASVs_high_med) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity") +
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_phyla_high_med.png",barplot_phyla_high_med)


# Genus level comparison
genus_sheetwash_sigASVs_high_med <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_genus_high_med = ggplot(genus_sheetwash_sigASVs_high_med) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_genus_high_med.png",barplot_genus_high_med)


# Species level comparison
species_sheetwash_sigASVs_high_med  <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_high_med) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

barplot_species_high_med = ggplot(species_sheetwash_sigASVs_high_med) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + 
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggsave(filename="barplot_species_high_med.png", barplot_species_high_med)
