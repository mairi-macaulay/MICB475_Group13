
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
load("Lab_Notebook/DESEQ/Sheetwashing_deseq/dorms_final_sheetwashfreq_deseq.RData")

#### DESeq Object Creation ####
#adding +1 to all counts in the OTU table to correct for zero's that DESeq cant handle
phyloseq_object_plus1 <- transform_sample_counts(dorms_final, function(x) x+1)
#turning phloseq object to deseq object
sheetwash_deseq <- phyloseq_to_deseq2(phyloseq_object_plus1, ~`sex`)
#running DESeq
DESEQ_sheetwash <- DESeq(sheetwash_deseq)
#viewing DESeq results
#high group is the comparison group and low group is reference
res <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sex","female","male"))


### Creating the Volcano plot: effect size VS significance ###
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#saving file
#ggsave(filename="volcano_plot.png",volcano_plot)

### Getting a table of Results ###
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
#Significant ASVs
sigASVs_vec <- sigASVs %>%
  pull(ASV)


### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned <- prune_taxa(sigASVs_vec,dorms_final)

# Phlyum level comparison
phylum_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

barplot_phylum_female_male = ggplot(phylum_sheetwash_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1))


# Genus level comparison
genus_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_genus_female_male = ggplot(genus_sheetwash_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1))


# Species level comparison
species_sheetwash_sigASVs  <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

barplot_species_female_male = ggplot(species_sheetwash_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1))
