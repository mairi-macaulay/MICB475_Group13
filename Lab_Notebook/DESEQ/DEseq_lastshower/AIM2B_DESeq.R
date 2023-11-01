#AIM2B_DEseq Analysis

# Load packages
library(phyloseq) 
library(ape)
library(tidyverse)
library(vegan)
library(FSA)
library(DESeq2)

# Load data 
load("dorms_final_showerrecency_deseq.RData")

#Filter Males? (Females?)

# DESeq for "recent" and "not_recent" last shower day
dorms_last_shower_deseq <- phyloseq_to_deseq2(dorms_final, ~`last_shower_binned`)
DESEQ_dorms_last_shower <- DESeq(dorms_last_shower_deseq)

res_last_shower <- results(DESEQ_dorms_last_shower, tidy=TRUE)
View(res_last_shower)

# Look at results 

# Volcano plot: effect size VS significance
# Make variable to color by whether it is significant + large change
vol_plot_last_shower <- res_last_shower %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="shower_recency_vol_plot.png",vol_plot_last_shower)

# To get table of results
sigASVs_last_shower <- res_last_shower %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_last_shower)

# Get only significant ASV names
sigASVs_last_shower_vec <- sigASVs_last_shower %>%
  pull(ASV)
#There are 1 ASV significantly different between two groups

# Prune phyloseq file
dorm_last_shower_DESeq <- prune_taxa(sigASVs_last_shower_vec,dorms_final)

# Phlyum level comparison
dorms_last_shower_sigASVs <- tax_table(dorm_last_shower_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_last_shower) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

ggplot(dorms_last_shower_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Genus level comparison
dorms_last_shower_sigASVs <- tax_table(dorm_last_shower_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_last_shower) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(dorms_last_shower_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))

# Species level comparison
dorms_last_shower_sigASVs <- tax_table(dorm_last_shower_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_last_shower) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

ggplot(dorms_last_shower_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))
#Species: NA
