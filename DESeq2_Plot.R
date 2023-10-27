#!/usr/bin/env Rscript
library(phyloseq) 
library(ape)
library(tidyverse)
library(vegan)
library(FSA)
library(DESeq2)
#### Loading data #### 
load("C:/Users/cwjle/OneDrive - UBC/Desktop/GradWork/Classes/TAMICB475_2022-2023/team repos/JAAK_STAT_Pathway/dysautonomia_final.RData")

#### DESeq for Mild vs Severe FD ####
deseq_dysautonomia_mil_severe <- phyloseq_to_deseq2(dysautonomia_final,~FD.severity)
DESEQ_mild_severe <- DESeq(deseq_dysautonomia_mil_severe)
res_4 <- results(DESEQ_mild_severe, tidy=TRUE)
View(res_4)

# Volcano plot: effect size VS significance
ggplot(res_4) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
gg_volcano <- res_4 %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# Table of results
sigASVs_mild_severe <- res_4 %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_mild_severe)
# Significant ASVs
sigASVs_msevere_vec <- sigASVs_mild_severe %>%
  pull(ASV)
# There are 5 ASVs significantly different between the mild and severe groups: 721bbde09abf2f51fc7d8caeab5b22f1, 0d29a1eb24f264423251f3b2a92f7914, 273fa0191072af3d33e32271b35c8f18, acf76b1f22c9536ca982df6f9b0219da, cdb79eecfa38ac58e99bfdef45fabfce    

# Prune phyloseq file
msevere_FD_DESeq <- prune_taxa(sigASVs_msevere_vec,dysautonomia_final)
# Phlyum level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Genus level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))

# Species level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))
