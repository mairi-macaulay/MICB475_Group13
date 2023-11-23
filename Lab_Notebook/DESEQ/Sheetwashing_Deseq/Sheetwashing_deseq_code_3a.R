### Gender Comparison of "high" and "low" sheetwashing
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

##setting a seed##
set.seed(1) 

#Load dorms_final (filtered data)
load("Lab_Notebook/DESEQ/Sheetwashing_deseq/dorms_final_sheetwashfreq_deseq.RData")

##Adding combined sex and sheet washing frequency column##
sample_data(dorms_final)$sex_sheetwashfreq <- paste(sample_data(dorms_final)$sex, sample_data(dorms_final)$sheetwashfreq_binned)

##Create new phyloseq for each condition##
#female high vs male high#
p_malehigh_vs_femalehigh <- subset_samples(dorms_final, sex_sheetwashfreq %in% c("male high", "female high"))
#female low vs male low#
p_malelow_vs_femalelow <- subset_samples(dorms_final, sex_sheetwashfreq %in% c("male low", "female low"))

View(p_malehigh_vs_femalehigh)

##DESeq Object Creation##
#adding +1 to all counts in the OTU table to correct for zero's that DESeq cant handle
phyloseq_object_plus1_gender_high <- transform_sample_counts(p_malehigh_vs_femalehigh, function(x) x+1)
phyloseq_object_plus1_gender_low <- transform_sample_counts(p_malelow_vs_femalelow, function(x) x+1)
#turning phloseq object to deseq object
sheetwash_deseq_gender_high <- phyloseq_to_deseq2(phyloseq_object_plus1_gender_high, ~`sex_sheetwashfreq`)
sheetwash_deseq_gender_low <- phyloseq_to_deseq2(phyloseq_object_plus1_gender_low, ~`sex_sheetwashfreq`)
#running DESeq
DESEQ_sheetwash_gender_high <- DESeq(sheetwash_deseq_gender_high)
DESEQ_sheetwash_gender_low <- DESeq(sheetwash_deseq_gender_low)



#### Male high vs Female High ###
###viewing DESeq results- comparison group 1 
#high group is the comparison group and low group is reference
res <- results(DESEQ_sheetwash_gender_high, tidy=TRUE, contrast= c("sex_sheetwashfreq","male high","female high"))
View(res)

### Creating the Volcano plot: effect size VS significance ###
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#saving file
#ggsave(filename="volcano_plot_high_low.png",volcano_plot)

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
sheetwash_DESeq_pruned <- prune_taxa(sigASVs_vec,dorms_final)

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

#ggsave(filename="barplot_phyla_high_low.png",barplot_phyla_high_low)


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

#ggsave(filename="barplot_genus_high_low.png",barplot_genus_high_low)


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

#ggsave(filename="barplot_species_high_low.png", barplot_species_high_low)











#### Male low vs Female low ###
###viewing DESeq results- comparison group 1 
#high group is the comparison group and low group is reference
res_gender_low <- results(DESEQ_sheetwash_gender_low, tidy=TRUE, contrast= c("sex_sheetwashfreq","male low","female low"))


### Creating the Volcano plot: effect size VS significance ###
ggplot(res_gender_low) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot_gender_low =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))


### Getting a table of Results ###
sigASVs_gender_low <- as.data.frame(res_gender_low) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
#Significant ASVs
sigASVs_vecs_gender_low <- sigASVs_gender_low %>%
  pull(ASV)
#There are 45 significant ASV's


### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned_gender_low <- prune_taxa(sigASVs_vecs_gender_low,dorms_final)

# Phlyum level comparison
phylum_sheetwash_sigASVs_gender_low <- tax_table(sheetwash_DESeq_pruned_gender_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))

barplot_phyla_gender_low = ggplot(phylum_sheetwash_sigASVs_gender_low) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 


# Genus level comparison
genus_sheetwash_sigASVs_gender_low <- tax_table(sheetwash_DESeq_pruned_gender_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

barplot_genus_gender_low = ggplot(genus_sheetwash_sigASVs_gender_low) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 


# Species level comparison
species_sheetwash_sigASVs_gender_low  <- tax_table(sheetwash_DESeq_pruned_gender_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

barplot_species_gender_low = ggplot(species_sheetwash_sigASVs_gender_low) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) 