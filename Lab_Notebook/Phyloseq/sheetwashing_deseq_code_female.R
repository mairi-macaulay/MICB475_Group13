
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
load("Lab_Notebook/Phyloseq/AIM_2B_phyloseq/dorms_final_sheetwashfreq_deseq_female.RData")


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


### Creating the Volcano plot: effect size VS significance ###
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

### Getting a table of Results ###
sigASVs <- as.data.frame(res) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  dplyr::rename(ASV=row)

#Significant ASVs
sigASVs_vec <- sigASVs %>%
  pull(ASV)



### Creating Bar plots ###
#Prune phyloseq file
sheetwash_DESeq_pruned <- prune_taxa(sigASVs_vec,dorms_final_sheetwashfreq_deseq_female)

sheetwash_sigASVs_female <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange)

sheetwash_sigASVs_female <- sheetwash_sigASVs_female[,-1]
sheetwash_sigASVs_female$Phylum <- gsub("^...","",sheetwash_sigASVs_female$Phylum)
sheetwash_sigASVs_female$Class <- gsub("^...","",sheetwash_sigASVs_female$Class)
sheetwash_sigASVs_female$Order <- gsub("^...","",sheetwash_sigASVs_female$Order)
sheetwash_sigASVs_female$Family <- gsub("^...","",sheetwash_sigASVs_female$Family)
sheetwash_sigASVs_female$Genus <- gsub("^...","",sheetwash_sigASVs_female$Genus)
sheetwash_sigASVs_female$Species <- gsub("^...","",sheetwash_sigASVs_female$Species)

female_df = sheetwash_sigASVs_female

female_list = female_df$Genus

#Gathering unique genus names from both lists
all_genus_list = na.omit(unique(female_list))


###Low group
low_df <- low_df %>%
  mutate(Shared = ifelse((grepl(paste(shared_genus, collapse = "|"), Genus)), "Shared",
                         ifelse((grepl(paste(low_genus, collapse = "|"), Genus)), "Unique","")))
low_df = low_df[-c(5,21),]

low_df_merged = low_df %>%
  group_by(Genus, Shared) %>%
  summarize(log2FoldChange_avg = mean(log2FoldChange),lfcSE_avg = mean(lfcSE))

low_df_merged <- low_df_merged[order(low_df_merged$log2FoldChange_avg),]
femalelow_vs_malelow_barplot <- ggplot(low_df_merged) +
  geom_bar(aes(y=reorder(Genus, sort(as.numeric(log2FoldChange_avg))), x=log2FoldChange_avg, fill = Shared), stat="identity") +
  geom_errorbar(aes(y=Genus, xmin=log2FoldChange_avg-lfcSE_avg, xmax=log2FoldChange_avg+lfcSE_avg)) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))+ ylab('Genus') +
  scale_fill_manual(values = c("orange","purple"))








# # Phlyum level comparison
# phylum_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Phylum = make.unique(Phylum)) %>%
#   mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
# 
# barplot_phyla_high_low = ggplot(phylum_sheetwash_sigASVs) +
#   geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# 
# 
# # Genus level comparison
# genus_sheetwash_sigASVs <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Genus = make.unique(Genus)) %>%
#   mutate(Genus = factor(Genus, levels=unique(Genus)))
# 
# barplot_genus_high_low = ggplot(genus_sheetwash_sigASVs) +
#   geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 


# # Species level comparison
# species_sheetwash_sigASVs  <- tax_table(sheetwash_DESeq_pruned) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Species = make.unique(Species)) %>%
#   mutate(Species = factor(Species, levels=unique(Species)))
# 
# barplot_species_high_low = ggplot(species_sheetwash_sigASVs) +
#   geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 













# ###viewing DESeq results- comparison group 2
# #high group is the comparison group and low group is reference
# res_med_low <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sheetwashfreq_binned","medium","low"))
# 
# 
# ### Creating the Volcano plot: effect size VS significance ###
# ggplot(res_med_low) +
#   geom_point(aes(x=log2FoldChange, y=-log10(padj)))
# 
# volcano_plot_med_low =  res_med_low %>%
#   mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
#   ggplot() +
#   geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# 
# #saving file
# #ggsave(filename="volcano_plot_med_low_female.png",volcano_plot_med_low)
# 
# ### Getting a table of Results ###
# sigASVs_med_low <- as.data.frame(res_med_low) %>% 
#   filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
#   dplyr::rename(ASV=row)
# #Significant ASVs
# sigASVs_vec_med_low <- sigASVs_med_low %>%
#   pull(ASV)
# 
# 
# ### Creating Bar plots ###
# #Prune phyloseq file
# sheetwash_DESeq_pruned_med_low <- prune_taxa(sigASVs_vec_med_low,dorms_final_sheetwashfreq_deseq_female)
# 
# # Phlyum level comparison
# phylum_sheetwash_sigASVs_med_low <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_med_low) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Phylum = make.unique(Phylum)) %>%
#   mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
# 
# barplot_phyla_med_low = ggplot(phylum_sheetwash_sigASVs_med_low) +
#   geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_phyla_med_low_female.png",barplot_phyla_med_low)
# 
# 
# # Genus level comparison
# genus_sheetwash_sigASVs_med_low <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_med_low) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Genus = make.unique(Genus)) %>%
#   mutate(Genus = factor(Genus, levels=unique(Genus)))
# 
# barplot_genus_med_low = ggplot(genus_sheetwash_sigASVs_med_low) +
#   geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_genus_med_low_female.png",barplot_genus_med_low)
# 
# 
# # Species level comparison
# species_sheetwash_sigASVs_med_low  <- tax_table(sheetwash_DESeq_pruned_med_low) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_med_low) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Species = make.unique(Species)) %>%
#   mutate(Species = factor(Species, levels=unique(Species)))
# 
# barplot_species_med_low = ggplot(species_sheetwash_sigASVs_med_low) +
#   geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + 
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_species_med_low_female.png", barplot_species_med_low)













# 
# ###viewing DESeq results- comparison group 3
# #high group is the comparison group and low group is reference
# res_high_med <- results(DESEQ_sheetwash, tidy=TRUE, contrast= c("sheetwashfreq_binned","high","medium"))
# 
# ### Creating the Volcano plot: effect size VS significance ###
# ggplot(res_high_med) +
#   geom_point(aes(x=log2FoldChange, y=-log10(padj)))
# 
# volcano_plot_high_med =  res_high_med %>%
#   mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
#   ggplot() +
#   geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
# 
# #saving file
# #ggsave(filename="volcano_plot_high_med_female.png",volcano_plot_high_med)
# 
# ### Getting a table of Results ###
# sigASVs_high_med <- as.data.frame(res_high_med) %>% 
#   filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
#   dplyr::rename(ASV=row)
# #Significant ASVs
# sigASVs_vec_high_med <- sigASVs_high_med %>%
#   pull(ASV)
# 
# 
# 
# ### Creating Bar plots ###
# #Prune phyloseq file
# sheetwash_DESeq_pruned_high_med <- prune_taxa(sigASVs_vec_high_med,dorms_final_sheetwashfreq_deseq_female)
# 
# # Phlyum level comparison
# phylum_sheetwash_sigASVs_high_med <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_high_med) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Phylum = make.unique(Phylum)) %>%
#   mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
# 
# barplot_phyla_high_med = ggplot(phylum_sheetwash_sigASVs_high_med) +
#   geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity") +
#   geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_phyla_high_med_female.png",barplot_phyla_high_med)
# 
# 
# # Genus level comparison
# genus_sheetwash_sigASVs_high_med <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_high_med) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Genus = make.unique(Genus)) %>%
#   mutate(Genus = factor(Genus, levels=unique(Genus)))
# 
# barplot_genus_high_med = ggplot(genus_sheetwash_sigASVs_high_med) +
#   geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_genus_high_med_female.png",barplot_genus_high_med)
# 
# 
# # Species level comparison
# species_sheetwash_sigASVs_high_med  <- tax_table(sheetwash_DESeq_pruned_high_med) %>% as.data.frame() %>%
#   rownames_to_column(var="ASV") %>%
#   right_join(sigASVs_high_med) %>%
#   arrange(log2FoldChange) %>%
#   mutate(Species = make.unique(Species)) %>%
#   mutate(Species = factor(Species, levels=unique(Species)))
# 
# barplot_species_high_med = ggplot(species_sheetwash_sigASVs_high_med) +
#   geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
#   geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + 
#   theme(text = element_text(size=8),
#         axis.text.x = element_text(angle=90, hjust=1)) 
# 
# #ggsave(filename="barplot_species_high_med_female.png", barplot_species_high_med)
