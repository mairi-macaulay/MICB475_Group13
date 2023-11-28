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

####Creating bar graphs with color!##########
##Create new phyloseq for each condition##
#female high vs male high#
p_femalehigh_vs_malehigh <- subset_samples(dorms_final, sex_sheetwashfreq %in% c("female high", "male high"))
#female low vs male low#
p_femalelow_vs_malelow <- subset_samples(dorms_final, sex_sheetwashfreq %in% c("female low", "male low"))

##DESeq Object Creation##
#adding +1 to all counts in the OTU table to correct for zero's that DESeq cant handle
phyloseq_object_plus1_gender_high <- transform_sample_counts(p_femalehigh_vs_malehigh, function(x) x+1)
phyloseq_object_plus1_gender_low <- transform_sample_counts(p_femalelow_vs_malelow, function(x) x+1)
#turning phloseq object to deseq object
sheetwash_deseq_gender_high <- phyloseq_to_deseq2(phyloseq_object_plus1_gender_high, ~`sex_sheetwashfreq`)
sheetwash_deseq_gender_low <- phyloseq_to_deseq2(phyloseq_object_plus1_gender_low, ~`sex_sheetwashfreq`)
#running DESeq
DESEQ_sheetwash_gender_high <- DESeq(sheetwash_deseq_gender_high)
DESEQ_sheetwash_gender_low <- DESeq(sheetwash_deseq_gender_low)


#### Female high vs Male High ###
res <- results(DESEQ_sheetwash_gender_high, tidy=TRUE, contrast= c("sex_sheetwashfreq","female high","male high"))
#View(res)

### Creating the Volcano plot: effect size VS significance ###
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot =  res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

### Getting a table of Results ###
sigASVs_gender_high <- as.data.frame(res) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  dplyr::rename(ASV=row)
#View(sigASVs)
#Significant ASVs
sigASVs_vec_gender_high <- sigASVs_gender_high %>%
  pull(ASV)

#Prune phyloseq file
sheetwash_DESeq_pruned_gender_high <- prune_taxa(sigASVs_vec_gender_high,dorms_final)

sheetwash_sigASVs_gender_high <- tax_table(sheetwash_DESeq_pruned_gender_high) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_gender_high) %>%
  arrange(log2FoldChange)

sheetwash_sigASVs_gender_high <- sheetwash_sigASVs_gender_high[,-1]
sheetwash_sigASVs_gender_high$Phylum <- gsub("^...","",sheetwash_sigASVs_gender_high$Phylum)
sheetwash_sigASVs_gender_high$Class <- gsub("^...","",sheetwash_sigASVs_gender_high$Class)
sheetwash_sigASVs_gender_high$Order <- gsub("^...","",sheetwash_sigASVs_gender_high$Order)
sheetwash_sigASVs_gender_high$Family <- gsub("^...","",sheetwash_sigASVs_gender_high$Family)
sheetwash_sigASVs_gender_high$Genus <- gsub("^...","",sheetwash_sigASVs_gender_high$Genus)
sheetwash_sigASVs_gender_high$Species <- gsub("^...","",sheetwash_sigASVs_gender_high$Species)

high_df = sheetwash_sigASVs_gender_high

high_list = high_df$Genus

## Female Low vs. Male Low ##
res_gender_low <- results(DESEQ_sheetwash_gender_low, tidy=TRUE, contrast= c("sex_sheetwashfreq","female low","male low"))

### Creating the Volcano plot: effect size VS significance ###
ggplot(res_gender_low) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

volcano_plot_gender_low = res_gender_low %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

### Getting a table of Results ###
sigASVs_gender_low <- as.data.frame(res_gender_low) %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2 & baseMean > 1) %>%
  dplyr::rename(ASV=row)
#Significant ASVs
sigASVs_vecs_gender_low <- sigASVs_gender_low %>%
  pull(ASV)
#Prune phyloseq file
sheetwash_DESeq_pruned_gender_low <- prune_taxa(sigASVs_vecs_gender_low,dorms_final)

sheetwash_sigASVs_gender_low <- tax_table(sheetwash_DESeq_pruned_gender_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_gender_low) %>%
  arrange(log2FoldChange)

sheetwash_sigASVs_gender_low <- sheetwash_sigASVs_gender_low[,-1]
sheetwash_sigASVs_gender_low$Phylum <- gsub("^...","",sheetwash_sigASVs_gender_low$Phylum)
sheetwash_sigASVs_gender_low$Class <- gsub("^...","",sheetwash_sigASVs_gender_low$Class)
sheetwash_sigASVs_gender_low$Order <- gsub("^...","",sheetwash_sigASVs_gender_low$Order)
sheetwash_sigASVs_gender_low$Family <- gsub("^...","",sheetwash_sigASVs_gender_low$Family)
sheetwash_sigASVs_gender_low$Genus <- gsub("^...","",sheetwash_sigASVs_gender_low$Genus)
sheetwash_sigASVs_gender_low$Species <- gsub("^...","",sheetwash_sigASVs_gender_low$Species)

low_df = sheetwash_sigASVs_gender_low

#list of genus from low
low_list = low_df$Genus

#Gathering all unique genus names from both lists
all_genus_list = na.omit(unique(append(low_list,high_list)))

#Now for the shared/unique part
shared_genus = matrix()
low_genus = matrix()
high_genus = matrix()
for(genus in all_genus_list){
  #genus = "Kocuria"
  if(genus %in% low_list & genus %in% high_list){
    print(paste(genus,"is shared!"))
    shared_genus = append(shared_genus, genus)
    
  }else if(genus %in% low_list){
    print(paste(genus,"is in Low list!"))
    low_genus = append(low_genus, genus)
    
  }else{
    print(paste(genus,"is in High list!"))
    high_genus = append(high_genus, genus)
    
  }
  
}

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
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1))+ ylab('Genus') +
  scale_fill_manual(values = c("red","purple"))


###High group
high_df <- high_df %>%
  mutate(Shared = ifelse((grepl(paste(shared_genus, collapse = "|"), Genus)), "Shared",
                         ifelse((grepl(paste(high_genus, collapse = "|"), Genus)), "Unique","")))
high_df_merged = high_df %>%
  group_by(Genus, Shared) %>%
  summarize(log2FoldChange_avg = mean(log2FoldChange),lfcSE_avg = mean(lfcSE))

high_df_merged <- high_df_merged[order(high_df_merged$log2FoldChange_avg),]
femalehigh_vs_malehigh_barplot <- ggplot(high_df_merged) +
  geom_bar(aes(y= reorder(Genus, sort(as.numeric(log2FoldChange_avg))), x=log2FoldChange_avg, fill = Shared), stat="identity" )+
  geom_errorbar(aes(y=Genus, xmin=log2FoldChange_avg-lfcSE_avg, xmax=log2FoldChange_avg+lfcSE_avg)) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1))+ ylab('Genus') +
  scale_fill_manual(values = c("orange","blue"))





