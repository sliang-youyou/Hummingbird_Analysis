#### DeSeq for GSEA analysis
## Anusha Shankar, github: nushiamme
## Code started April 5, 2022

## Initially followed this tutorial:
## https://www.youtube.com/watch?v=OzNzO8qwwp0
## And then for visualizations
## https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
## For GSEA: https://www.youtube.com/watch?v=KY6SS4vRchY

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('DESeq2')
# BiocManager::install('EnhancedVolcano') 
# BiocManager::install('apeglm') 
# Enhanced Volcano plot help: 
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(here)
library(viridis)
library(gridExtra)
library(apeglm)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(dendextend) # For making gene trees for heatmaps
library(ComplexHeatmap)
library(cowplot) # for plot_grid function
library(enrichR)

## If you opened the .Rproj file from the OneDrive folder, reading in these files will work- 
## the relative paths should be the same
data <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "RNASeq_rawcounts_data.csv"))
meta <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "AnushaShankar_RNASeq_metadata.csv"))
star_alignment <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "star_alignment_plot.csv"))
star <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "star_gene_counts.csv"))
foldchange <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "TopFoldChanges.csv"))
## Made this one in this script
norm_counts_df <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "RNASeq_NormCounts.csv"))

## Just for prepping for GO analysis on https://usadellab.github.io/GeneExpressionPlots/#/data
#dat_exp <- read_table(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "RNASeq_rawcounts_data_ForGONet.txt"))
#dat_meta <- read_table(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "AnushaShankar_RNASeq_metadata_GONet.txt"))


my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

my_theme2 <- theme_classic(base_size = 20) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

## Viridis colors
my_gradient <- c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")
my_col_rainbows <- c("#f94144", "#f3722c", "#f8961e", "#f9844a",
                     "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1")
mycols <- c("orange", "navy", "springgreen4", "mediumorchid", "gold4", "plum1", "springgreen")

#### Don't need to re-do - this is a check, not used in analyses really.                     
## Checking on the star alignment
# head(star)
# star <- star %>%
#   rename(Sample = ï..Sample) %>%
#   mutate(Sample = gsub("c", "", Sample))
# 
meta2 <- meta %>%
  rename(Sample = X)

# starlong_meta <- star %>%
#   gather(key="Mapped", value="count", -Sample) %>%
#   left_join(., meta2, by="Sample")
# 
# starlong_meta %>%
#   mutate(Mapped_relevel = 
#          fct_relevel(Mapped, 
#                     "Overlapping_Genes", "Multimapping", "Unmapped", "Ambiguous_Features", "No_Feature")) %>%
#   ggplot(., aes(Tissue, count, fill=Mapped_relevel)) + 
#   geom_bar(position="fill", stat="identity") +
#   scale_fill_viridis(discrete = T) + my_theme + ylab("Percent counts")
# 
# 
# starlong_meta %>%
#   ggplot(., aes(Tissue, count, fill=Mapped)) + 
#   geom_bar(position="fill", stat="identity") +
#   scale_fill_viridis(discrete = T)
# 
# starlong_meta %>%
#   ggplot(., aes(Tissue, count, fill=Mapped)) + 
#   geom_bar(stat="identity", position="fill") +
#   scale_fill_viridis(discrete = T)

rownames(meta) <- unique(meta$X)
rownames(data) <- unique(data$X)
data <- subset(data, select = -c(X))
meta <- subset(meta, select = -c(X))

meta <- meta[match(colnames(data), rownames(meta)),]

## Check that all column names from the data set are present as rownames in the metadata file
all(colnames(data) %in% rownames(meta))


## Check that the cols are in the same order in the data as the rownames in the metadata
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = Tissue ~ Metabolic_State)

## Pre-filtering
## Taking out rows that have counts less than 10 reads total - recommended, not required, step
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

## Set the factor level to compare against (in our case, normothermy)
dds$Metabolic_State <- relevel(dds$Metabolic_State, ref="N")

## Run DESeq
dds <- DESeq(dds)


## PCAs
vst_dds <- vst(dds,blind=TRUE, fitType='local')


## Don't use this
# tissue_pca <- plotPCA(vst_dds, 
#         intgroup = "Tissue")  
# 
# tissue_pca + geom_point(size=4) + my_theme + scale_color_manual(values = mycols, name = "Tissue")

## Use this
plotPCA_jh = function(pp1=1, pp2=2, 
                      object, intgroup="condition", 
                      ntop=1000, returnData=FALSE) {
  
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pp1], PC2=pca$x[,pp2], group=group, intgroup.df, 
                  name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pp1:pp2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", label = "name")) + 
    geom_point(size=3) + 
    xlab(paste0("PC", pp1, ": ",round(percentVar[pp1] * 100),"% variance")) +
    ylab(paste0("PC", pp2, ": ",round(percentVar[pp2] * 100),"% variance")) +
    #scale_color_brewer(type = "qual", palette = "Dark2", direction = 1)+
    colorspace::scale_color_discrete_qualitative(palette = "Dark 3", rev = TRUE)+
    cowplot::theme_cowplot()+
    coord_fixed()
}

pc12 <- plotPCA_jh(pp1 = 1, pp2 = 2, object=vst_dds, intgroup = "Tissue") + geom_point(size=4) + 
  guides(col="none") + my_theme2 + scale_color_manual(values = mycols)
pc13 <- plotPCA_jh(pp1 = 1, pp2 = 3, object=vst_dds, intgroup = "Tissue") + 
  geom_point(size=4) + my_theme2 + scale_color_manual(values = mycols, name = "Tissue")

grid.arrange(pc12, pc13, ncol=2, nrow=1, widths=c(1.5,2), heights = c(1,1))
library(cowplot)
plot_grid(pc12, pc13, align = "h", rel_widths = c(0.45, 0.55))

## Sample tissue one- just heart
heart_dds <- dds[,dds$Tissue=="Heart"]
vstdds_heart <- vst(heart_dds,blind=TRUE, fitType='local')
pc12_heart <- plotPCA_jh(pp1 = 1, pp2 = 2, object=vstdds_heart, intgroup = "Metabolic_State") + geom_point(size=4) + 
  guides(col="none") + my_theme2 + scale_color_manual(values = mycols)
pc13_heart <- plotPCA_jh(pp1 = 1, pp2 = 3, object=vstdds_heart, intgroup = "Metabolic_State") + 
  geom_point(size=4) + my_theme2 + scale_color_manual(values = mycols, name = "Metabolic_State", 
                                                      labels=c("Normothermy", "Transition", "Deep torpor"))

#grid.arrange(pc12, pc13, ncol=2, nrow=1, widths=c(1.5,2), heights = c(1,1))
plot_grid(pc12_heart, pc13_heart, align = "h", rel_widths = c(1.5, 2))


## Get normalized data, save it RE-RUN only if necessary
normalized_counts <- counts(dds, normalized=TRUE)
# norm_counts_df <- as.data.frame(normalized_counts)

#### Needed for later analyses
## Make data long-form and merge with metadata file
meta2$Metabolic_State <- as.factor(meta2$Metabolic_State)
datlong <- norm_counts_df %>%
  #rename(gene = X) %>%
  gather(key = 'Sample', value= 'counts', 3:ncol(norm_counts_df)) %>% 
  left_join(., meta2, by="Sample") %>%
  mutate(Metabolic_State = 
           fct_relevel(Metabolic_State, 
                       "N", "T", "D")) 

## Make data long-form and merge with metadata file
foldlong <- foldchange %>%
  select(Tissue, gene, D_vs_N_FoldChange, D_vs_N_log2FoldChange, D_vs_N_padj,
         T_vs_D_FoldChange, T_vs_D_log2FoldChange, T_vs_D_padj, T_vs_N_FoldChange, T_vs_N_log2FoldChange, T_vs_N_padj) %>%
  gather(key = 'Measure', value= 'value', -c(gene,Tissue))


#### Don't re-run unless you need to make the norm counts df again
## write.csv(norm_counts_df, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//RNASeq_NormCounts.csv"))

## For GSEA, make a file that has all the groups in a row in the order they are in in the columns of the the normalized 
## counts dataset. Then in Excel, add a top row that has e.g. <119 2 1> where 
## 119 is the number of samples; 2 is the number of groups and 1 is a constant for all files.
## Add a second row that has <# Trt1 Trt2 Cntrl> with spaces in between where each substring is the name of the treatment group 
## Encompassing all the cell values below
meta2$Tissue_State <- paste0(meta2$Tissue, "_", meta2$Metabolic_State)

phenotype_labs <- data.frame()

phenotype_labs_state <- rbind(phenotype_labs, meta2$Metabolic_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_state.csv"))

phenotype_labs_tissue <- rbind(phenotype_labs, meta2$Tissue)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_tissue.csv"))

phenotype_labs_state <- rbind(phenotype_labs, meta2$Tissue_State)
write.csv(phenotype_labs_state, file = here("..//DESeq_Data_Mar2022_all tissues//all tissues//final//Pheno_labs_tissue_state.csv"))


## If you have technical replicates, collapse them now. We don't.. we only have biological replicates. 
# Do not collapse those.


### can run here onwards every time
## Take a look at the data a bit
## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## Save results
res_unshrunken <- results(dds)

res <- results(dds)


## Take a quick look at the results
res_unshrunken

## Summary of results
summary(res)

## Compare different pairs
res_ND <- results(dds, contrast=c("Metabolic_State", "D", "N"))
summary(res_ND)
res_NT <- results(dds, contrast=c("Metabolic_State", "T", "N"))
summary(res_NT)
res_TD <- results(dds, contrast=c("Metabolic_State", "D", "T"))
summary(res_TD)


### Subsetting results per tissue type
res_Liver_ND <- results(dds[dds$Tissue=="Liver",], contrast=c("Metabolic_State", "D", "N"))

res_Heart_ND <- results(dds[dds$Tissue=="Heart",], contrast=c("Metabolic_State", "D", "N"))

res_Pect_ND <- results(dds[dds$Tissue=="Pect",], contrast=c("Metabolic_State", "D", "N"))

res_Lungs_ND <- results(dds[dds$Tissue=="Lungs",], contrast=c("Metabolic_State", "D", "N"))


## Trying to make a res file per tissue type
results(dds, contrast=c("Metabolic_State", "D", "N"),)

# ### Adjusted p values of 0.01
# res_ND_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "N"), alpha=0.01)
# summary(res_ND_0.01)
# res_NT_0.01 <- results(dds, contrast=c("Metabolic_State", "T", "N"), alpha=0.01)
# summary(res_NT_0.01)
# res_TD_0.01 <- results(dds, contrast=c("Metabolic_State", "D", "T"), alpha=0.01)
# summary(res_TD_0.01)

## Shrinking
res <- lfcShrink(dds, coef = 3, res = res)
#res


## Fold change results
mcols(res, use.names=T)
res %>% data.frame() %>% View()
## Summarize results
summary(res)

## We can easily subset the results table to only include 
## those that are significant using the filter() function, but first we will convert the results table into a tibble:
res_tb <- res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

## CHECK THIS for inflection point
ggplot(res_tb, aes(log2FoldChange)) + geom_density() + geom_vline(xintercept = c(-0.5, 0.5), col="red") + my_theme


### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58 ## 0.58 is equal to fold change of 1.5


## Just subset significantly different genes
sig <- res_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


# Create tibbles including row names of normalized counts
normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# ## Order results by padj values
# top30_sig_genes <- res_tb %>% 
#   arrange(padj) %>% 	#Arrange rows by padj values
#   pull(gene) %>% 		#Extract character vector of ordered genes
#   head(n=30) 		#Extract the first 30 genes

# ## Order results by padj values, just ND and just upreg genes
# top30_sig_genes_upreg <- res_ND_tb %>% 
#   arrange(padj) %>% 	#Arrange rows by padj values
#   filter(log2FoldChange>1) %>% ### NOT WORKING YET
#   pull(gene) %>% 		#Extract character vector of ordered genes
#   head(n=30) 		#Extract the first 30 genes

## For enrichr analysis, just filtering out genes that have basemean >10, logFC >1 or < -1,
# And adjusted p-val of < 0.05

## Aug 25 trying out L2FC +/- 0.5

upreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## For Transition vs. Normo
upreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_NT <-  res_NT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## For deep torpor vs. transition
upreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_TD <-  res_TD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange< -lfc.cutoff & padj < padj.cutoff & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

## normalized counts for top significantly upregulated genes
top_upreg_ND_norm <- normalized_counts %>%
  filter(gene %in% upreg_ND)

# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_ND <- top_upreg_ND_norm %>%
  gather(colnames(top_upreg_ND_norm)[2:120], key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_upreg_ND)

gathered_top_upreg_ND <- inner_join(meta2, gathered_top_upreg_ND, by="Sample")

## Plot this subset of these upreg genes
ggplot(gathered_top_upreg_ND) + #facet_grid(.~Metabolic_State) +
  geom_point(aes(x = gene, y = normalized_counts, color = Tissue)) +
  scale_y_log10() + 
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top sig upregulated genes Deep Torpor vs. Normo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d()


### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% upreg_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Upreg and downreg per tissue, just for ND comparison
upreg_Liver_ND <-  res_Liver_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_Liver_ND <-  res_Liver_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -1 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_Heart_ND <-  res_Heart_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_Heart_ND <-  res_Heart_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_Lungs_ND <-  res_Lungs_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_Lungs_ND <-  res_Lungs_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

upreg_Pect_ND <-  res_Pect_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)

downreg_Pect_ND <-  res_Pect_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
  arrange(padj) %>%
  pull(gene)


liver_samples <- meta2$Sample[meta2$Tissue=="Liver"]
heart_samples <- meta2$Sample[meta2$Tissue=="Heart"]
lungs_samples <- meta2$Sample[meta2$Tissue=="Lungs"]
pect_samples <- meta2$Sample[meta2$Tissue=="Pect"]


## Liver UP
## normalized counts for top significantly upregulated genes
top_upreg_Liver_ND_norm <- normalized_counts[,c("gene",liver_samples)] %>%
  filter(gene %in% upreg_Liver_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_Liver_ND <- top_upreg_Liver_ND_norm %>%
  gather(colnames(top_upreg_Liver_ND_norm)[2:ncol(top_upreg_Liver_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_upreg_Liver_ND)

gathered_top_upreg_Liver_ND <- inner_join(meta2, gathered_top_upreg_Liver_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_Liver <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% upreg_Liver_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Liver DOWN
## normalized counts for top significantly downregulated genes
top_downreg_Liver_ND_norm <- normalized_counts[,c("gene",liver_samples)] %>%
  filter(gene %in% downreg_Liver_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_Liver_ND <- top_downreg_Liver_ND_norm %>%
  gather(colnames(top_downreg_Liver_ND_norm)[2:ncol(top_downreg_Liver_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top_downreg_Liver_ND)

gathered_top_downreg_Liver_ND <- inner_join(meta2, gathered_top_downreg_Liver_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_Liver <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% downreg_Liver_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

## Heart UP
## normalized counts for top significantly upregulated genes
top_upreg_Heart_ND_norm <- normalized_counts[,c("gene",heart_samples)] %>%
  filter(gene %in% upreg_Heart_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_Heart_ND <- top_upreg_Heart_ND_norm %>%
  gather(colnames(top_upreg_Heart_ND_norm)[2:ncol(top_upreg_Heart_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_Heart_ND)

gathered_top_upreg_Heart_ND <- inner_join(meta2, gathered_top_upreg_Heart_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_Heart <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% upreg_Heart_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Heart DOWN
## normalized counts for top significantly downregulated genes
top_downreg_Heart_ND_norm <- normalized_counts[,c("gene",heart_samples)] %>%
  filter(gene %in% downreg_Heart_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_Heart_ND <- top_downreg_Heart_ND_norm %>%
  gather(colnames(top_downreg_Heart_ND_norm)[2:ncol(top_downreg_Heart_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_Heart_ND)

gathered_top_downreg_Heart_ND <- inner_join(meta2, gathered_top_downreg_Heart_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_Heart <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% downreg_Heart_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Lungs UP
## normalized counts for top significantly upregulated genes
top_upreg_Lungs_ND_norm <- normalized_counts[,c("gene",lungs_samples)] %>%
  filter(gene %in% upreg_Lungs_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_Lungs_ND <- top_upreg_Lungs_ND_norm %>%
  gather(colnames(top_upreg_Lungs_ND_norm)[2:ncol(top_upreg_Lungs_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_Lungs_ND)

gathered_top_upreg_Lungs_ND <- inner_join(meta2, gathered_top_upreg_Lungs_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_Lungs <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% upreg_Lungs_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Lungs DOWN
## normalized counts for top significantly downregulated genes
top_downreg_Lungs_ND_norm <- normalized_counts[,c("gene",lungs_samples)] %>%
  filter(gene %in% downreg_Lungs_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_Lungs_ND <- top_downreg_Lungs_ND_norm %>%
  gather(colnames(top_downreg_Lungs_ND_norm)[2:ncol(top_downreg_Lungs_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_Lungs_ND)

gathered_top_downreg_Lungs_ND <- inner_join(meta2, gathered_top_downreg_Lungs_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_Lungs <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% downreg_Lungs_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

## Pect UP
## normalized counts for top significantly upregulated genes
top_upreg_Pect_ND_norm <- normalized_counts[,c("gene",pect_samples)] %>%
  filter(gene %in% upreg_Pect_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_upreg_Pect_ND <- top_upreg_Pect_ND_norm %>%
  gather(colnames(top_upreg_Pect_ND_norm)[2:ncol(top_upreg_Pect_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_upreg_Pect_ND)

gathered_top_upreg_Pect_ND <- inner_join(meta2, gathered_top_upreg_Pect_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_upreg_Pect <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% upreg_Pect_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Pect DOWN
## normalized counts for top significantly downregulated genes
top_downreg_Pect_ND_norm <- normalized_counts[,c("gene",pect_samples)] %>%
  filter(gene %in% downreg_Pect_ND)


# Gathering the columns to have normalized counts to a single column
gathered_top_downreg_Pect_ND <- top_downreg_Pect_ND_norm %>%
  gather(colnames(top_downreg_Pect_ND_norm)[2:ncol(top_downreg_Pect_ND_norm)],
         key = "Sample", value = "normalized_counts")

## check the column header in the "gathered" data frame
#View(gathered_top_downreg_Pect_ND)

gathered_top_downreg_Pect_ND <- inner_join(meta2, gathered_top_downreg_Pect_ND, by="Sample")

### Extract normalized expression for significant genes from the samples 
# (in e.g. it was 4:9), and set the gene column (1) to row names
norm_sig_ND_downreg_Pect <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% downreg_Pect_ND) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


## Trying out heatmaps
## Tried a bit from but that didn't work 
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/


# my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(4,2)))
# row.names(my_sample_col) <- colnames(data_subset)

annotation_col <- data.frame(
  Tissue = factor(meta2$Tissue), 
  Metab = factor(meta$Metabolic_State))
row.names(annotation_col) <- colnames(norm_sig_ND_upreg)

pheatmap(norm_sig_ND_upreg, annotation_row = my_gene_col, annotation_col = annotation_col)


### Annotate our heatmap (optional)
annotation <- meta2 %>% 
  select(Sample, Metabolic_State, Tissue) %>% 
  data.frame(row.names = "Sample") %>%
  arrange(Metabolic_State)

### Set a color palette
heat_colors <- brewer.pal(9, "YlOrRd")

### Run pheatmap without col clusters
pheatmap(norm_sig_ND_upreg, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         cluster_cols = FALSE)


# create a `df` with the samples grouped in the same way you want to show
anno <- data.frame(SampleID = unique(gathered_top_upreg_ND$Sample), 
                   TissueState= unique(gathered_top_upreg_ND$Tissue_State))

# set rownames so that anno and your data can be matched
rownames(anno) <- gathered_top_upreg_ND$Sample


# set colours
anno_colors <- list(
  TissueState = c(Liver_N = "#1B9E77", Gut1_N = "#D95F02", Gut2_N = "#823de9", Gut3_N = "#7855ce",
                  Heart_N = "#6e6eb2", Lungs_N = "#648697", Pect_N = "#599e7c", Gut3_D = "#f94144", 
                  Gut1_D = "#f3722c", Gut2_D = "#f8961e", Gut3_T = "#f9844a", Lungs_T = "#f9c74f",
                  Lungs_D = "#90be6d", Pect_D = "#43aa8b", Heart_D = "#4d908e", Gut2_T = "#577590", 
                  Gut1_T = "#277da1", Liver_D = "#004e64", Liver_T = "#ffba08", Pect_T = "#f7b2bd", Heart_T = "#c60f7b"))


## THIS is the one I used in the lab meeting presentation on 4/21/22
# set colours
anno_colors <- list(
  Tissue = c(Liver = "#f3722c", Gut1 = "#D95F02", Gut2 = "#f8961e", Gut3 = "#577590",
             Heart = "#ffba08", Lungs = "#004e64", Pect = "#599e7c"),
  Metabolic_State = c(N = "#f9c74f", D = "#43aa8b", T = "#277da1"))

pheatmap(norm_sig_ND_upreg, 
         annotation_col = annotation, 
         annotation_colors = anno_colors,
         #color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         border_color = NA, 
         scale = "row", 
         fontsize_row = 10,
         fontsize_col = 10,
         height = 20,
         cluster_cols = F,
         legend = T,
         fontsize = 20)


annotation_col <- data.frame(
  Tissue = factor(meta2$Tissue), 
  Metab = factor(meta$Metabolic_State))

#rownames(annotation_col) = paste("Test", 1:10, sep = "")

# annotation_row <- data.frame(
#   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
# )
# rownames(annotation_row) = paste("Gene", 1:20, sep = "")

ann_colors = list(
  Tissue = c(Liver = "#1B9E77", Gut1 = "#D95F02", Gut2 = "#823de9", Gut3 = "#7855ce",
             Heart = "#6e6eb2", Lungs = "#648697", Pect = "#599e7c"),
  Metab = c(N = "red", T = "black", D = "purple")
)

#c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")

## Tissue-specific plots


# set colours
anno_colors <- list(
  # Tissue = c(Liver = "#f3722c", Gut1 = "#D95F02", Gut2 = "#f8961e", Gut3 = "#577590",
  #            Heart = "#ffba08", Lungs = "#004e64", Pect = "#599e7c"),
  Metabolic_State = c(N = "#f9c74f", T = "violet", D = "#599e7c"))

## Liver heatmap
annotation_liver <- meta2 %>%
  filter(Tissue == "Liver") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## Ordering cols by metabolic state
## USE THIS
# 1) reorder the matrix based in the annotation
liver_upND_ordered <- top_upreg_Liver_ND_norm[, rownames(annotation_liver)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(liver_upND_ordered,
                   annotation_col = annotation_liver, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   #cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Downreg Liver ND
# 1) reorder the matrix based in the annotation
liver_dnND_ordered <- norm_sig_ND_downreg_Liver[, rownames(annotation_liver)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(liver_dnND_ordered,
                   annotation_col = annotation_liver, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Heart heatmaps
## Ordering cols by metabolic state
annotation_heart <- meta2 %>%
  filter(Tissue == "Heart") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
heart_upND_ordered <- norm_sig_ND_upreg_Heart[, rownames(annotation_heart)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(heart_upND_ordered,
                   annotation_col = annotation_heart, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg heart ND
# 1) reorder the matrix based in the annotation
heart_dnND_ordered <- norm_sig_ND_downreg_Heart[, rownames(annotation_heart)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(heart_dnND_ordered,
                   annotation_col = annotation_heart, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Lungs heatmaps
## Ordering cols by metabolic state
annotation_lungs <- meta2 %>%
  filter(Tissue == "Lungs") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
lungs_upND_ordered <- norm_sig_ND_upreg_Lungs[, rownames(annotation_lungs)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(lungs_upND_ordered,
                   annotation_col = annotation_lungs, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg lungs ND
# 1) reorder the matrix based in the annotation
lungs_dnND_ordered <- norm_sig_ND_downreg_Lungs[, rownames(annotation_lungs)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(lungs_dnND_ordered,
                   annotation_col = annotation_lungs, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Pect heatmaps
## Ordering cols by metabolic state
annotation_pect <- meta2 %>%
  filter(Tissue == "Pect") %>%
  select(Sample, Metabolic_State) %>% 
  data.frame(row.names = "Sample") %>%
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)

## USE THIS
# 1) reorder the matrix based in the annotation
pect_upND_ordered <- norm_sig_ND_upreg_Pect[, rownames(annotation_pect)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(pect_upND_ordered,
                   annotation_col = annotation_pect, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)

## Downreg pect ND
# 1) reorder the matrix based in the annotation
pect_dnND_ordered <- norm_sig_ND_downreg_Pect[, rownames(annotation_pect)]

# 2) plot heatmap with no row clusters
pheatmap::pheatmap(pect_dnND_ordered,
                   annotation_col = annotation_pect, 
                   annotation_colors = anno_colors,
                   #color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = T,
                   border_color = NA, 
                   scale = "row", 
                   fontsize_row = 15,
                   fontsize_col = 15,
                   height = 20,
                   cluster_cols = F,
                   legend = T,
                   fontsize = 20)


## Volcano plot normothermy vs deep torpor
#### These are probab
ND_volcano <- EnhancedVolcano(res_ND,
                              lab = rownames(res_ND),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Deep torpor vs. Normothermy',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)


## Volcano plot normothermy vs deep torpor
NT_volcano <- EnhancedVolcano(res_NT,
                              lab = rownames(res_NT),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Transition vs. Normothermy',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)

## Volcano plot normothermy vs deep torpor
TD_volcano <- EnhancedVolcano(res_TD,
                              lab = rownames(res_TD),
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Deep Torpor vs. Transition',
                              subtitle = "Enhanced volcano, p cutoff = 0.05",
                              pCutoff = 0.05)

grid.arrange(ND_volcano, NT_volcano, TD_volcano, ncol=3)

## Use adjusted p-vals to make plots instead of nominal p-values
ND_adj_volcano <- EnhancedVolcano(res_ND, lab = rownames(res_ND),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Deep Torpor vs. Normothermy',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")
#legend=c("NS","Log2 FC","Adjusted p-value",
#        "Adjusted p-value & Log2 FC"),
#legendPosition = "bottom",
#legendLabSize = 10,
#legendIconSize = 3.0)

NT_adj_volcano <- EnhancedVolcano(res_NT, lab = rownames(res_NT),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Transition vs. Normothermy',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")

TD_adj_volcano <- EnhancedVolcano(res_TD, lab = rownames(res_TD),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  #xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(-Log[10]~adjusted~italic(P)),
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  #transcriptLabSize = 3.0,
                                  colAlpha = 1,
                                  title = 'Deep Torpor vs. Transition',
                                  subtitle = "Enhanced volcano, adj p cutoff = 0.05")

grid.arrange(ND_adj_volcano, NT_adj_volcano, TD_adj_volcano, ncol=3)

ND_adj_volcano
NT_adj_volcano
TD_adj_volcano


## Tissue specific
## Use adjusted p-vals to make plots instead of nominal p-values
ND_adj_volcano_liver <- EnhancedVolcano(res_Liver_ND, lab = rownames(res_Liver_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'Liver: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_liver


ND_adj_volcano_heart <- EnhancedVolcano(res_Heart_ND, lab = rownames(res_Heart_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'Heart: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_heart

ND_adj_volcano_lungs <- EnhancedVolcano(res_Lungs_ND, lab = rownames(res_Lungs_ND),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        #xlab = bquote(~Log[2]~ "fold change"),
                                        ylab = bquote(-Log[10]~adjusted~italic(P)),
                                        pCutoff = 0.05,
                                        FCcutoff = 1.0,
                                        #transcriptLabSize = 3.0,
                                        colAlpha = 1,
                                        title = 'Lungs: Deep Torpor vs. Normothermy',
                                        subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_lungs

ND_adj_volcano_pect <- EnhancedVolcano(res_Pect_ND, lab = rownames(res_Pect_ND),
                                       x = 'log2FoldChange',
                                       y = 'padj',
                                       #xlab = bquote(~Log[2]~ "fold change"),
                                       ylab = bquote(-Log[10]~adjusted~italic(P)),
                                       pCutoff = 0.05,
                                       FCcutoff = 1.0,
                                       #transcriptLabSize = 3.0,
                                       colAlpha = 1,
                                       title = 'Pect: Deep Torpor vs. Normothermy',
                                       subtitle = "Enhanced volcano, adj p cutoff = 0.05")
ND_adj_volcano_pect


## geom_density
datlong %>%
  filter(gene == "RSRP1") %>%
  ggplot(., aes(x=counts, fill=Metabolic_State)) +
  geom_density(alpha=0.5) + my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states")

# Boxplot/Violin plot
datlong %>%
  filter(gene == "RSRP1") %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("RSRP1 gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

#old, manually selected from each tissue type
#genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')

## pull together all the top upreg and downreg genes from each pairwise comparison
## and remove the ones which are not annotated
top_genes <- c(upreg_ND, upreg_NT, upreg_TD, downreg_ND, downreg_NT, downreg_TD)
top_genes_annotated <- top_genes[!grepl("LOC", top_genes)]

upreg_ND[!grepl("LOC", upreg_ND)]
upreg_NT[!grepl("LOC", upreg_NT)]
upreg_TD[!grepl("LOC", upreg_TD)]

downreg_ND[!grepl("LOC", downreg_ND)]
downreg_NT[!grepl("LOC", downreg_NT)]
downreg_TD[!grepl("LOC", downreg_TD)]


upreg_Liver_ND[!grepl("LOC", upreg_Liver_ND)]
upreg_Heart_ND[!grepl("LOC", upreg_Heart_ND)]
upreg_Lungs_ND[!grepl("LOC", upreg_Lungs_ND)]
upreg_Pect_ND[!grepl("LOC", upreg_Pect_ND)]

downreg_Liver_ND[!grepl("LOC", downreg_Liver_ND)]
downreg_Heart_ND[!grepl("LOC", downreg_Heart_ND)]
downreg_Lungs_ND[!grepl("LOC", downreg_Lungs_ND)]
downreg_Pect_ND[!grepl("LOC", downreg_Pect_ND)]


## Trying out enrichR package
enrichr(upreg_ND, databases = c("BioPlanet_2019", "KEGG_2021_Human"))
enrichr(downreg_ND, databases = c("BioPlanet_2019", "KEGG_2021_Human"))

if (getOption("enrichR.live")) {
  dbs <- listEnrichrDbs()
  enrichRLive <- TRUE
  if (is.null(dbs)) enrichRLive <- FALSE
  dbs <- c("BioPlanet_2019", "KEGG_2021_Human")
  enriched <- enrichr(upreg_ND, dbs)
  # Plot top 20 GO-BP results ordered by P-value
  if (enrichRLive) {
    plotEnrich(enriched[[2]], showTerms = 10, numChar = 50, y = "Ratio",
               orderBy = "P.value")
  }
}


head(upreg_ND)

## Make numeric version of metabolic states
datlong$MetabState_numeric <- NA
datlong$MetabState_numeric[datlong$Metabolic_State=="N"] <- 1
datlong$MetabState_numeric[datlong$Metabolic_State=="T"] <- 2
datlong$MetabState_numeric[datlong$Metabolic_State=="D"] <- 3

datlong %>%
  filter(gene %in% top_genes_annotated) %>%
  ggplot(., aes(x = Tissue, y = gene, fill = counts)) +
  geom_tile() + my_theme +
  scale_fill_gradient(low = 'white', high = 'red') + facet_grid(.~Metabolic_State, scales = "free")

clock_all <- c("CLOCK", "CRY1", "CRY2", "PER2", "PER3")
clock_col <- c("black", "pink", "red", "green", "blue")
clockgenes1 <- c('CRY1', 'CRY2')
clockgenes2 <- c('PER2', 'PER3')


datlong %>%
  filter(gene %in% clock_all, Tissue=="Heart") %>%
  ggplot(., aes(y=log(counts), x=MetabState_numeric)) +
  geom_smooth(aes(col=gene)) + 
  #geom_violin() +
  my_theme + #facet_wrap(.~Tissue, scales = "free") +
  ggtitle("Clock genes' expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Normothermy", "Transition", "Deep torpor")) +
  ylab("Gene counts") + xlab("Metabolic state")

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }
preproc <- normalize(datlong$counts[datlong$gene =="CLOCK" & datlong$Tissue=="Heart"])

datlong$NormCounts <- NA
datlong <- datlong %>% group_by(gene,Tissue) %>% mutate(NormCounts = counts/max(counts))

datlong %>%
  filter(!is.na(gene), gene=='CLOCK') %>%
  ggplot(., aes(y=NormCounts, x=MetabState_numeric)) +
  geom_line(aes(col=Tissue), method = 'lm', stat = "smooth", size=2) + 
  #geom_violin() +
  my_theme2 + #facet_wrap(.~Tissue, scales = "free") +
  ggtitle("CLOCK expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #scale_color_viridis_d(begin = 0, end=1,option = "C", direction = -1) + 
  scale_color_manual(values = mycols) +
  theme(axis.text.x = element_text(size=15)) +
  scale_x_continuous(labels = c("Normothermy", "Transition", "Deep torpor"), breaks = c(1,2,3)) +
  ylab("Normalized gene count") + xlab("Metabolic state")

datlong %>%
  filter(!is.na(gene), gene  == 'CRY2') %>%
  ggplot(., aes(y=NormCounts, x=MetabState_numeric)) +
  geom_line(aes(col=Tissue), method = 'lm', stat = "smooth", size=2) + 
  #geom_violin() +
  my_theme2 + facet_wrap(.~gene, scales = "free") +
  ggtitle("CRY2 expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.3)) +
  #scale_color_viridis_d(begin = 0, end=1,option = "C", direction = -1) + 
  scale_color_manual(values = mycols) +
  theme(axis.text.x = element_text(size=15)) +
  scale_x_continuous(labels = c("Normothermy", "Transition", "Deep torpor"), breaks = c(1,2,3)) +
  ylab("Normalized gene count") + xlab("Metabolic state")


datlong %>%
  filter(!is.na(gene), gene %in% clockgenes1) %>%
  ggplot(., aes(y=NormCounts, x=MetabState_numeric)) +
  geom_line(aes(col=gene),method = 'lm', stat = "smooth", size=2) + 
  #geom_violin() +
  my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("Clock genes' expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = clock_col) +
  scale_x_discrete(labels = c("Normothermy", "Transition", "Deep torpor")) +
  ylab("Gene counts") + xlab("Metabolic state")

datlong %>%
  filter(gene == 'CLOCK', Tissue %in% c("Heart", "Lungs")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("CLOCK gene expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Normothermy", "Transition", "Deep torpor")) +
  ylab("Gene counts") + xlab("Metabolic state")


datlong %>%
  filter(gene == clockgenes1, Tissue %in% c("Heart", "Lungs")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_grid(gene~Tissue, scales = "free") +
  ggtitle("Clock genes' expression across tissues and metabolic states") +
  scale_x_discrete(labels = c("Normothermy", "Transition", "Deep torpor")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")


datlong %>%
  filter(gene == clockgenes2, Tissue %in% c("Heart", "Lungs")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + 
  facet_grid(gene~Tissue, scales = "free") +
  ggtitle("Clock genes' expression across tissues and metabolic states") +
  scale_x_discrete(labels = c("Normothermy", "Transition", "Deep torpor")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")


## Metabolism genes from Figure 3 in https://www.nature.com/articles/s41598-018-31506-2/figures/3
# This paper also has clock genes
metabgenes <- c('PPARA', 'SIRT1', 'LEPR')
# Boxplot/Violin plot of metabolism genes
datlong %>%
  filter(gene == metabgenes) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(gene~Tissue, scales = "free") +
  ggtitle("Metabolism genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

# Boxplot/Violin plot of metabolism genes just gut
datlong %>%
  filter(gene == metabgenes & Tissue %in% c("Gut1", "Gut2", "Gut3")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(gene~Tissue, scales = "free") +
  ggtitle("Metabolism genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")

# Boxplot/Violin plot of metabolism genes just non-gut
datlong %>%
  filter(gene == metabgenes & Tissue %in% c("Liver", "Pect", "Heart", "Lungs")) %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  #geom_violin() +
  my_theme + facet_wrap(gene~Tissue, scales = "free") +
  ggtitle("Metabolism genes expression across tissues and metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state")


## Examples
genes.of.interest <- c('NR1D1', 'SGK1', 'RSRP1', 'DUSP1', 'NME5', 'TMEM39B', 'BHLHE40', 'LOC115598688')
foldlong %>%
  filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')

foldlong %>%
  filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


sig_ND <-  res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(log2FoldChange>1 | log2FoldChange < -1) %>%
  filter(padj < 0.05 & baseMean > 10) %>%
  left_join(., meta2, by="Sample") %>%
  mutate(Metabolic_State = 
           fct_relevel(Metabolic_State, 
                       "N", "T", "D")) 
arrange(padj) 

ggplot(sig_ND, aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


res_ND %>%
  data.frame() %>%
  rownames_to_column(var="gene")
filter(Measure %in% c('D_vs_N_log2FoldChange', 'T_vs_N_log2FoldChange', 'T_vs_D_log2FoldChange')) %>%
  ggplot(., aes(x = Measure, y = gene, fill = value)) +
  geom_tile() + my_theme + facet_grid(.~Tissue) +
  scale_fill_gradient(low = 'yellow', high = 'red')


sig_top_all_comparisons  <- datlong %>%
  filter(gene %in% top_genes_annotated)

sig_top_all_comparisons <- normalized_counts[,c(1, 2:120)] %>% 
  filter(gene %in% top_genes_annotated) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

pheatmap(sig_top_all_comparisons, 
         annotation_col = annotation, 
         annotation_colors = anno_colors,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         border_color = NA, 
         scale = "row", 
         fontsize_row = 15,
         fontsize_col = 10,
         height = 20,
         cluster_cols = F,
         legend = T,
         fontsize = 20)


