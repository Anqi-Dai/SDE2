# this script will do the filtering on the salmon counts
# There has to be at least 1 read normalized count in all three replicates

library(Biobase)
library(tidyverse)
library(VennDiagram)
library(ggpubr)
setwd(dir = '~/Downloads/git_repo/SDE2/analysis/09_mitochondria/scripts/')

# subset to have each SAMPLE together
eset <- read_rds('../../08_findingGenes/data/Norm_salmon_cts_replicate1_eset_unfiltered.RDS')
eset14 <- eset[,grepl('14', sampleNames(eset))]
eset89 <- eset[,grepl('89', sampleNames(eset))]
esetCtrl <- eset[,grepl('Control', sampleNames(eset))]


# the filtering is to remove genes with any zero in the three replicates
removeAnyZero <- function(eset) {
  keep_idx = apply(exprs(eset), 1, function(row) all(row !=0 ))
  return(eset[keep_idx,])
}

# the filtering
eset14_f <- removeAnyZero(eset14)
eset89_f <- removeAnyZero(eset89)
esetCtrl_f <- removeAnyZero(esetCtrl)

# overlap with genes that are signicant in any of the 3 AS tool runs

# loading the sig genes
sig_gene <- read_table('../../08_findingGenes/data/union_sig_genes_from_3_analysis.txt', col_names = F)





# plot to check the how is the normalized counts like
# take the control for example

ctrl_haveE <- setdiff(featureNames(esetCtrl_f), sig_gene$X1)
ctrl_intersect <- intersect(featureNames(esetCtrl_f), sig_gene$X1)
ctrl_sig <- setdiff(sig_gene$X1, featureNames(esetCtrl_f) )


# could be the difference in the ensembl gene IDs
length(featureNames(eset)[!sig_gene$X1 %in% featureNames(eset)])
# 45 sig IDs are not in the feature names of the whole eset.
# the different annotation each tool uses could also affect the ...


 
# get the freaking counts
DF <- exprs(esetCtrl[ctrl_haveE,]) %>%
  as.data.frame %>%
  rownames_to_column %>%
  mutate(Group = 'haveE')

DF_2 <- exprs(esetCtrl[ctrl_intersect,]) %>%
  as.data.frame %>%
  rownames_to_column %>%
  mutate(Group = 'Intersect')


DF_3 <- exprs(esetCtrl[ctrl_sig[ctrl_sig %in% featureNames(esetCtrl)],]) %>%
  as.data.frame %>%
  rownames_to_column %>%
  mutate(Group = 'Sig')


# bind the three up and plot
df  <- bind_rows(DF, DF_2) 
df  <- bind_rows(df, DF_3) 

df <- df %>%
  mutate(Group = factor(Group)) %>%
  gather(key = Replicate, value = value, RF_Control_1:RF_Control_3) %>%
  mutate(Replicate = factor(Replicate)) 

df %>%
  ggplot(aes(x = Group, y = value)) +
  geom_boxplot()+
  scale_y_log10() + 
  facet_wrap(Replicate~.)+
  labs(y =  'log10(normalized counts)',
       title = 'Genes expression in RF-control three replicates') +
  ggsave('../figs/Genes expression in RF-control three replicates.png')


# the detk-filter
# the pheno table
pheno <- pData(eset)
rownames(pheno) <- gsub('_','-', rownames(pheno))
write.csv(pheno, '../../08_findingGenes/data/pheno_data.csv', quote = F)


# running detk filter to filter the normalized count matrix to filter and reduce the pink circle a bit
#  detk-filter '(nonzero(Status[KD1]) >= 3 and nonzero(Status[CTRL]) >= 3) or (nonzero(Status[KD2]) >= 3 and nonzero(Status[CTRL]) >= 3)' Norm_salmon_cts_replicate1_matrix.csv pheno_data.csv -o detk_fil_Norm_salmon_cts_replicate1_matrix.csv

detkF <- read_csv('../../08_findingGenes/data/detk_fil_Norm_salmon_cts_replicate1_matrix.csv') %>%
  mutate(Name = str_replace_all(string = Name, pattern = '\\..*$', replacement = '')) 
  
# below is the filtered sig genes.
Fsig <- intersect(detkF$Name, sig_gene$X1)

# redo the venn diagram
# between 14 and sig
diff_list <- list(haveExpre14 = featureNames(eset14_f),
                  significant = Fsig)
fill <- c("light blue", "pink")
size  <- rep(0.5, 2)
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "haveExpression in 14 VS sig genes");
png('../figs/haveExpression in 14 VS sig genes filtered.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()

# between 89 and sig 
diff_list <- list(haveExpre89 = featureNames(eset89_f),
                  significant = Fsig)
fill <- c("light blue", "pink")
size  <- rep(0.5, 2)
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "haveExpression in 89 VS sig genes");
png('../figs/haveExpression in 89 VS sig genes filtered.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()

# between control and sig 
diff_list <- list(haveExpreCtrl = featureNames(esetCtrl_f),
                  significant = Fsig)
fill <- c("light blue", "pink")
size  <- rep(0.5, 2)
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "haveExpression in Control VS sig genes");
png('../figs/haveExpression in Control VS sig genes filtered.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()



