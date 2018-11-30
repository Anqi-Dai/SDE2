# This script is for the filtering of the IRfinder result:
# input: raw result from the IRfinder. two comparisons: 14-control and 89-control
# output: a combined DF of both 14 and 89, including the TRUE significant that are significant
# locus in both, and FALSE significant that are significant in neither but are overlapped
# between the two
# (command + shift + C for multiline comments)

setwd('../angelD/SDE2/analysis/05_IRfinder/script/')

library(tidyverse)

##################################################################################

# clean the table 
# the non-filtered raw result of the 14 and 89
data_path <- '../data/rawIRfinderResult'
fns <- file.path(data_path, list.files(data_path, pattern = 'txt'))

# format the table like below (generate the locus ID, change the colnames, format the geneID string, separate that one column to three, generate the padj value.)
raw_df_list <- lapply(fns, function(fn){
  res = read_tsv(fn,comment = '#',col_names = T) %>%
    mutate(Locus =paste(Chr,paste(Start, End, sep = '-'), sep = ':' ) ) 
  
  colnames(res) <- colnames(res) %>%
    gsub('-', '_', .)
    
  res = res %>%
    separate('Intron_GeneName/GeneID', into = c('Symbol', 'GeneID', 'Category'), sep = '/') %>%
    # get rid of the . and digits after that in GeneID
    mutate(GeneID = str_replace_all(string = GeneID, pattern = '\\..*$', replacement = '')) %>%
    mutate(padj = p.adjust(p_diff, method = 'BH')) 
})

##################################################################################

# filter the result like below:
# * Remove records marked with "known-exon". "known-exon" indicates that there are known transcripts both including and splicing out this region. IRFinder quantifies the relative abundance of these occurrences.
# * Remove records marked with "MinorIsoform" either in experiment or control. "MinorIsoform" indicates IRFinder has detected there is additional alternate-splicing occurring, as such the detection of intron-retention may be confounded by other changes in alternate splicing.
# * Adjust the P-value with BH method (The Audic and Claverie test is used for low replicates situations (<= 3 replicates). The assumption is that the IR ratio (log transformed) follows a Poisson distribution)
# * Select the records with padj < 0.05

fil_df_list <- lapply(raw_df_list, function(df) {
  res = df  %>%
    filter(! Category == 'known-exon') %>%
    filter( A_IRok != 'MinorIsoform' &  B_IRok != 'MinorIsoform' ) %>%
    filter(padj < 0.05) 
  return(res)
} )  


##################################################################################

# TRUE: sig in both
# FALSE: sig in neither but overlapped in both

true_locus <- intersect(fil_df_list[[1]]$Locus, fil_df_list[[2]]$Locus)
overlapped_locus <- intersect(raw_df_list[[1]]$Locus, raw_df_list[[2]]$Locus)
false_locus <- setdiff(overlapped_locus, fil_df_list[[1]]$Locus) 
false_locus <- setdiff(false_locus, fil_df_list[[2]]$Locus)

# subset the df to include only the true and false locus

subset_df_list <- lapply(raw_df_list, function(df){
  res = df %>%
    filter(Locus %in% true_locus | Locus %in% false_locus) 
})

RES <- bind_rows(subset_df_list[[1]] %>%
                   mutate(Group = 'RF14'),
                 subset_df_list[[2]] %>%
                   mutate(Group = 'RF89')) %>%
  mutate(Group =  factor(Group)) %>%
  mutate(Significance = ifelse(Locus %in% true_locus, TRUE, FALSE))


# output the df
write_csv(RES, '../data/IRfinder.res.filtered.14n89.csv')

