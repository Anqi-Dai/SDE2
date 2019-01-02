# relax FDR to 0.25 and get a union of the unique genes that have the locus that are significant from all three tool runs.
library(tidyverse)
# set the wd to be the SDE2 project folder
setwd("/Users/angelD/Downloads/git_repo/SDE2")

################################################################################
# THE WHIPPET SIG GENES (absDeltaPsi > 0.1 & **Probability > 0.8** & Psi_A > 0.1)
################################################################################

RF14.ctrl <- read_tsv('analysis/03_ASprofile/data/whippetResult/RF14_vs_cntrl_output.diff') %>%
  mutate(absDeltaPsi = abs(DeltaPsi)) %>%
  mutate(Gene = str_replace_all(string = Gene, pattern = '\\..+$', replacement = '')) %>%
  mutate(Group = 'RF14.ctrl') 


RF89.ctrl <- read_tsv('analysis/03_ASprofile/data/whippetResult/RF89_vs_cntrl_output.diff') %>%
  mutate(absDeltaPsi = abs(DeltaPsi)) %>%
  mutate(Gene = str_replace_all(string = Gene, pattern = '\\..+$', replacement = '')) %>%
  mutate(Group = 'RF89.ctrl') 

ALL <- bind_rows(RF14.ctrl, RF89.ctrl) %>%
  mutate(Group = factor(Group))

ALL.sig <- ALL %>%
  filter(absDeltaPsi > 0.1 & Probability > 0.8 & Psi_A > 0.1) 

## The unique genes of those locus

sig_genes_whippet <- ALL.sig %>%
  filter(!duplicated(Gene)) %>%
  pull(Gene)

length(sig_genes_whippet)
  

################################################################################
# THE IRFINDER SIG GENES (FDR < 0.25)
################################################################################

data_path <- 'analysis/05_IRfinder/data/rawIRfinderResult'
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
# * Remove records marked with "known-exon". 
# * Remove records marked with "MinorIsoform" either in experiment or control. 
# * Adjust the P-value with BH method 
# * Select the records with padj < 0.25

true_df_list <- lapply(raw_df_list, function(df) {
  res = df  %>%
    filter(! Category == 'known-exon') %>%
    filter( A_IRok != 'MinorIsoform' &  B_IRok != 'MinorIsoform' ) %>%
    filter(padj < 0.25) 
  return(res)
} )

false_df_list <- lapply(raw_df_list, function(df) {
  res = df  %>%
    filter(! Category == 'known-exon') %>%
    filter( A_IRok != 'MinorIsoform' &  B_IRok != 'MinorIsoform' ) %>%
    filter(padj >= 0.25) 
  return(res)
} )

##################################################################################

# TRUE: sig in both
# FALSE: sig in neither but overlapped in both

true_locus <- intersect(true_df_list[[1]]$Locus, true_df_list[[2]]$Locus)
false_locus <- intersect(false_df_list[[1]]$Locus, false_df_list[[2]]$Locus)

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

## The unique genes of those sig locus

sig_genes_IRfinder <- RES %>%
  filter(Significance == TRUE) %>%
  filter(!duplicated(GeneID)) %>%
  pull(GeneID)

length(sig_genes_IRfinder)

################################################################################
# THE RMATS SIG GENES (FDR < 0.25)
################################################################################

path <- 'analysis/06_rMATS/data/rMATS_result/RF14_control'
fns <- file.path(path, list.files(path, 'JCEC.txt'))
JCEC_sig_RF14 <- lapply(fns, function(fn){
  ret = read_tsv(fn, col_types = cols(ID = col_number())) %>%
    filter(FDR < .25) %>%
    mutate(Type = fn) %>%
    mutate(Type = str_replace_all(Type, '../data/rMATS_result/RF14_control/','')) %>%
    mutate(Type = str_replace_all(Type, '.MATS.JCEC.txt','')) %>%
    mutate(Position = paste(.[[6]], .[[7]], sep = '-')) %>%
    mutate(Locus = paste(chr, Position,sep = ':')) %>%
    dplyr::select(ID, Locus, GeneID, Type, IncLevelDifference, FDR) %>%
    mutate(GeneID = str_replace_all(string = GeneID, pattern = '\\..*$', replacement = '')) 
  return(ret)
})  %>% bind_rows()


path <- 'analysis/06_rMATS/data/rMATS_result/RF89_control'
fns <- file.path(path, list.files(path, 'JCEC.txt'))
JCEC_sig_RF89 <- lapply(fns, function(fn){
  ret = read_tsv(fn, col_types = cols(ID = col_number())) %>%
    filter(FDR < .25) %>%
    mutate(Type = fn) %>%
    mutate(Type = str_replace_all(Type, '../data/rMATS_result/RF89_control/','')) %>%
    mutate(Type = str_replace_all(Type, '.MATS.JCEC.txt','')) %>%
    mutate(Position = paste(.[[6]], .[[7]], sep = '-')) %>%
    mutate(Locus = paste(chr, Position,sep = ':')) %>%
    dplyr::select(ID, Locus, GeneID, Type, IncLevelDifference, FDR) %>%
    mutate(GeneID = str_replace_all(string = GeneID, pattern = '\\..*$', replacement = '')) 
  return(ret)
})  %>% bind_rows()

sig_genes_rMATS <- bind_rows(JCEC_sig_RF14, JCEC_sig_RF89) %>%
  filter(!duplicated(GeneID)) %>%
  pull(GeneID)

length(sig_genes_rMATS)

################################################################################
# A Union of those genes
################################################################################

union_sig_genes <- union(sig_genes_rMATS, union(sig_genes_whippet, sig_genes_IRfinder))

## output those genes
write.table(union_sig_genes, 'analysis/08_findingGenes/data/union_sig_genes_from_3_analysis.txt', quote = F, col.names = F, row.names = F)
