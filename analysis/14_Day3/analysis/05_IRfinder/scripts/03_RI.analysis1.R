# This script is to do analysis on RI locus and append information of below to the df
# * the length of RI
# * the TPM of the genes
# * the locus of the RI at the genes

library(tidyverse)

# load the data
RES <- read_csv('../data/IRfinder.res.filtered.14n89.csv')

##################################################################################

# 1. The length of the RI (abs difference between the start and end position of the RI)

RES <- RES %>%
  mutate(RI_len = abs(End-Start))

##################################################################################

# 2. The TPM of the genes. Using TPM result from Whippet

cts <- read_csv('../data/meanTPM.genes.whippet.csv', col_types = 'cddd')

# join the tpm of the genes in 14 and 89 separately
RES14 <- RES %>%
  filter(Group == 'RF14') %>%
  left_join(cts %>%
              dplyr::select(GeneID, RF14, RFControl), by = 'GeneID') %>%
  rename(RFExper = RF14)
  

RES89 <- RES %>%
  filter(Group == 'RF89') %>%
  left_join(cts %>%
              dplyr::select(GeneID, RF89, RFControl), by = 'GeneID')%>%
  rename(RFExper = RF89)

RES <- bind_rows(RES14, RES89)

##################################################################################

# 3. The locus of the RI at the genes (Compare the proportion of the middle point of RI at the whole gene length)(The closer to 1, the closer to 3' end).

# get the start and end and description info of those genes 
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# the ensembl_gene_id of the genes that you are interested in 
genes <- RES %>%
  distinct(GeneID) %>%
  pull(GeneID)

# retrieve the feature information of the listed attributes
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c('ensembl_gene_id',
                               'chromosome_name',
                               'start_position',
                               'end_position', 
                               'description'),
                
                values = genes, 
                mart = mart)

# there are genes that I currently don't have feature info in this table
RES <- RES %>%
  left_join(symbol %>%
               rename(GeneID = ensembl_gene_id), by = 'GeneID')  %>%
  mutate(midRI = floor((Start+End)/2),
         geneLen = end_position - start_position,
         RI_prop = (midRI - start_position)/geneLen)

##################################################################################

# write out the table

write_csv(RES, '../data/IRfinder.res.filtered.14n89.analysis1.csv')

