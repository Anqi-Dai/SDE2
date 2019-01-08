# process the sh*t out of the gtf using my freaking tidyverse skills

library(tidyverse)

setwd('~/SDE2/analysis/09_mitochondria/scripts/')

gtf <- read_tsv('../../../reference/v29/gencode.v29.annotation.gtf' , comment = '#', col_names = F) %>%
  filter(X3 == 'exon') %>%
  dplyr::select(X9)

col9 <- gtf %>%
  separate(X9, into = c('gene_id','transcript_id','gene_type','gene_name','transcript_type','transcript_name','exon_number','exon_id','level','transcript_support_level','tag','havana_gene','havana_transcript'), sep = '; ')

col9_sub <- col9 %>%
  dplyr::select(gene_id,transcript_id, exon_number,exon_id) 

#group by the transcript to have the exon number and
#then group by the geneID to have he median exon number of that

col9_sum <- col9_sub %>%
  group_by(gene_id,transcript_id) %>%
  summarise(exon_total_num = n()) 

col9_median <- col9_sum %>%
  group_by(gene_id) %>%
  summarise(med_exon_num = median(exon_total_num))

med_exon_number <- col9_median %>%
  mutate(gene_id = str_replace_all(gene_id, 'gene_id "','')) %>%
  mutate(gene_id = str_replace_all(gene_id, '"','')) %>%
  mutate(med_exon_num = ceiling(med_exon_num)) %>%
  mutate(gene_id = str_replace_all(string = gene_id, pattern = '\\..*$', replacement = ''))  %>%
  mutate(med_intron_num = med_exon_num - 1)

# write out the table 
write_csv(med_exon_number %>%
            dplyr::select(gene_id, med_intron_num), '../data/gene_median_intron_num.csv')


