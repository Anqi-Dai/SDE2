# This script is to do analysis on RI locus and append information of below to the df
# * GC content 
# * MaxEnt score

library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(seqinr)


# load the data
RES <- read_csv('../data/IRfinder.res.filtered.14n89.analysis1.csv', col_types = cols(Group = col_factor()))

##################################################################################

# 1. The GC content of the RI

RES <- RES %>%
  mutate(Start = Start + 1, Chr = factor(Chr))

genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlevelsStyle(genome) = "UCSC"
grange = makeGRangesFromDataFrame(RES, seqnames.field=c("Chr"), start.field="Start",
                                  end.field=c("End"), strand.field="Direction")
seq = BSgenome::getSeq(genome, grange)

RES <- RES %>%
  mutate(Seq = as.character(seq))  %>%
  mutate(GC = stringr::str_count(Seq, "G") + stringr::str_count(Seq, "C")) %>%
  mutate(GC_per = GC / stringr::str_length(Seq) * 100)

##################################################################################

# 2. The score of the 5' ss and 3' ss  according to the Maximum Entropy Model

# check whether the splice site is canonical ss

RES <- RES %>%
  mutate(Donor = stringr::str_sub(Seq, 1, 2),
         Acceptor = stringr::str_sub(Seq, -2),
         Splicesite_pair = paste(Donor, Acceptor, sep = '-')) 

# * Score 5' splice site: Each sequence must be 9 bases long. 3 bases in exon and 6 bases in intron
# * Score 3' splice site: Each sequence must be 23 bases long. 20 bases in the intron and 3 base in the exon

RES <-  RES %>%
  mutate(ss5_start = Start - 3,
         ss5_end = Start + 5,
         ss3_start = End -20,
         ss3_end = End + 2 )

grange5 <-  makeGRangesFromDataFrame(RES, 
                                     seqnames.field=c("Chr"),  
                                     start.field="ss5_start", 
                                     end.field=c("ss5_end"), 
                                     strand.field="Direction")


grange3 <-  makeGRangesFromDataFrame(RES, 
                                     seqnames.field=c("Chr"),  
                                     start.field="ss3_start", 
                                     end.field=c("ss3_end"), 
                                     strand.field="Direction")

RES <- RES %>%
  mutate(ss5_seq = as.character(BSgenome::getSeq(genome, grange5)),
         ss3_seq = as.character(BSgenome::getSeq(genome, grange3)),
         id = paste(Chr,paste(Start, End, sep = '-'), sep = ':' ) )

# write out to fasta file
write.fasta(sequences = as.list(RES$ss5_seq), 
            names = RES$id,
            file.out =  '../data/ss5.seq.fasta')


write.fasta(sequences = as.list(RES$ss3_seq), 
            names = RES$id,
            file.out =  '../data/ss3.seq.fasta')

# get the results from online page http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html and append the score to the df

maxent_res <- lapply(file.path('../output', list.files('../output/', pattern = 'txt$')), function(fn) {
  res = read_delim(fn,comment = '>',delim = '\t',
                   col_names = F) %>%
    mutate(score = as.numeric(stringr::str_replace_all(string = X2, pattern = 'MAXENT: ', replacement = '')))
})

RES$ss5_score <- maxent_res[[2]]$score
RES$ss3_score <- maxent_res[[1]]$score

##################################################################################

# write out the table

write_csv(RES, '../data/IRfinder.res.filtered.14n89.analysis2.csv')


# select certain columns and sort and output the table


OUT <-  RES %>%
  dplyr::select(Locus, GeneID, Symbol, description, Direction, A_IRratio, B_IRratio,A_IntronDepth, B_IntronDepth, p_diff, p_increased, p_decreased, padj, Group, Significance, RI_prop, GC_per, Splicesite_pair, ss5_score, ss3_score, RFExper, RFControl)  %>%
  filter(Significance == TRUE  ) %>%
  mutate(absDelta = abs(A_IRratio - B_IRratio)) %>%
  arrange(desc(A_IntronDepth), desc(absDelta)) %>%
  write_csv('../output/IRfinder.true.locus.sorted.both.comparison.csv')
