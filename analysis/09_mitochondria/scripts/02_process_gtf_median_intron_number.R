# this script will process the GTF file and get the median intron number of the genes
library(refGenome)
library(Rsubread)

setwd('~/Downloads/git_repo/newFlynn1802/flynn_altstatus/reference/')

# create ensemblGenome object for storing Ensembl genomic annotation data(It has to be the ensembl gtf file not the gencode one!)
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "Homo_sapiens.GRCh38.94.chr.gtf")

class(ens)

#group by the exon to have the exon number and
#then group by the geneID to have he median exon number of that

tableSeqids(ens)
