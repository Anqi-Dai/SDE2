# this script will process the GTF file and get the median intron number of the genes
library(refGenome)
library(Rsubread)

setwd('~/Downloads/git_repo/newFlynn1802/flynn_altstatus/reference/')

# create ensemblGenome object for storing Ensembl genomic annotation data(It has to be the ensembl gtf file not the gencode one!)
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "Homo_sapiens.GRCh38.94.chr.gtf")

# now the annotation is read into the ens obj
class(ens)

# extract in the primary assembly
enpa <- extractSeqids(ens,ensPrimAssembly())

#group by the exon to have the exon number and
#then group by the geneID to have he median exon number of that

# Extract Exons
enex <- refExons(enpa)
gt <- getExonData(enex)
DF <- enex
  dplyr::select(id, seqid, start, end ,exon_id, transcript_id, exon_number, gene_name,gene_id)
