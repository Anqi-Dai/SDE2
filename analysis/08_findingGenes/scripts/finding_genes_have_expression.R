## Finding genes that have expression (should be around 20,000)
setwd('~/Downloads/git_repo/SDE2/analysis/08_findingGenes/scripts/')

### BUILDING THE ESET

# loading the normalized salmon counts
cts <- read_csv('../data/Norm_salmon_cts_replicate1_matrix.csv', col_types = cols(Name = col_character())) %>%
  rename_all(
    funs(
      stringr::str_replace_all(., '-', '_') 
    )) %>%
  mutate(Name = str_replace_all(string = Name, pattern = '\\..*$', replacement = '')) %>%
  filter(!duplicated(Name)) %>%
  column_to_rownames(var = 'Name')  


# finding the feature data of the genes
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


genes <- rownames(cts)


symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name',
                               'chromosome_name','start_position','end_position', 'description'),
                values = genes, 
                mart = mart)

featureInfo <- symbol %>%
  filter(!duplicated(ensembl_gene_id)) %>%
  column_to_rownames(var = 'ensembl_gene_id')  


# the pheno data
pheno <- data.frame(
  row.names = colnames(cts),
  Status = c(rep('KD1',3), rep('KD2',3),rep('CTRL',3))
)

# build eset
library(Biobase)

phenoData <- new("AnnotatedDataFrame", data = pheno)
featureData <- new('AnnotatedDataFrame', data = featureInfo)
CTS <- cts[rownames(cts) %in% rownames(featureInfo),]

# build the eset object
rawEset <- ExpressionSet(assayData=as.matrix(CTS), 
                         phenoData=phenoData, 
                         featureData= featureData)

write_rds(rawEset, '../data/Norm_salmon_cts_replicate1_eset_unfiltered.RDS')


### FILTERING ON THE ESET

## remove those genes without at least 1 read per million in at least 'n' samples 
## n = least amount of samples in a condition (4 in example)
removeLowExpression <- function(eset, class_id)
{
  groups <- pData(eset)[,class_id]
  min.samples <- min( sapply(levels(groups), function(x){length(which(groups %in% x))}) )
  rpm <- colSums(exprs(eset))/1000000
  filter_ind <- t(apply(exprs(eset), 1,function(x) {x >rpm}))
  filter_ind_rowsums <- apply(filter_ind, 1, sum)
  return(eset[filter_ind_rowsums > min.samples,])
}

feset <- removeLowExpression(eset = rawEset, class_id = "Status")


write_rds(feset, '../data/FILTERED_salmon_cts_replicate1_eset.RDS')

# the remaining number of features
length(featureNames(feset))



