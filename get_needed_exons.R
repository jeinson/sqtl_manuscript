library(stringr)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
file = args[2]

data = read.table(file, header=T, sep='\t')

# Gene name extraction

gene_name = unlist(str_split(data$phenotype_id[j],pattern = "_"))[1]
gene_name
bash = c("bash retrieve_nucleotide_sequence_from_genename.sh -g", gene_name)
bash =  paste(bash, collapse=" ")
print(bash)
## all gene exons extraction
system(bash)
