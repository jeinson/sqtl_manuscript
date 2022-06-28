library(stringr)

args <- commandArgs(trailingOnly = TRUE)
j = as.numeric(args[1])
file = args[2]

data = read.table(file, header=T, sep='\t')
## parsing bedtools output
bedtools_parse = function(input){
        sta = character(0)
        sto = character(0)
        temp = character(0)
        exo = character(0)
        for (i in 1:nrow(input)){
                if (substring(input[i,], 1,nchar(">chr")) == ">chr"){
                        pars = unlist(str_split(input[i,], pattern = ":"))[2]
                        pars = unlist(str_split(pars, pattern = "-"))
                        exo = c(exo, temp)
                        sta = c(sta, pars[1])
                        sto = c(sto, pars[2])
                        temp = character(0)
                }else{
                        temp = c(temp, as.character(input[i,1]))
                        temp = paste(temp, collapse = '')
                }
                if (i == nrow(input)){
                        exo = c(exo, temp) 
                }
        }
        return(list(sta, sto, exo))
}

## bedtools parsing
output = read.table("output.txt", sep="\n")
res = bedtools_parse(output)
starts = unlist(res[1])
stops = unlist(res[2])
exons = unlist(res[3])

pos = as.numeric(unlist(str_split(data$phenotype_id[j],pattern = "_"))[2])
pos

## needed exon extraction
print(pos)
print(starts[pos])
print(stops[pos])

write(paste(starts[pos], stops[pos], sep = "-"), "exon_meta.txt")
needed_exon = exons[pos]
print(needed_exon)
needed_exon = paste0(">", unlist(str_split(data$phenotype_id[j],pattern = "_"))[1], "\n", needed_exon)
write(needed_exon, "exon_seq.fa")
