library(data.table)
library(stringr)


seq  = as.character(fread("exon_seq.fa"))
coord = readLines("exon_meta.txt")
name = str_remove(readLines("exon_seq.fa")[1], ">")
print(name)
if (file.exists("data/top_sQTLs_with_top_coloc_with_AF.csv")){
        input = read.csv("data/top_sQTLs_with_top_coloc_with_AF.csv")
	colnames(input) = c("NC SEQ","COORDS", "GENE ID", "PROT REFSEQ", "START", "END", "SEQ", "ALIGN LENGTH", "EVAL", "SCORE", "ASN #", "CYS #")
} else {
        input = data.frame()
}
line = tryCatch(
	expr = {
		data = fread("try_blastx.txt", header=F)
		data.table(seq, coord, data[1,], str_count(data[1, "V5"], pattern="C"), str_count(data[1, "V5"], pattern="N"))
	},
	error = function(e){ 
		data.table(seq, coord, name, NA, NA, NA, NA, NA ,NA, NA, NA, NA)
	}
)
print(line)
colnames(line) = c("NC SEQ", "COORDS", "GENE ID", "PROT REFSEQ", "START", "END", "SEQ", "ALIGN LENGTH", "EVAL", "SCORE", "ASN #", "CYS #")


line = rbind(input, line)
write.csv(line, "data/top_sQTLs_with_top_coloc_with_AF.csv", row.names = F)
