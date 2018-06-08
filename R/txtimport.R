##########################################
# Salmon matrix generation by David Brandao
#Z:/ is a generic name for a local server

#load packages
library(tximport)
library(readr)

# Read a single file of patient's expression to get colnames
FILES<-list.files("Z:/quants/")
matrix_expression = read.table(paste("Z:/quants/",FILES[1],"/quant.sf",sep = ""),header = T)
head(matrix_expression)

#---------------------from transcripts to genes---------------------------#
transcript_name = matrix_expression$Name
transcript_name.list = lapply(transcript_name,function(x){unlist(strsplit(as.character(x),split="[|]"))[6]}) 
transcript_name.list = unlist(transcript_name.list)
t2xgene = data.frame(transcript_name,transcript_name.list)
colnames(t2xgene) = c("target_id","GENEID")
matrix_expression$GENEID<-transcript_name.list

# do the matrix
dir_in=list.files("Z:/quants/")
dir_in.path = file.path("Z:/quants",dir_in,"quant.sf")
dir_in.path<-dir_in.path[file.exists(dir_in.path)]
length(unique(dir_in))

# run it
SalmonGENEXP <- tximport(file= dir_in.path, type="salmon", tx2gene=t2xgene,dropInfReps=TRUE ) # reader=read_tsv

colnames(SalmonGENEXP$abundance)<-dir_in
colnames(SalmonGENEXP$counts)<-dir_in
colnames(SalmonGENEXP$length)<-dir_in

# save it
save(SalmonGENEXP,file = "SalmonGENEXP.RData")
