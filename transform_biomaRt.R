####################################################################################################
##                                               LIBRARIES
####################################################################################################

library(biomaRt)

####################################################################################################
##                                                 DATA
####################################################################################################

load("~/Documents/Annotations/Mouse2HumanTable")

list_files = list.files(".",pattern = ".gct", full.names = F)

for(file in list_files) {
  print(file)
  data = read.table(file,h =T, sep='\t', stringsAsFactors = FALSE)
  i = 1
  for (gene in data$Name) {
    gene_name = unlist(strsplit(gene,".",fixed = T))[1]
    pos = which(gene_name == Mouse2HumanTable$Mouse.Gene_ID)
    if (length(pos) == 0 ) {
      data[i,"human"] = NA
    }
    else {
      data[i,"human"] = Mouse2HumanTable[pos[1],"Human.Gene_ID"]
    }
    i = i+1
  }

  data$Name = data$human
  data$human = NULL
  write.table(data,paste0("human_",file),col.names = T, row.names = F,quote = F,sep = '\t')
}
