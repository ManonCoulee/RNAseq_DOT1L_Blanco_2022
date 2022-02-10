####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('UpSetR','ggplot2')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

dir_romain_data = "~/Documents/Mouse/RNAseq/RNA_seq_Romain"
rout = file.path(".")
samples = c("RS","SC","SCII")

####################################################################################################
##                                                 DATA
####################################################################################################

list_files = list.files(dir_romain_data),
  pattern = "cells_volcano_table.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})

####################################################################################################
##                                      UPSET GENES DEREGULATED
####################################################################################################

gene = unique(c(list_data$SC$gene_ENS, list_data$RS$gene_ENS, list_data$SCII$gene_ENS))

RS = list_data$RS[list_data$RS$PValue <= 0.05,]
SC = list_data$SC[list_data$SC$PValue <= 0.05,]
SCII = list_data$SCII[list_data$SCII$PValue <= 0.05,]

RS_down = RS[RS$gl == "Down-regulated","gene_ENS"]
RS_up = RS[RS$gl == "Up-regulated","gene_ENS"]
SC_down = SC[SC$gl == "Down-regulated","gene_ENS"]
SC_up = SC[SC$gl == "Up-regulated","gene_ENS"]
SCII_down = SCII[SCII$gl == "Down-regulated","gene_ENS"]
SCII_up = SCII[SCII$gl == "Up-regulated","gene_ENS"]

up_down_table = data.frame("gene" = gene, "SC_up" = 0, "SC_down" = 0,
  "SCII_up" = 0, "SCII_down" = 0, "RS_up" = 0, "RS_down" = 0)

up_down_list = list(SC_up,SC_down,SCII_up,SCII_down,RS_up,RS_down)
names(up_down_list) = colnames(up_down_table)[-1]

for(type in colnames(up_down_table)[-1]){
  print(type)
  pos = unlist(lapply(up_down_list[[type]], function(x) {
    which(up_down_table$gene == x)
  }))
  up_down_table[pos,type] = 1
}

png("SC_SCII_RS_upset.png",width = 1500,height = 800)
upset(up_down_table, sets = c("SC_up","SC_down", "SCII_up", "SCII_down", "RS_up", "RS_down"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  sets.bar.color = rep(c("red","blue"), 3), text.scale = c(3,3,3,3,4,5),
  number.angles = 0,
  point.size = 6, line.size = 2)
dev.off()

RS_deg = RS[RS$gl != "Not regulated","gene_ENS"]
SC_deg = SC[SC$gl != "Not regulated","gene_ENS"]
SCII_deg = SCII[SCII$gl != "Not regulated","gene_ENS"]

deg_table = data.frame("gene" = gene,"SC" = 0, "SCII" = 0, "RS" = 0)

deg_list = list(SC_deg,SCII_deg,RS_deg)
names(deg_list) = colnames(deg_table)[-1]

for(type in colnames(deg_table)[-1]){
  print(type)
  pos = unlist(lapply(deg_list[[type]], function(x) {
    which(deg_table$gene == x)
  }))
  deg_table[pos,type] = 1
}

png("SC_SCII_RS_upset_reduce.png",width = 1080,height = 1080)
upset(deg_table, sets = c("SC_deg","SCII_deg","RS_deg"),
  order.by = "freq", keep.order = T, mainbar.y.label = "number of genes",
  sets.x.label = "Number gene\nper celltype",
  text.scale = c(3,3,3,3,4,5),
  number.angles = 0,
  point.size = 6, line.size = 2)
dev.off()
