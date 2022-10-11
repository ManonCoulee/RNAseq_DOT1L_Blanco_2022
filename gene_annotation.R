####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('edgeR', 'DESeq2', 'FactoMineR', 'reshape', 'TCseq', 'statmod','GenomicFeatures')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

wd = '.'
sample_dir = "~/Documents/Mouse/RNAseq/Samples"
dirs = list.dirs(sample_dir, full.names = T, recursive = F)
files = list.files(dirs,"RPG")

rout = wd

LFC = log2(1.5)

CT = list.dirs(sample_dir, full.names = F, recursive = F)

## Load SamplePlan data descriptor
SamplePlan = read.table(paste0(wd,"/SamplePlan.tsv"),sep = "\t", h = T, row.names = 1)
SamplePlan$SampleType = factor(SamplePlan$SampleType)
SamplePlan$CellType = factor(SamplePlan$CellType,levels = CT)
SamplePlan$SamplePool = factor(SamplePlan$SamplePool)
SamplePlan$CellTypeRed = factor(SamplePlan$CellTypeRed)

####################################################################################################
##                                                  ANNOTATION
####################################################################################################

## Gencode gene counts from STAR
x = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)), function(path){
  read.table(paste0(path, '.RPG.tsv'), sep = "\t", header = F, skip = 4)[,4]
}, simplify = T),
  row.names = read.table(paste0(SamplePlan[1,"SamplePath"],"/",rownames(SamplePlan)[1],".RPG.tsv"),
  sep = "\t", header = F, skip = 4)[,1])

## Gene counts data normalization
min_sample = 2
min_expression = 1
xx = cpm(y = x, normalized.lib.sizes = T, log = F)
expressed_genes = rowSums(xx > min_sample) >= min_expression ## Filter-out low expressed genes

genecounts = list(
  ## Raw read counts (as estimated by STAR)
  raw = x[expressed_genes, ],
  ## For expression data visualization
  cpm = cpm(y = x[expressed_genes, ], lib.size = colSums(x), normalized.lib.sizes = T, log = F),
  ## For performing PCA analysis
  rlog = assay(rlog(DESeqDataSetFromMatrix(countData = x,
    colData = SamplePlan, design = ~1)))[expressed_genes, ]
)

write.table(genecounts$raw,"genecounts_SC_SCII_RS_raw",row.names = T, col.names = T, quote = F,
  sep = "\t")

## Load gene annotations (from Genecode transcripts.fa)
gene_annotation = data.frame(do.call(rbind, strsplit(x = gsub(
  pattern = "^>", replacement = '', x = unlist(system(paste('zgrep -P "^>"',
    file.path('~/Documents/Annotations/gencode.vM19.transcripts.fa.gz')), intern = T)),
  perl = T), split = '|', fixed = T)), stringsAsFactors = F)
names(gene_annotation) = c('tx_ENS', 'gene_ENS', 'gene_OTT', 'tx_OTT', 'tx_name', 'gene_name',
  'tx_length', 'gene_type')
gene_annotation[gene_annotation == '-'] = NA
gene_annotation$tx_length = as.numeric(gene_annotation$tx_length)

## Estimating gene_length by using GenomicFeatures
txdb = makeTxDbFromGFF(file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'), format = 'gtf')
exons.list.per.gene = exonsBy(x = txdb, by = 'gene')
## Same but much faster
exonic.gene.sizes = sum(width(reduce(exons.list.per.gene)))
genecounts$rpkm = rpkm(y = x, gene.length = exonic.gene.sizes[match(names(exonic.gene.sizes),
  rownames(x))], normalized.lib.sizes = T, log = F)[expressed_genes, ]

table = data.frame(genecounts$cpm)
colnames(table) = rownames(SamplePlan)
table$gene_ENS = rownames(genecounts)

table_CTL = table[,grepl("Ctl",colnames(table))]

for(dir in CT) {
  x = table_CTL[,grepl(paste0(dir,"cell"),colnames(table_CTL))]
  x$gene_expression = (x[,1] + x[,2])/2
  x$gene_ENS = rownames(x)

  ## Create file usable for GSEA
  xx = table[,grepl(paste0(dir,"cell"),colnames(table))]
  write.table(xx,paste0(dir,"_genecounts.gct"),row.names = F, col.names = T, quote = F,
    sep = "\t")
}
