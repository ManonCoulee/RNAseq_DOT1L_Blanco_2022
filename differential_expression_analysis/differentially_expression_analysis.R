####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'FactoMineR', 'reshape', 'TCseq', 'statmod')
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
expressed_genes = rowSums(xx > min_expression) >= min_sample ## Filter-out low expressed genes

genecounts = read.table("genecounts_filtered_raw.tsv", sep = "\t", h = T, stringsAsFactors = F)

####################################################################################################
##                                              DE ANALYSIS
####################################################################################################

## Perform differential gene expression analysis
cell = SamplePlan$CellType
phenotype = relevel(SamplePlan$SampleType, ref = "CTL")
		 
## Create a model with all multifactor combinations
design = model.matrix(~phenotype+cell+phenotype:cell)
colnames(design) = c("all","KOvsCTL.RS","CTL.SC","CTL.SCII","KOvsCTL.SC","KOvsCTL.SCII")

## Creating DGEList object
y = DGEList(counts = genecounts)
## Calculating TMM-based scaling factors
y = calcNormFactors(object = y)
## Estimating data dispersion
y = estimateDisp(y = y, design = design)

## Plot MDS
n_top = 500
o = plotMDS(x = y, plot = FALSE, top = n_top)
p = ggplot(data = data.frame(x = o$x, y = o$y), aes(x = x, y = y, shape = SamplePlan$SampleType,
    color = SamplePlan$CellType)) +
  geom_point(size = 10, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
  geom_hline(yintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
  geom_text(aes(x = x, y = y), label = SamplePlan$SamplePool, fontface = 'bold', size = 4,
    alpha = .5, hjust = .5, vjust = .5, show_guide = FALSE) +
  scale_color_manual(name = 'CellType', values = c(brewer.pal(9, 'Set1'))) +
  scale_shape_discrete(name = 'SampleType') +
  labs(title = 'Multi-dimensional scaling plot (MDS)',
    subtitle = paste0('Pairwize gene selection / cpm>',min_expression,' / n=',n_top,' top genes'),
    x = 'Leading LogFC dim. 1', y = 'Leading LogFC dim. 2') +
  theme(
    plot.title = element_text(size = 35, lineheight = 2, vjust = 1, face = 'bold'),
    plot.subtitle = element_text(size = 25, face = 'bold'),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    panel.border = element_rect(colour = 'grey50', fill = NA)
)
ggsave(file.path(rout, 'MDS.png'), plot = p, width = 12, height = 8, device = 'png', dpi = 300)

		 
## DEG analysis with glm model
fit = glmQLFit(y = y, design = design, robust = T)
		 		 
cmp = makeContrasts(
  'DOT1L KO effect within RS cells' = KOvsCTL.RS,
  'DOT1L KO effect within SC cells' = KOvsCTL.RS + KOvsCTL.SC,
  'DOT1L KO effect within SCII cells' = KOvsCTL.RS + KOvsCTL.SCII,
  levels = design
)
		 
thresholds = data.frame(
  FC = rep(1.5,3),
  P = rep(.05,3),
  row.names = colnames(cmp)
)

DEG = sapply(colnames(cmp), function(contrast){
  FC = log2(thresholds[contrast,'FC'])
  qlf = glmTreat(glmfit = fit, contrast = cmp[, contrast], lfc = LFC)
  tt = with(topTags(object = qlf, n = NULL, sort.by = 'none'), table)
  tt$gl = 0
  tt$gl[tt$PValue <= thresholds[contrast, 'P']  & tt$logFC < -FC] = 1
  tt$gl[tt$PValue <= thresholds[contrast, 'P']  & tt$logFC > FC] = 2
  tt$gl = factor(tt$gl, levels = c(1, 0, 2), labels = c('Down-regulated', 'Not regulated',
    'Up-regulated'))
  tt
}, simplify = F)

invisible(lapply(names(DEG), function(x){
  write.table(file.path(rout,paste0(gsub(pattern = ' ', replacement = '_', x = x),
    '_volcano_table.tsv')),
    x = data.frame(m_export,DEG[[x]], check.names = F), sep = "\t", row.names = F, quote = F)
}))

## DEG / MD plots
invisible(lapply(names(DEG), function(x){
  tt = DEG[[x]]

  p = ggplot(data = tt) +
    geom_point(aes(x = logCPM, y = logFC, color = gl), size = 5, alpha = .5) +
    geom_hline(yintercept = c(-log2(thresholds[x, 'FC']), log2(thresholds[x, 'FC'])), alpha = .5,
      color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = 0, alpha = .5, color = '#FFFFFF') +
    scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    labs(subtitle = paste(x, paste0('(*n=', sum(tt$gl != 'NR'), ')')),
      title = paste(paste0('FC>',thresholds[x,'FC']),paste0('& P<',100*thresholds[x,'P'],'%'),
        'MD plot'),
      x = 'Average LogCPM estimates',y = 'Log Fold-change') +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )+
    ylim(-6,6)

  ggsave(file.path(rout, paste0(gsub(pattern = ' ', replacement = '_', x = x), '_MD_plot.png')),
    plot = p, width = 10, height = 8, device = 'png', dpi = 300)
}))
