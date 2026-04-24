library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
library('NMF')
library('pheatmap')
## maftools-based mutation landscape analysis (region level).
##
## This script loads region-level MAF and clinical annotations and produces common
## maftools visualizations:
## - Ti/Tv summary
## - Oncoplots for short vs long survival groups
## - Pathway and VAF plots
##
## Inputs (working directory)
## - `clinical1.tsv` : sample -> outcome group annotation (short/long)
## - `cnTable.tsv`   : gene-level CN status table (optional)
## - `allSample.maf` : merged MAF for all regions
##
## Notes
## - The script subsets to clonal mutations (`clonality == 'clone'`) and high-impact calls
##   using `IMPACT1` (project-specific field).
driver = c('KRAS', 'IDH1', 'NRAS', 'BAP1', 'TP53', 'ARID1A', 'PBRM1',
           'STK11', 'EPHA2', 'SMARCA4', 'VARS', 'OBSL1', 'TGFBR1')

### outcome  group
rmFlags = TRUE
setwd("/home_2//analysisOutcome/maftools/regionLevel/")

clinicalData = read.table('clinical1.tsv', sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
cnTable = read.table('cnTable.tsv', sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
allSample = read.maf(maf = 'allSample.maf', cnTable = cnTable, clinicalData = 'clinical1.tsv', rmFlags = rmFlags, )

allSample = subsetMaf(maf = allSample, query = "clonality == 'clone'")
CHOL.titv = titv(maf = allSample, plot = FALSE, useSyn = TRUE)
plotTiTv(res = CHOL.titv, sampleOrder = as.character(clinicalData$Tumor_Sample_Barcode), showBarcodes = TRUE)

a = clinicalData[clinicalData$group == 'short', 'Tumor_Sample_Barcode']
b = CHOL.titv$fraction.contribution
t1 = b[!b$Tumor_Sample_Barcode %in% a, ]$`C>T`  
t2 = b[b$Tumor_Sample_Barcode %in% a, ]$`C>T` 
t.test(t1, t2)

tmp = allSample@data[allSample@data$Hugo_Symbol == 'MUC4']$IMPACT1
tmp[is.na(tmp)] = 'HIGH'
allSample@data[allSample@data$Hugo_Symbol == 'MUC4']$IMPACT1 = tmp
tmp = allSample@data[allSample@data$Hugo_Symbol == 'CDKN2A']$IMPACT1
tmp[is.na(tmp)] = 'HIGH'
allSample@data[allSample@data$Hugo_Symbol == 'CDKN2A']$IMPACT1 = tmp
allSample = subsetMaf(maf = allSample, query = "IMPACT1 == 'HIGH'")

oncoplot(maf = shortSur, top = 10, removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE)
oncoplot(maf = longSur, top = 10, removeNonMutated = FALSE, showTumorSampleBarcodes = TRUE)
st.vs.lg <- mafCompare(minMut = 5, m1 = shortSur, m2 = longSur, m1Name = 'shortSur', m2Name = 'longSur')
print (st.vs.lg)
genes = c('TP53', 'KRAS')
genes = c('TP53', 'KRAS', 'CDKN2A')
coOncoplot(showSampleNames = TRUE, genes = genes, m1 = shortSur, m2 = longSur, m1Name = 'shortSur', m2Name = 'longSur', removeNonMutated = FALSE)
coOncoplot(genes = driver, m1 = shortSur, m2 = longSur, m1Name = 'shortSur', m2Name = 'longSur', removeNonMutated = FALSE)
lollipopPlot2(gene = "TP53", m1 = longSur, m2 = shortSur, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m2_name = "shortSur", m1_name = "longSur")



plotmafSummary(maf = allSample, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = TRUE, )
oncoplot(maf = allSample, clinicalFeatures = 'group', removeNonMutated = FALSE,
         sortByAnnotation = TRUE, drawColBar = TRUE, altered = TRUE, 
         drawRowBar = TRUE, showTumorSampleBarcodes = TRUE, top=15, 
         #genes = c('TP53', 'KRAS', 'CDKN2A', 'MUC4')
         )

oncoplot(maf = allSample, top = 20,removeNonMutated = FALSE)
oncoplot(maf = allSample, genes = driver, removeNonMutated = FALSE)


plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE) 
plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE, genes = genes)
plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE, genes = driver)
OncogenicPathways(maf = shortSur)
PlotOncogenicPathways(maf = shortSur, pathways = "RTK-RAS")
PlotOncogenicPathways(maf = shortSur, pathways = "Cell_Cycle")
PlotOncogenicPathways(maf = shortSur, pathways = "TP53")
OncogenicPathways(maf = longSur)

pdf('coOncoplot.pdf', width=8, height=4)

filein1 = '/home_2//analysisOutcome/gistic2/all/all_lesions.conf_99.txt'
filein2 = '/home_2//analysisOutcome/gistic2/all/amp_genes.conf_99.txt'
filein3 = '/home_2//analysisOutcome/gistic2/all/del_genes.conf_99.txt'
filein4 = '/home_2//analysisOutcome/gistic2/all/scores.gistic'
tumor = readGistic(filein1, filein2, filein3, filein4)
qvalue = 1e-5
gisticBubblePlot(tumor, fdrCutOff = qvalue, log_y = TRUE, markBands = '9p21.3')
gisticBubblePlot(tumor, fdrCutOff = qvalue, log_y = TRUE, markBands = tumor@cytoband.summary[tumor@cytoband.summary$qvalues <= qvalues]$Cytoband[1:10])
gisticBubblePlot(tumor, fdrCutOff = qvalue, log_y = TRUE, markBands = tumor@cytoband.summary[tumor@cytoband.summary$qvalues <= qvalues]$Cytoband)
gisticChromPlot(gistic = tumor, fdrCutOff = qvalue, ref.build = 'hg38', markBands = "all")

