library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
library('NMF')
library('pheatmap')
## maftools-based mutation landscape analysis (TLS grouping, region level).
##
## This script loads a region-level MAF with clonality annotations and compares
## mutation patterns across TLS categories (e.g., Mature TLS vs others).
##
## Inputs (working directory)
## - `clinical.tsv`            : sample -> TLS group annotation
## - `cnTable.tsv`             : CN status table for maftools (optional)
## - `allSampleClonality.maf`  : merged MAF with a `clonality` column
##
## Common outputs
## - Ti/Tv summaries and oncoplots
## - Statistical comparisons of selected mutation classes
driver = c('KRAS', 'IDH1', 'NRAS', 'BAP1', 'TP53', 'ARID1A', 'PBRM1',
           'STK11', 'EPHA2', 'SMARCA4', 'VARS', 'OBSL1', 'TGFBR1')

### TLS group
rmFlags = TRUE
setwd("/home_2//analysisOutcome/maftools/regionLevel/")
clinicalData = read.table('clinical.tsv', sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
cnTable = read.table('cnTable.tsv', sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
allSample = read.maf(maf = 'allSampleClonality.maf', cnTable = cnTable, clinicalData = 'clinical1.tsv', rmFlags = rmFlags, )

allSample = subsetMaf(maf = allSample, query = "clonality == 'clone'")
CHOL.titv = titv(maf = allSample, plot = FALSE, useSyn = TRUE)
plotTiTv(res = CHOL.titv, sampleOrder = as.character(clinicalData$Tumor_Sample_Barcode), showBarcodes = TRUE)

a = clinicalData[clinicalData$group == 'short', 'Tumor_Sample_Barcode']
b = CHOL.titv$fraction.contribution
t1 = b[!b$Tumor_Sample_Barcode %in% a, ]$`C>T`
t2 = b[b$Tumor_Sample_Barcode %in% a, ]$`C>T`
t.test(t1, t2)



#longSur = subsetMaf(maf = longSur, query = "(IMPACT1 == 'HIGH') | (IMPACT1 == 'HIGH')")
#longSur = subsetMaf(maf = longSur, query = "IMPACT1 == 'HIGH'")
#shortSur = subsetMaf(maf = shortSur, query = "IMPACT1 == 'HIGH'")

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
         annotationOrder=c('Mature TLS', 'Immature TLS', 'No TLS','NA')
        )

oncoplot(maf = allSample, top = 20)
oncoplot(maf = allSample, genes = driver)


plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE)
plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE, genes = genes)
plotVaf(maf = allSample, vafCol = 'VAF', flip = FALSE, genes = driver)
OncogenicPathways(maf = shortSur)
PlotOncogenicPathways(maf = shortSur, pathways = "RTK-RAS")
PlotOncogenicPathways(maf = shortSur, pathways = "Cell_Cycle")
PlotOncogenicPathways(maf = shortSur, pathways = "TP53")
OncogenicPathways(maf = longSur)

pdf('coOncoplot.pdf', width=8, height=4)




