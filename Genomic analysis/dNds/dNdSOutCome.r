library('dndscv')
## dN/dS analysis by outcome grouping using dndscv.
##
## This script runs `dndscv` on patient-level TSV files (prepared by `dNdS.py`) and writes:
## - global dN/dS summaries
## - gene-level significant results (`sel_cv`) filtered by q-value threshold
##
## Groups
## - long survival (long.tsv)
## - short survival (short.tsv)
## Each group is also split by clonality labels (clone vs subclone).
setwd("/home_2/wzt/analysisOutcome/dNdS/")
refdb = '/home/wzt/database/dNdS/RefCDS_human_GRCh38.p12.rda'

### outcome grouping


### longSur group
datlong = read.table('long.tsv', sep = '\t', header = 1)
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/longoutGene.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'longout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('long.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'clone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/longCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'longCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('long.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'subclone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/longSubCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'longSubCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)

### shortSur group
datlong = read.table('short.tsv', sep = '\t', header = 1)
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/shortoutGene.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'shortout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('short.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'clone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/shortCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'shortCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('short.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'subclone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'GeneLevel/shortSubCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'shortSubCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)

### all Sample analysis
datlong = read.table('allSample.tsv', sep = '\t', header = 1)
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'allSampleoutGene.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'allSampleout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('allSample.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'clone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'allSampleCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'allSampleCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)
datlong = read.table('allSample.tsv', sep = '\t', header = 1)
datlong = datlong[datlong$clonality == 'subclone', ]
dndsout_long = dndscv(datlong, refdb = refdb, cv = NULL)
sel_cv_long = dndsout_long$sel_cv
signif_genes_long = sel_cv_long[sel_cv_long$qglobal_cv<0.1, ]
write.table(signif_genes_long, file = 'allSampleSubCloneGeneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(dndsout_long$globaldnds, file = 'allSampleSubCloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsout_long$globaldnds)

