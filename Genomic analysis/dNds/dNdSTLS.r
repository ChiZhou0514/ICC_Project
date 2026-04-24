library('dndscv')
## dN/dS analysis by TLS grouping using dndscv.
##
## This script runs `dndscv` on region-level TSV files (prepared by `dNdS.py`) and writes
## global dN/dS summaries for:
## - No TLS
## - Immature TLS
## Each group is additionally split by clonality labels (clone vs subclone).
##
## Inputs (working directory)
## - `NoTLS.tsv`, `ImmatureTLS.tsv`, ... (dndscv-style mutation tables)
##
## Outputs
## - `*_out.tsv`, `*_Cloneout.tsv`, `*_Subcloneout.tsv` global dN/dS summaries
setwd("/home_2/wzt/analysisOutcome/dNdS/RegionLevel/")
refdb = '/home/wzt/database/dNdS/RefCDS_human_GRCh38.p12.rda'


### TLS grouping

### NoTLS
datNoTLS = read.table('NoTLS.tsv', sep = '\t', header = 1)
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'NoTLS_out.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('NoTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'subclone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'NoTLS_Subcloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('NoTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'clone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'NoTLS_Cloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)


## Immature
datNoTLS = read.table('ImmatureTLS.tsv', sep = '\t', header = 1)
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'ImmatureTLS_out.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('ImmatureTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'subclone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'ImmatureTLS_Subcloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('ImmatureTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'clone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'ImmatureTLS_Cloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

### mature
datNoTLS = read.table('MatureTLS.tsv', sep = '\t', header = 1)
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'MatureTLS_out.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('MatureTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'subclone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'MatureTLS_Subcloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)

datNoTLS = read.table('MatureTLS.tsv', sep = '\t', header = 1)
datNoTLS = datNoTLS[datNoTLS$clonality == 'clone', ]
dndsOut_NoTLS = dndscv(datNoTLS, refdb = refdb, cv = NULL)
write.table(dndsOut_NoTLS$globaldnds, file = 'MatureTLS_Cloneout.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
print(dndsOut_NoTLS$globaldnds)





