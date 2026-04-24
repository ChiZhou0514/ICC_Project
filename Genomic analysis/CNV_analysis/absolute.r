## ABSOLUTE copy-number inference from CNV segmentation.
##
## This script runs the ABSOLUTE R package on a precomputed segment file produced from
## CNVkit outputs (e.g., `cnv/tumor2absolute.seg`).
##
## Inputs (relative to the working directory set by the caller script)
## - `cnv/tumor2absolute.seg` : segment means (log2 ratios) in ABSOLUTE format
## - `gatk/4pyclone_snpindels/output_pass2absolute.maf` : simplified MAF-like table
##
## Outputs
## - ABSOLUTE results under `results.dir = 'absolute'`
##
## Notes
## - Parameter ranges (ploidy/sigma) are project-specific; tune them if ABSOLUTE fails.
library(ABSOLUTE)
RunAbsolute('cnv/tumor2absolute.seg', min.ploidy = 0.5, max.ploidy = 10, 
            max.sigma.h = 0.2, platform = 'Illumina_WES', 
            copy_num_type = 'total', sigma.p = 0, results.dir = 'absolute', 
            primary.disease = 'CHOL', sample.name = 'tumor', 
            max.as.seg.count = 4000, max.non.clonal = 0.05, 
            max.neg.genome = 0.005, 
            maf.fn = 'gatk/4pyclone_snpindels/output_pass2absolute.maf', 
            min.mut.af = 0.05
)
