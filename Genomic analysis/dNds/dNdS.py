"""
Prepare dNdS (dndscv) input tables from MAF files.

This script converts MAF-like mutation tables into the simplified tabular format expected
by `dndscv` in downstream R scripts.

Two modes are supported:
- Patient-level: convert `patientLevel/<group>Clonality.maf` into `<group>.tsv`
- Region-level: split a region-level MAF by TLS grouping and write multiple TSVs

Output columns (dndscv convention)
- sampleID : sample identifier
- chr      : chromosome
- pos      : 1-based genomic position
- ref      : reference allele
- mut      : alternate allele
- clonality: project-specific clonality label
"""

from util import *

doMultiProcess = RunMultiProcess()

def MafOutcome1():
    """
    Patient-level export:
    Convert MAFs to `dndscv` TSV format for short/long/allSample cohorts.
    """
    os.chdir('/home_2/wzt/analysisOutcome/dNdS/PatientLevel')
    for i in ['short', 'long', 'allSample']:
        filein = '/home_2/wzt/analysisOutcome/maftools/patientLevel/{}Clonality.maf'.format(i)
        alldat = pd.read_csv(filein, sep='\t')
        # Keep only the columns needed for dndscv and rename them to its expected schema.
        alldat = alldat[['Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'clonality']]
        alldat.columns = ['sampleID', 'chr', 'pos', 'ref', 'mut', 'clonality']
        alldat.to_csv('{}.tsv'.format(i), sep='\t', header=True, index=False)

        
        
def MafOutcome2():
    """
    Region-level export:
    Subset a region-level MAF into TLS-defined groups and write one TSV per group.
    """
    os.chdir('/home_2/wzt/analysisOutcome/dNdS/RegionLevel')
    filein = '/home_2/wzt/analysisOutcome/maftools/regionLevel/allSampleClonality.maf'.format(i)
    clinical = pd.read_csv('/home_2/wzt/analysisOutcome/maftools/regionLevel/clinical.tsv', sep='\t')
    for i, j in zip(['allSample', 'MatureTLS', 'NoTLS', 'ImmatureTLS'], ["allSample", "Mature TLS", "No TLS", "Immature TLS"]):
        alldat = pd.read_csv(filein, sep='\t')
        if i != 'allSample':
            # Use the clinical table to select samples belonging to the target TLS category.
            Tumor_Sample_Barcode = clinical[clinical['group'] == j]['Tumor_Sample_Barcode']
            alldat = alldat[alldat['Tumor_Sample_Barcode'].isin(Tumor_Sample_Barcode)]   
        alldat = alldat[['Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'clonality']]
        alldat.columns = ['sampleID', 'chr', 'pos', 'ref', 'mut', 'clonality']
        alldat.to_csv('{}.tsv'.format(i), sep='\t', header=True, index=False)

if __name__ == '__main__':
    #MafOutcome1()
    MafOutcome2()
    
