"""
Cross-cohort mutation frequency comparison.

This script collects MAF files from multiple cohorts and compares mutation frequencies
for a curated gene list across:
- Our cohort (all / short survival / long survival)
- CPTAC iCCA cohort
- TCGA cohort
- Zou cohort (project-specific external cohort)

Key utilities
- `TCGA_vcfProcess`: convert raw TCGA VCFs to MAF (PASS + AF/depth filter, liftOver, vcf2maf)
- `TCGA_merge` / `NC_merge`: merge per-sample MAFs into cohort-level MAFs
- `calFre`: compute per-gene mutation frequencies and export a bar plot

Notes
- Paths are hard-coded and assume an HPC environment.
- Frequencies use fixed denominators (e.g., 58, 32, 26, 253, 44, 103) matching cohort sizes.
"""

from itertools import chain
from util import *
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import stats

from util import getOutcome

def TCGA_vcfProcess():
    """
    Convert TCGA tumor/normal VCFs into filtered MAF files per patient.

    Output
    - `lifted_overPASS.maf` written in each `TCGA*/` directory.
    """
    dirNames = glob.glob('/home_2//ICC_TCGA/SNVIndel/TCGA*')
    for dirName in tqdm(dirNames):
        os.chdir(dirName)
        try:
            patient = os.path.basename(dirName)
            if os.path.isfile('lifted_over.vep.vcf'): os.remove('lifted_over.vep.vcf')
            # Keep PASS records and apply minimal AF/AD filters to remove low-confidence calls.
            cmd1 = "bcftools view  -s  NORMAL,TUMOR  -o output_filterPro.vcf  -M2 -O v -i '(FILTER=\"PASS\")' output_filter.vcf"
            cmd2 = "bcftools filter  -i '(FMT/AF[1:0] >= 0.05 && SUM(FMT/AD[1:])>= 25  &&  FMT/AD[1:1]>=5)'   output_filterPro.vcf   -o output_pass.vcf   -O v"
            cmd3 = 'java -jar /home//software/picard.jar  LiftoverVcf  I=output_pass.vcf   O=lifted_over.vcf   CHAIN=//home//database/UCSC/hg38ToHg19.over.chain   REJECT=rejected_variants.vcf   R=/home//database/NCBI/hg19/hg19.fa  WARN_ON_MISSING_CONTIG=true'
            cmd4 = 'perl /home//software/vcf2maf-1.6.19/vcf2maf.pl  --input-vcf  lifted_over.vcf    --output-maf   lifted_over.maf    --vep-path /home/zhouchi/software/ensembl-vep/  --vep-data ' \
            '/home//database/ensembl/hg19/vep-cache/   --ref-fasta   /home//database/NCBI/hg19/hg19.fa  --tumor-id TUMOR  --normal-id NORMAL'
            with open('log.txt', 'a') as fout:
                p = subprocess.Popen('&&'.join([cmd1, cmd2, cmd3, cmd4]), shell=True, stderr=subprocess.STDOUT, stdout=fout)
                p.wait()
            dat = pd.read_csv('lifted_over.maf', sep='\t', skiprows=1)
            dat['VAF'] = round(dat['t_alt_count'] / dat['t_depth'], 4)
            columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Variant_Classification',
            'Strand', 'Variant_Type', 'HGVSp_Short', 'Consequence', 'VAF', 'dbSNP_RS', 'Existing_variation', 'BIOTYPE',
            'SIFT', 'PolyPhen', 'IMPACT']
            dat = dat[columns]
            dat['Tumor_Sample_Barcode'] = patient
            dat.to_csv('lifted_overPASS.maf', sep='\t', index=False)
        except Exception as e:
            print (e); print (dirName)


### merge MAF files
def TCGA_merge():
    """Merge all TCGA per-sample MAFs into `TCGA.maf`."""
    os.chdir('/home_2//analysisOutcome/cohortCompare')
    alldat = []
    mafs = glob.glob('/home_2//ICC_TCGA/SNVIndel/TCGA*/lifted_overPASS.maf')
    for i in mafs:
        dat = pd.read_csv(i, sep='\t')
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('TCGA.maf', sep='\t', header=True, index=False)



def NC_merge():
    """Merge all external 'NC' cohort per-sample MAFs into `NC.maf`."""
    os.chdir('/home_2//analysisOutcome/cohortCompare')
    alldat = []
    mafs = glob.glob('/home_2//ICC_103/p*/WES/T/gatk/lifted_overPASS.maf')
    for i in mafs:
        dat = pd.read_csv(i, sep='\t')
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('NC.maf', sep='\t', header=True, index=False)

def calFre():
    """
    Compute and plot mutation frequencies for the selected gene list across cohorts.

    Outputs
    - `Gene_fre.tsv` : frequency table (genes x cohorts)
    - `Gene_fre.pdf` : bar plot
    """
    os.chdir('/NFS_home/NFS_home_3//analysisOutcome/cohortCompare')
    dat = pd.DataFrame(index=cohorts, columns=genes)
    shortSur, longSur = getOutcome()


    ourCohort = pd.read_csv('ourCohort_patient.maf', sep='\t')
    ourCohort = ourCohort[ourCohort['Variant_Classification'].isin(NonSyn)]

    shortSurCohort = ourCohort[ourCohort['Tumor_Sample_Barcode'].isin(shortSur)]
    longSurCohort = ourCohort[ourCohort['Tumor_Sample_Barcode'].isin(longSur)]


    CPTAC = pd.read_csv('253_samples_CPTAC_iCCA.txt', sep='\t')
    CPTAC = CPTAC[CPTAC['Mutation_Type'] != 'synonymous_variant']

    Zou = pd.read_csv('NC.maf', sep='\t')
    Zou = Zou[Zou['Variant_Classification'].isin(NonSyn)]
    
    TCGA = pd.read_csv('TCGA.maf', sep='\t')
    TCGA = TCGA[TCGA['Variant_Classification'].isin(NonSyn)]



    for gene in genes:
        frq1 = round(len(ourCohort[ourCohort['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique()) / 58, 3)
        frq2 = round(len(shortSurCohort[shortSurCohort['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique()) / 32, 3)
        frq3 = round(len(longSurCohort[longSurCohort['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique()) / 26, 3)
        
        frq4 = round(len(CPTAC[CPTAC['Gene'] == gene]['Sample_ID'].unique()) / 253, 3)
        frq5 = round(len(TCGA[TCGA['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique()) / 44, 3)
        frq6 = round(len(Zou[Zou['Hugo_Symbol'] == gene]['Tumor_Sample_Barcode'].unique()) / 103, 3)
        dat[gene] = [frq1, frq2, frq3, frq4, frq5, frq6]



    dat.to_csv('Gene_fre.tsv', sep='\t')

    dat = dat.T
    fig, axes = plt.subplots(figsize=(10, 7))
    dat.plot.bar(ax = axes)
    fig.savefig("Gene_fre.pdf", transparent=True, bbox_inches='tight')
    
        
        
genes = ['TP53', 'KRAS', 'ARID1A', 'ARID1B', 'BAP1', 'PBRM1', 'IDH1', 'FGFR2',  'LAMA3', 'SPEN',   'FRG1', 'KMT2A',  
 'NRAS', 'STK11', 'PIK3CA', 'BRCA1', 'BRCA2', 'ATM',  'ADAP1']

NonSyn = ['Missense_Mutation', 'Splice_Region', 'Splice_Site', 'Nonsense_Mutation',
'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']

cohorts = ['ourCohort', 'shortSurCohort', 'longSurCohort', 'CPTAC', 'TCGA', 'Zou']



if __name__ == '__main__':
    print ('hello, world')
    calFre()
