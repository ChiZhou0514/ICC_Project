"""
VCF post-processing and conversion to MAF.

This script wraps common steps for tumor-normal VCF processing:
- Apply allele frequency / read depth / allele depth filters with bcftools
- LiftOver VCF coordinates if needed (hg38 -> hg19 chain in this workflow)
- Convert VCF to MAF via vcf2maf + VEP
- Generate additional per-tumor VCF/MAF files for downstream clonal analysis (pyclone)

Directory assumptions
- Operates inside `.../WES/<tissue>/gatk/` folders.
- Expects `output_filter.vcf.gz` produced by the somatic caller.

Outputs
- `lifted_overPASS.maf` (hg19 MAF, filtered)
- `4pyclone_snpindels/` folder containing filtered VCF + indexes and MAF variants for pyclone.

Notes
- Commands are executed with `shell=True` and output appended to `log.txt`.
- This script is tuned for a specific environment; many paths are hard-coded.
"""

from util import *
from itertools import chain
import shutil, os

doMultiProcess = RunMultiProcess()


def f_subprocess(cmd):
    """
    Run a list of command segments joined by '&&' and append output to `log.txt`.
    """
    with open('log.txt', 'a') as fout:
        p = subprocess.Popen('&&'.join(cmd), shell=True, stderr=subprocess.STDOUT, stdout=fout)
        p.wait()
# conda  activate 2020plus
# conda  deactivate
# extract_gene_seq  -i  /home//database/NCBI/hg19/hg19.fa -b  snvboxGenes.bed  -o  snvboxGenes.fa

def fun2(line):
    """
    Collapse detailed MAF impact annotations into a simple `IMPACT1` label.

    The input `line` is a pandas Series or list-like row. This function implements
    a project-specific rule:
    - Keep `HIGH` as-is
    - Promote some `MODERATE` variants to `HIGH` if they are predicted deleterious
    - Otherwise return `LOW`
    """
    if line[20] == 'HIGH':
        return 'HIGH'
    elif line[20] == 'MODERATE':
        #tmp = ['probably_damaging', 'possibly_damaging']
        tmp = ['probably_damaging']
        if (line[18] is not np.nan and line[18].split('(')[0] == 'deleterious') or (line[19] is not np.nan and line[19].split('(')[0] in tmp): ### 
            return 'HIGH'
        else:
            return 'LOW'
    else:
        return 'LOW'

def vcfProcess(dirName):
    """
    Process one tumor tissue `gatk/` directory:
    1) Filter the VCF by PASS status and allele support
    2) LiftOver to hg19 (project convention) and convert to MAF
    3) Create a pyclone-friendly VCF subset (SNVs+indels) and index it

    Parameters
    ----------
    dirName : str
        Full path to a `.../gatk` directory.
    """
    try:
        os.chdir(dirName)
        caseid = dirName.split('/')[4]
        codeid = dirName.split('/')[6]
        if os.path.isfile('lifted_over.vep.vcf'): os.remove('lifted_over.vep.vcf')
        codeid_normal = 'NC'

        # 1) Keep only PASS records; keep both normal and tumor samples.
        cmd1 = "bcftools view  -s  {}{},{}{}  -o output_filterPro.vcf  -M2 -O v -i '(FILTER=\"PASS\")' output_filter.vcf.gz".format(caseid, codeid_normal, caseid, codeid)

        # 2) Filter by AF, total depth, and tumor alt depth thresholds.
        cmd2 = "bcftools filter  -i '(FMT/AF[1:0] >= 0.05  && SUM(FMT/AD[1:])>=50 &&  FMT/AD[1:1] >=7)'   output_filterPro.vcf   -o output_pass.vcf   -O v"

        # 3) LiftOver from hg38 -> hg19 to match downstream references for this project.
        cmd3 = 'java -jar /home//software/picard.jar  LiftoverVcf  I=output_pass.vcf   O=lifted_over.vcf   CHAIN=//home//database/UCSC/hg38ToHg19.over.chain   REJECT=rejected_variants.vcf   R=/home//database/NCBI/hg19/hg19.fa'

        # 4) Convert VCF to MAF using VEP annotations.
        cmd4 = 'perl /home//software/vcf2maf-1.6.19/vcf2maf.pl  --input-vcf  lifted_over.vcf    --output-maf   lifted_over.maf    --vep-path /home/zhouchi/software/ensembl-vep/  --vep-data ' \
        '/home//database/ensembl/hg19/vep-cache/   --ref-fasta   /home//database/NCBI/hg19/hg19.fa  --tumor-id {}{}  --normal-id  {}{}'.format(caseid, codeid, caseid, codeid_normal)
        
        f_subprocess([cmd1, cmd2, cmd3, cmd4])

        # Load the MAF, compute VAF, and keep a project-specific column subset.
        dat = pd.read_csv('lifted_over.maf', sep='\t', skiprows=1)
        dat['VAF'] = round(dat['t_alt_count'] / dat['t_depth'], 4)
        columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Variant_Classification', 
        'Strand', 'Variant_Type', 'HGVSp_Short', 'Consequence', 'VAF', 'dbSNP_RS', 'Existing_variation', 'BIOTYPE', 
        'SIFT', 'PolyPhen', 'IMPACT', 'vcf_pos']
        dat = dat[columns]
        if dat.shape[0] >=1:
            dat['IMPACT1'] = dat.apply(fun2, axis=1)
        else:
            dat['IMPACT1'] = 'Nan'
        dat.to_csv('lifted_overPASS.maf', sep='\t', index=False)     

        # Prepare a separate VCF subset for pyclone (SNVs + indels, tumor sample only).
        if not os.path.isdir('4pyclone_snpindels'): os.makedirs('4pyclone_snpindels')
        cmd5 = "bcftools view  -s  {}{} -v snps,indels -o 4pyclone_snpindels/output_filterPro.vcf.gz  -M2 -O z output_filter.vcf.gz".format(caseid, codeid)
        cmd6 = "bcftools filter  -i '(FMT/AF[1:0] >= 0.05  && SUM(FMT/AD[1:])>=50 &&  FMT/AD[1:1] >=7)'   output_filterPro.vcf    |  bcftools view   -v  snps,indels  -o 4pyclone_snpindels/output_pass2.vcf  -M2 -O v -"
        cmd7 = 'bcftools index  4pyclone_snpindels/output_filterPro.vcf.gz'
        f_subprocess([cmd5, cmd6, cmd7])
    except:
        print (dirName)



def f_vcfProcess():
    """
    Batch hg19-converted pipeline for all tumor tissues in the single-cell exome dataset.
    """
    dirNames = glob.glob('/home_2//outcomeSingleCellExome/*/WES/*/gatk')
    dirNames = [i for i in dirNames if os.path.basename(os.path.dirname(i)) != 'NC']
    doMultiProcess.myPool(vcfProcess, dirNames, processes=2)

def vcfProcess_hg38(dirName):
    """
    Alternative processing path that stays on hg38 and uses GRCh38 references in vcf2maf.

    This produces:
    - `4pyclone_snpindels/output_pass2.maf`
    - `4pyclone_snpindels/output_pass2absolute.maf` (absolute-compatible format)
    """
    try:
        os.chdir(dirName)
        caseid = dirName.split('/')[4]
        codeid = dirName.split('/')[6]
        codeid_normal = 'NC'
        if os.path.isfile('4pyclone_snpindels/output_pass2.vep.vcf'): os.remove('4pyclone_snpindels/output_pass2.vep.vcf')

        # Convert the already-filtered pyclone VCF to MAF on GRCh38.
        cmd4 = 'perl /home//software/vcf2maf-1.6.19/vcf2maf.pl  --input-vcf  4pyclone_snpindels/output_pass2.vcf    --output-maf   4pyclone_snpindels/output_pass2Pro.maf    --vep-path /home/zhouchi/software/ensembl-vep/  --vep-data ' \
        '/home/zhouchi/vep_data  --ref-fasta   /home//database/NCBI/hg38/Allhs_ref_GRCh38.p7.fa  --ncbi-build  GRCh38  --vep-forks 4 --tumor-id {}{}  --normal-id  {}{}'.format(caseid, codeid, caseid, codeid_normal)
        f_subprocess([cmd4])

        # Keep a consistent column subset and compute VAF + IMPACT1 for downstream analyses.
        dat = pd.read_csv('4pyclone_snpindels/output_pass2Pro.maf', sep='\t', skiprows=1)
        dat['VAF'] = round(dat['t_alt_count'] / dat['t_depth'], 4)
        columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Variant_Classification', 
        'Strand', 'Variant_Type', 'HGVSp_Short', 'Consequence', 'VAF', 'dbSNP_RS', 'Existing_variation', 'BIOTYPE', 
        'SIFT', 'PolyPhen', 'IMPACT', 'vcf_pos']
        dat = dat[columns]
        if dat.shape[0] >=1:
            dat['IMPACT1'] = dat.apply(fun2, axis=1)
        else:
            dat['IMPACT1'] = 'Nan'
        dat.to_csv('4pyclone_snpindels/output_pass2.maf', sep='\t', index=False)

        # Build a compact MAF-like table used by ABSOLUTE (chromosomes without 'chr' prefix).
        dat = pd.read_csv('4pyclone_snpindels/output_pass2Pro.maf', sep='\t', skiprows=1)
        columns = ['t_ref_count', 't_alt_count', 'dbSNP_Val_Status', 'Start_Position', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome']
        dat = dat[columns]
        dat.columns = ['t_ref_count', 't_alt_count', 'dbSNP_Val_Status', 'Start_position', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome']
        dat['Chromosome'] = dat['Chromosome'].apply(lambda x: x.strip('chr'))
        dat.to_csv('4pyclone_snpindels/output_pass2absolute.maf', sep='\t', index=False)
    except:
        print (dirName)







def f_vcfProcess_hg38():
    """Batch GRCh38 pipeline for all tumor tissues under `/home_2//outcome/*/WES/*/gatk`."""
    dirNames = glob.glob('/home_2//outcome/*/WES/*/gatk')
    dirNames = [i for i in dirNames if os.path.basename(os.path.dirname(i)) != 'NC']
    doMultiProcess.myPool(vcfProcess_hg38, dirNames, processes=10)

if __name__ == '__main__':
    print ('hello, world')
    f_vcfProcess()
    f_vcfProcess_hg38()
