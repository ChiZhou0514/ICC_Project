"""
Somatic variant calling pipeline (Sentieon TNhaplotyper2 / TNfilter).

This script contains two major stages:
1) Mapping + preprocessing (BWA-MEM, sorting, duplicate marking, recalibration table)
2) Tumor-normal somatic calling (TNhaplotyper2 + TNfilter)

It also includes an optional germline calling stage (Haplotyper) for normal samples.

Assumptions
- Directory layout resembles:
  `<patient>/WES/<tissue>/` where tumor tissues are non-"NC" and normal is "NC".
- Input FASTQs are named `A_R1.fastq.gz`, `A_R2.fastq.gz` (and optionally `B_*` for lanes).
- Outputs are written under the local `gatk/` (and `germline/`) folders.
- Sentieon installation and reference files are configured via `RunMultiProcess` in `util.py`.

Operational notes
- Commands are executed via `subprocess.Popen(..., shell=True)` and logged to `gatk/log.txt`.
- This script is designed for an HPC environment; paths are hard-coded.
"""

from util import *
import psutil

doMultiProcess = RunMultiProcess()
SENTIEON_LICENSE = doMultiProcess.SENTIEON_LICENSE
def check_license(SENTIEON_LICENSE):
    """
    Ensure Sentieon license server is running.

    Sentieon uses a local license server process (`licsrvr`). If it is not detected,
    start it using the license file path.
    """
    flag = 0
    for proc in psutil.process_iter():
        if 'licsrvr' in proc.name():
            flag = 1
    if flag == 0:
        print('start licsrvr!')
        cmd = r'sentieon-genomics-201911/libexec/licsrvr --start {}'.format(SENTIEON_LICENSE)
        subprocess.call(cmd, shell=True)


def f_subprocess(cmd):
    """
    Run a list of shell commands joined by '&&' and append output to `log.txt`.

    Parameters
    ----------
    cmd : list[str]
        Command segments; they will be joined with '&&' so failure stops the chain.
    """
    with open('log.txt', 'a') as fout:
        p = subprocess.Popen('&&'.join(cmd), shell=True, stderr=subprocess.STDOUT, stdout=fout)
        p.wait()


def removeFiles():
    """
    Remove intermediate BAM files to save space after deduplication.

    Keeps only `gatk/deduped.bam` and its index.
    """
    files = glob.glob('gatk/*bam*')
    keep_files = ['gatk/deduped.bam', 'gatk/deduped.bam.bai']
    for i in files:
        if i not in keep_files:
            os.remove(i)

def MappingAndProcess(dirName):
    """
    Map FASTQs and generate preprocessing artifacts for one sample directory.

    Parameters
    ----------
    dirName : str
        Sample directory, expected to contain FASTQs and where `gatk/` will be created.
    """
    check_license(SENTIEON_LICENSE)
    os.environ['SENTIEON_LICENSE'] = SENTIEON_LICENSE
    os.chdir(dirName)
    if not os.path.isdir('gatk'): os.makedirs('gatk')
    if os.path.isfile('gatk/deduped.bam'): return
    patient = os.path.basename(os.path.dirname(os.path.dirname(dirName)))
    tissue = os.path.basename(dirName)
    caseid = patient + tissue
    cmd = []
    try:
        # If lane B exists, map A and B separately then merge; otherwise map only A.
        if os.path.isfile('B_R1.fastq.gz'):
            cmd1 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  bwa  mem  -M -R "@RG\tID:{caseid}_A\tSM:{caseid}\tPL:ILLUMINA"  -t {nt}  -K 10000000  {fasta} A_R1.fastq.gz  A_R2.fastq.gz  ' \
                r'| {SENTIEON_INSTALL_DIR}/bin/sentieon  util sort -o  gatk/sorted_A.bam -t 64 --sam2bam -i -'.format(caseid = caseid, **doMultiProcess.__dict__) 
            cmd2 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  bwa  mem  -M -R "@RG\tID:{caseid}_B\tSM:{caseid}\tPL:ILLUMINA"  -t {nt}  -K 10000000  {fasta} B_R1.fastq.gz  B_R2.fastq.gz  ' \
                r'| {SENTIEON_INSTALL_DIR}/bin/sentieon  util sort -o  gatk/sorted_B.bam -t 64 --sam2bam -i -'.format(caseid = caseid, **doMultiProcess.__dict__) 
            cmd3 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  util  merge -r {fasta} -o  gatk/sorted.bam -i  gatk/sorted_A.bam  gatk/sorted_B.bam'.format(**doMultiProcess.__dict__)
            cmd += [cmd1, cmd2, cmd3]
        else:
            cmd1 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  bwa  mem  -M -R "@RG\tID:{caseid}_A\tSM:{caseid}\tPL:ILLUMINA"  -t {nt}  -K 10000000  {fasta} A_R1.fastq.gz  A_R2.fastq.gz  ' \
                r'| {SENTIEON_INSTALL_DIR}/bin/sentieon  util sort -o  gatk/sorted.bam -t {nt}  --sam2bam -i -'.format(caseid = caseid, **doMultiProcess.__dict__)
            cmd += [cmd1]

        # Sentieon best-practice preprocessing:
        # - LocusCollector: gather duplication scoring info
        # - Dedup: mark duplicates (produces deduped BAM)
        # - QualCal: generate recalibration table
        cmd4 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  driver -t {nt}  -i gatk/sorted.bam  --algo LocusCollector --fun score_info gatk/score.txt'.format(**doMultiProcess.__dict__)
        cmd5 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  driver -t {nt}  -i gatk/sorted.bam  --algo Dedup --score_info  gatk/score.txt --metrics gatk/dedup_metrics.txt  gatk/deduped.bam'.format(**doMultiProcess.__dict__)
        cmd6 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon  driver -r {fasta}  -t {nt}  -i gatk/deduped.bam --algo QualCal -k  {dbsnp}  -k  {known_Mills_indels}  -k {known_1000G_indels}  gatk/recal_data.table'.format(**doMultiProcess.__dict__)
        cmd7 = 'java  -jar /home//software/picard.jar  CollectWgsMetrics  R={fasta}  I=gatk/deduped.bam   O=gatk/WgsSummary.tsv  INTERVALS=/home_2//Allhs_ref_GRCh38.p7.dict'.format(**doMultiProcess.__dict__)
        cmd += [cmd4, cmd5, cmd6]
        with open('gatk/log.txt', 'a') as fout:
            p = subprocess.Popen('&&'.join(cmd), shell=True, stderr=subprocess.STDOUT, stdout=fout)
            p.wait()
        removeFiles()
    except Exception as e:
        print (e); print (dirName)


def callMutation(dirName):
    """
    Run tumor-normal somatic calling for one sample directory.

    Parameters
    ----------
    dirName : str
        Tumor tissue directory containing `gatk/deduped.bam` and `gatk/recal_data.table`.
    """
    check_license(SENTIEON_LICENSE)
    os.environ['SENTIEON_LICENSE'] = SENTIEON_LICENSE
    os.chdir(dirName)
    caseid = os.path.basename(os.path.dirname(os.path.dirname(dirName)))
    codeid = os.path.basename(dirName)
    try:
        if os.path.isfile('gatk/output_filter.vcf.gz'): return

        # TNhaplotyper2 produces an unfiltered tumor-normal VCF plus auxiliary models
        # (orientation bias and contamination) for robust filtering.
        cmd1 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t  {nt} --interval {bed}  -i  gatk/deduped.bam -i  ../NC/gatk/deduped.bam  -q  gatk/recal_data.table  -q ../NC/gatk/recal_data.table  --interval_padding 50  ' \
        ' --algo TNhaplotyper2  --min_base_qual  20  --callable_depth 20 --tumor_sample  {caseid}{codeid} --normal_sample {caseid}NC --germline_vcf {germline} --pon {pon}  gatk/output.vcf.gz' \
        ' --algo OrientationBias --tumor_sample  {caseid}{codeid}  gatk/orientation_data_file  ' \
        ' --algo ContaminationModel   --tumor_sample  {caseid}{codeid}  --normal_sample  {caseid}NC  --vcf {germline}  --tumor_segments  gatk/CONTAMINATION_DATA.segments  gatk/CONTAMINATION_DATA'.format(caseid=caseid, codeid = codeid, **doMultiProcess.__dict__)

        # TNfilter applies FP-rate control + orientation/contamination models to produce final VCF.
        cmd2 = r'{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} --algo TNfilter --tumor_sample {caseid}{codeid} --normal_sample {caseid}NC  -v  gatk/output.vcf.gz --max_fp_rate {max_fp_rate} --orientation_priors ' \
            '  gatk/orientation_data_file --contamination  gatk/CONTAMINATION_DATA   --tumor_segments  gatk/CONTAMINATION_DATA.segments   gatk/output_filter.vcf.gz '.format(caseid=caseid, codeid = codeid, **doMultiProcess.__dict__)
        with open('gatk/log.txt', 'a') as fout:
            p = subprocess.Popen(cmd1 + '&&' + cmd2, shell=True, stderr=subprocess.STDOUT, stdout=fout)
            p.wait()
    except Exception as e:
        print (e); print (dirName)

def f_MappingAndProcess():
    """Batch mapping + preprocessing for all samples (including NC) in the dataset."""
    dirNames = glob.glob('/NFS_home/NFS_home_3//outcomeSingleCellExome/*/WES/*')
    doMultiProcess.myPool(MappingAndProcess, dirNames, processes=1)

def f_callMutation():
    """Batch somatic calling for all tumor tissues (excluding NC)."""
    dirNames = glob.glob('/NFS_home/NFS_home_3//outcomeSingleCellExome/*/WES/*')
    dirNames = [i for i in dirNames if os.path.basename(i) != 'NC']
    doMultiProcess.myPool(callMutation, dirNames, processes=2)

### germline  mutation
def callMutation_G(dirName):
    """
    Call germline variants for a patient based on all available WES BAMs.

    Parameters
    ----------
    dirName : str
        Directory of the normal sample (`.../WES/NC`). Output is written into `germline/`.
    """
    os.chdir(dirName)
    if not os.path.isdir('germline'): os.makedirs('germline')
    check_license(SENTIEON_LICENSE)
    os.environ['SENTIEON_LICENSE'] = SENTIEON_LICENSE    
    patient = os.path.basename(os.path.dirname(os.path.dirname(dirName)))
    files = glob.glob('//NFS_home/NFS_home_3//outcomeSingleCellExome/{}/WES/*'.format(patient))
    files = [os.path.basename(i) for i in files]
    # Build repeated `-i <bam> -q <recal>` pairs for each tissue so Haplotyper sees all inputs.
    tmp = ['-i ../{}/gatk/deduped.bam -q ../{}/gatk/recal_data.table'.format(i, i) for i in files]
    cmd = r'{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t  {nt} --interval {bed}  {tmp} --algo Haplotyper --emit_conf=30 --call_conf=30  germline/germline.vcf'.format(tmp=' '.join(tmp), **doMultiProcess.__dict__)
    f_subprocess([cmd])
    # Apply a depth filter to keep well-supported germline calls.
    cmd = 'bcftools filter  -i "FMT/DP[*]>=50"  germline/germline.vcf  -o  germline/germlinePASS.vcf  -O  v'
    f_subprocess([cmd])


def f_callMutation_G():
    """Batch germline calling for all normals (`.../WES/NC`)."""
    dirNames = glob.glob('/NFS_home/NFS_home_3//outcomeSingleCellExome/*/WES/NC')
    doMultiProcess.myPool(callMutation_G, dirNames, processes=2)
    
if __name__ == '__main__':
    print ('hello, world')
    f_MappingAndProcess()
    f_callMutation()
    f_callMutation_G()
