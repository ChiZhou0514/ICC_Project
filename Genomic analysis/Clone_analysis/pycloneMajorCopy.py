# -*- coding: utf-8 -*- 
"""
Pyclone-based clonal deconvolution using parental copy number (major/minor CN).

This script is a project-specific wrapper that prepares inputs for PyClone, runs or
post-processes results, and converts outputs into formats used by downstream tools
such as clonevol and REVOLVER.

Main capabilities
1) Merge per-tissue VCFs into a combined VCF for each patient (`mergeVCF`)
2) Extract per-mutation ref/alt counts from merged VCF, and annotate each mutation
   with (normal_cn, minor_cn, major_cn) derived from CNVkit outputs (`prePyclone`)
3) Convert PyClone outputs into clonevol input tables (`preClonevol`)
4) Run clonevol plotting scripts across patients (`plotPyclone`)

Assumptions / environment
- Designed for an HPC/NFS directory layout with hard-coded paths under `/home_2//`.
- Requires upstream steps to generate filtered VCFs, CNVkit calls, and PyClone outputs.
"""
from util import *
import os, re, sys, glob, subprocess
import numpy as np, pandas as pd
from multiprocessing import Pool
from itertools import chain
from collections import defaultdict
import shutil, h5py, signal

### use pyclone to call tumor evolution



### use parental copy number
doMultiProcess = RunMultiProcess()

### conda activate pyclone
def f_subprocess(cmd):
    """
    Run shell command segments joined by '&&' and append output to `log.txt`.

    Parameters
    ----------
    cmd : list[str]
        Command segments that will be joined with '&&'. A failure in one segment
        prevents subsequent segments from running.
    """
    with open('log.txt', 'a') as fout:
        p = subprocess.Popen('&&'.join(cmd), shell=True, stderr=subprocess.STDOUT, stdout=fout)
        p.wait()

def f_subprocess_time(cmd):
    """
    Run a shell command with a 1-hour timeout; kill the process group on timeout.
    """
    p = subprocess.Popen(cmd, shell=True, preexec_fn=os.setsid)
    try:
        p.wait(timeout = 3600)
    except:
        print (subprocess.getoutput('pwd'))
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

def mergeVCF(dirName):
    """
    Merge per-tissue VCFs for one patient into a single `merge.vcf`.

    Parameters
    ----------
    dirName : str
        Patient WES root directory, e.g. `/home_2//outcome/<patient>/WES`.
    """
    try:
        patient = os.path.basename(os.path.dirname(dirName))
        wkdir = '/home_2//analysisOutcome/pyclone_majorCopy/{}'.format(patient)
        if not os.path.isdir(wkdir): os.makedirs(wkdir)
        os.chdir(wkdir)
        files = getPatientFiles(patient)
        if len(files) == 1:
            cmd = 'gunzip -c -k {}/{}/gatk/4pyclone_snpindels/output_filterPro.vcf.gz > merge.vcf'.format(dirName, files[0])
        else:
            files1 = ' '.join([dirName + '/' +i + '/gatk/4pyclone_snpindels/output_filterPro.vcf.gz' for i in files])
            #cmd = 'bcftools merge  {}  | bcftools view -o merge.vcf -M2 -O v - '.format(files1)
            cmd = 'bcftools merge  {}  | bcftools view -o merge.vcf  -O v - '.format(files1)
        subprocess.call(cmd, shell=True)
    except:
        print (dirName)    

def f_mergeVCF():
    """Batch `mergeVCF` across all patients under `/home_2//outcome/*/WES`."""
    dirNames = glob.glob('/home_2//outcome/*/WES')
    doMultiProcess.myPool(mergeVCF, dirNames, processes=6)


def getPos_info(patient):
    """
    Parse `merge.vcf` to extract per-sample ref/alt read counts for each variant.

    Returns
    -------
    dict
        Mapping `<chr>_<pos>` -> { sample_name: [ref_count, alt_count] }.

    Notes
    - This parser assumes the genotype field uses an AD-like string at index 1 after
      splitting by ':' (project-specific VCF layout).
    - Multi-allelic sites are handled by selecting the first non-missing alt depth.
    """
    filein = '/home_2//analysisOutcome/SCpyclone_majorCopy/{}/merge.vcf'.format(patient)
    mydict = defaultdict(dict)
    def getinfo(info):
        infos = info.split(':')[1].split(',')
        if infos[0] == '.':
            ref_count, alt_count = '.', '.'
        else:
            ref_count, alt_count = infos[0], infos[1]
        return ref_count, alt_count
    
    def getinfo1(info):
        infos = info.split(':')[1].split(',')
        if infos[0] == '.':
            ref_count, alt_count = '.', '.'
        else:
            ref_count = infos[0]
            alt_count = [i for i in infos[1:] if i !='.'][0]
        return ref_count, alt_count    

    with open(filein, 'r') as fin:
        for line in fin:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                lines = line.strip().split('\t')
                files = lines[9:]
            else:
                lines = line.strip().split('\t')
                if lines[0] in ['chrX', 'chrY', 'chrM']:
                    continue
                chr, pos = lines[0], lines[1]
                if len(lines[4].split(',')) == 1:
                    for i, j in zip(lines[9:], files):
                        ref_count, alt_count = getinfo(i)
                        mydict[chr + '_' + pos][j] = [ref_count, alt_count]
                else:
                    for i, j in zip(lines[9:], files):
                        ref_count, alt_count = getinfo1(i)
                        mydict[chr + '_' + pos][j] = [ref_count, alt_count]

    return mydict


def tree():
    """Recursive defaultdict factory used to store CNV intervals."""
    return defaultdict(tree)

def getCNV_info(patient, files):
    """
    Load CNVkit allele-specific CN calls for each tissue and store by genomic interval.

    Returns
    -------
    defaultdict
        Nested mapping: mydict[tissue][chr]["start_end"] -> [minor_cn, major_cn]
    """
    mydict = defaultdict(tree)
    for i in files:
        filein = '/home_2//outcomeSingleCellExome/{}/WES/{}/cnv/tumor.callGermline.cns'.format(patient,i)
        with open(filein, 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if lines[0] in ['chrX', 'chrY', 'chrM'] or not lines[7] or not lines[8]: continue
                if int(lines[7]) <= 0:  lines[7] = 1
                if int(lines[8]) <= 0 : lines[8] = 0
                if int(lines[7]) < int(lines[8]): continue
                mydict[i][lines[0]][lines[1] + '_' + lines[2]] = [lines[7], lines[8]]
    return mydict

def getTotal_cn(chr, pos, file, mydict):
    """
    Return (minor_cn, major_cn) for a given position by searching CNV intervals.

    Falls back to (1, 1) if no interval matches.
    """
    for i in mydict[file][chr]:
        start, end = i.split('_')
        if  int(start) <= int(pos) <= int(end):
            return mydict[file][chr][i]
    return 1, 1


def pos2Gene(patient, files):
    """
    Map `<chr>_<pos>` to a compact mutation identifier used by PyClone.

    The identifier is assembled from several MAF columns (project-specific).
    """
    mydict = {}
    for i in files:
        filein = '/home_2//outcomeSingleCellExome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2.maf'.format(patient,i)
        with open(filein, 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                chr, pos = lines[2], lines[-2]
                mydict[chr + '_' + pos]  = '_'.join([lines[0], pos, lines[9], lines[11], lines[-3]])
    return mydict

def prePyclone(dirName):
    """
    Prepare per-sample PyClone input tables for one patient.

    This function combines:
    - read count information from a merged VCF (`getPos_info`)
    - allele-specific copy numbers from CNVkit (`getCNV_info`)
    - mutation identifiers derived from a MAF file (`pos2Gene`)

    Output files
    - `<tissue>.txt` (one per tissue in `files`) with columns:
      mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn
    """
    try:
        mydict = defaultdict(dict)
        patient = os.path.basename(os.path.dirname(dirName))
        wkdir = '/home_2//analysisOutcome/SCpyclone_majorCopy/{}'.format(patient)
        os.chdir(wkdir)
        files = getPatientFiles(patient)
        files = ['PT']
        posinfo = getPos_info(patient)
        CNVinfo = getCNV_info(patient, files)
        for i in files:
            vcf = '/home_2//outcomeSingleCellExome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2.vcf'.format(patient, i)
            with open(vcf, 'r') as fin:
                for line in fin:
                    if not line.startswith('#'):
                        lines = line.strip().split('\t')
                        chr, pos = lines[0], lines[1]
                        if chr + '_' + pos in mydict:
                            continue
                        if chr + '_' + pos in posinfo:
                            for j in files:
                                ref_count, alt_count = posinfo[chr + '_' + pos][patient + j]
                                if ref_count == '.':
                                    vcf_tmp = '/home_2//outcomeSingleCellExome/{}/WES/{}/gatk/deduped.bam'.format(patient, j)
                                    cmd = 'samtools  depth  -r {}:{}-{}  {}'.format(chr, pos, pos, vcf_tmp)
                                    record = subprocess.getoutput(cmd).split('\t')
                                    if len(record) == 3: ref_count = record[2]
                                    else: ref_count = 0
                                    alt_count = np.random.randint(low=0, high=1)
                                minor, major = getTotal_cn(chr, pos, j, CNVinfo)
                                minor, major = min(int(minor), int(major)), max(int(minor), int(major))
                                mydict[chr + '_' + pos][j] = [str(ref_count), str(alt_count), '2', str(minor), str(major)]
        posToGene = pos2Gene(patient, files)
        for i in files:
            with open('{}.txt'.format(i), 'w') as fout:
                fout.write('mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n')
                for j in mydict:
                    fout.write('{}\t{}\n'.format(posToGene[j], '\t'.join(mydict[j][i])))
            
        with open('pyclone_vi.txt'.format(i), 'w') as fout:
            fout.write('mutation_id\tsample_id\tref_counts\talt_counts\tnormal_cn\tminor_cn\tmajor_cn\ttumour_content\n')
            for i in files:
                pass
                #tumour_content = tumorContentDict[patient][i][0]
                #for j in mydict:
                #    fout.write('{}\t{}\t{}\t{}\n'.format(posToGene[j], i, '\t'.join(mydict[j][i]), tumour_content))
    except Exception as e:
        print (e)
        print (dirName)


def f_prePyclone():
    """Batch `prePyclone` across all patients under `/home_2//outcome/*/WES`."""
    dirNames = glob.glob('/home_2//outcome/*/WES')
    doMultiProcess.myPool(prePyclone, dirNames, processes= 4)


### covert to clonevol software format
def preClonevol(dirName, ClusterNums = 10, ClusterNums1 = 10):
    """
    Convert PyClone outputs into clonevol input (`clonevol.txt`).

    Parameters
    ----------
    dirName : str
        Patient folder containing `pyclone_analysis/tables/`.
    ClusterNums : int
        Minimum number of mutations in a cluster to keep.
    ClusterNums1 : int
        Maximum number of clusters to keep after sorting by mean CCF.
    """
    try:
        driverGene1 = getDriveGene(value=.1); driverGene2 = getDriveInto(); driverGene3 = getDriver()
        driverGenes = driverGene1 + driverGene2
        mydict1 = {}; mydict2 = defaultdict(dict); mydict3 = {}
        os.chdir(dirName)
        if not os.path.isfile('pyclone_analysis/tables/loci.tsv'):
            print (dirName); return
        patient = os.path.basename(dirName)
        #files = getPatientFiles(patient)
        files = ['PT']
        with open('pyclone_analysis/tables/loci.tsv', 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                mydict2[lines[0]][lines[1]] = float(lines[3])   ### ccf
                mydict3[lines[0]] = lines[2]        ###  cluster
        
        dat = pd.read_csv('pyclone_analysis/tables/cluster.tsv', sep='\t')
        dat = dat[ (dat['size'] >= ClusterNums) & (dat['sample_id'] == files[0])]
        
        dat.sort_values(by='mean', ascending=False, inplace=True)
        if dat.shape[0] <=0:
            print (dirName); return
        if dat.shape[0] >= ClusterNums1:
            print (dirName)
            dat = dat.iloc[:ClusterNums1, :] 
        for index, i in enumerate(dat['cluster_id']):
            mydict1[str(i)] = index + 1
        with open('clonevol.txt', 'w') as fout:
            fout.write('ID\tcluster\tgene\tDriver\t{}\n'.format('\t'.join(files)))
            for key in mydict3:
                if mydict3[key] in mydict1:
                    cluster = mydict1[mydict3[key]]
                    gene = key.split('_')[0]; effect = key.split('_')[-1]
                    a = 'TRUE' if gene in driverGenes and effect in ['HIGH', 'MODERATE'] else 'FALSE'
                    ccf = [str(round(float(mydict2[key][i]) * 100, 2)) for i in files]
                    fout.write('{}\t{}\t{}\t{}\n'.format(key, cluster,'\t'.join([gene, a]), '\t'.join(ccf)))
        datRaw = pd.read_csv('clonevol.txt', sep='\t')
        datRaw.sort_values(by='cluster', ascending=True, inplace=True)
        datRaw.to_csv('clonevol.txt', sep='\t', index=False)
    except Exception as e:
        print (dirName, e)   

def f_preClonevol():
    """Batch clonevol conversion for all patient folders in SCpyclone outputs."""
    dirNames = glob.glob('/home_2//analysisOutcome/SCpyclone_majorCopy/*')
    for dirName in dirNames:
        preClonevol(dirName, ClusterNums=10)

def plotPyclone(patient):
    """Run the clonevol plotting R script for a single patient."""
    try:
        os.chdir('/home_2//analysisOutcome/pyclone_majorCopy/{}'.format(patient))
        cmd = '/usr/local/bin/Rscript  /home_2//clonevol.r'
        f_subprocess([cmd])
    except:
        print (patient)

def f_plotPyclone():
    """Batch plotting for patients with multiple tissues (>=2)."""
    patients = os.listdir('/home_2//outcome')
    patients = [patient for patient in patients if len(getPatientFiles(patient)) >=2]
    doMultiProcess.myPool(plotPyclone, patients, processes=20)

def f_plotPyclone1():
    """Batch plotting for patients with a single tissue (==1)."""
    patients = os.listdir('/home_2//outcome')
    patients = [patient for patient in patients if len(getPatientFiles(patient)) ==1]
    doMultiProcess.myPool(plotPyclone, patients, processes=20)

def preOneSample():
    """
    Utility: duplicate a sample column for single-sample patients in clonevol tables.

    Some clonevol workflows expect at least two samples. This duplicates `PT1-1` into
    `PT1-2` for patients with a single tissue to enable plotting.
    """
    patients = [patient for patient in patients if len(getPatientFiles(patient)) ==1]
    for patient in patients:
        filein = '/home_2//analysisOutcome/pyclone_majorCopy/{}/clonevol.txt'.format(patient)
        dat = pd.read_csv(filein, sep='\t')
        dat['PT1-2'] = dat['PT1-1']
        dat.to_csv(filein, sep='\t', index=False)


if __name__ == '__main__':
    print ('hello, world')
    f_mergeVCF()
    f_prePyclone()
    f_preClonevol()
    f_plotPyclone()
    f_plotPyclone1()
