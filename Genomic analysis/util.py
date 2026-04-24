"""Shared helper functions for genomic and transcriptomic analysis scripts."""

from itertools import chain
import os, re, sys, glob, subprocess, shutil
from collections import defaultdict
from numpy.lib.arraysetops import isin
from tqdm import tqdm
import numpy as np, pandas as pd
from multiprocessing import Pool


pd.options.display.float_format = '{:.3f}'.format

OUTCOME_COLUMN_ALIASES = {
    'patient_id': ['patient_id', '\u75c5\u4eba\u5e8f\u53f7'],
    'patient_group': ['patient_group', '\u60a3\u8005\u7c7b\u522b'],
    'patient_name': ['patient_name', '\u75c5\u4eba\u59d3\u540d'],
    'max_tumor_diameter': ['max_tumor_diameter', '\u80bf\u7624\u6700\u5927\u5f84'],
}

TLS_COLUMN_ALIASES = {
    'sample_id': ['sample_id', '\u7f16\u53f72'],
    'sample_id_level1': ['sample_id_level1', '\u7f16\u53f71'],
}

RAPID_PROGRESSION_VALUES = {'rapid_progression', '\u5feb\u901f\u8fdb\u5c55'}
LONG_TERM_SURVIVAL_VALUES = {'long_term_survival', '\u957f\u671f\u751f\u5b58'}
TLS_WORKBOOK_CANDIDATES = (
    'analysisOutcome/TLS_grouping.xlsx',
    'analysisOutcome/TLS\u5206\u7ec4.xlsx',
)


def _resolve_column(dataframe, aliases):
    for name in aliases:
        if name in dataframe.columns:
            return name
    raise KeyError('Missing expected columns: {}'.format(aliases))


def _rename_columns_by_alias(dataframe, alias_dict):
    rename_map = {}
    for standard_name, aliases in alias_dict.items():
        for alias in aliases:
            if alias in dataframe.columns:
                rename_map[alias] = standard_name
                break
    if rename_map:
        dataframe = dataframe.rename(columns=rename_map)
    return dataframe


def _normalize_outcome_dataframe(dataframe):
    return _rename_columns_by_alias(dataframe, OUTCOME_COLUMN_ALIASES)


def _normalize_tls_dataframe(dataframe):
    return _rename_columns_by_alias(dataframe, TLS_COLUMN_ALIASES)


def _resolve_tls_workbook():
    for path in TLS_WORKBOOK_CANDIDATES:
        if os.path.isfile(path):
            return path
    return TLS_WORKBOOK_CANDIDATES[0]


class RunMultiProcess(object):
    def __init__(self):
        self.SENTIEON_LICENSE="sentieon-genomics-201911/Biological_and_Medical_Big_Data_Mining_Group_of_Tongji_University_cluster.lic"
        self.SENTIEON_INSTALL_DIR="sentieon_202010/sentieon-genomics-202010/"
        self.dbsnp="database/NCBI/hg38/gatk/dbsnp_146.hg38.vcf.gz"
        self.fasta='database/NCBI/hg38/bwa/Allhs_ref_GRCh38.p7.fa'
        self.known_Mills_indels="database/NCBI/hg38/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        self.known_1000G_indels="database/NCBI/hg38/gatk/1000G_phase1.indels.hg38.sites.vcf.gz"
        self.pon="database/NCBI/hg38/gatk/mutect2/pon.hg38.vcf.gz"
        self.germline="database/NCBI/hg38/gatk/mutect2/af-only-gnomad.hg38.vcf.gz"
        self.bed="MedExome_hg38_capture_targets.bed"
        self.refFlat = 'database/UCSC/hg38/refFlat.txt'
        self.mappable_file = 'database/UCSC/hg38/access.hg38.bed'
        self.max_fp_rate=0.01
        self.nt = 64
        self.gc = 'database/hg38.gc50Base.wig.gz'
    
    def myPool(self, func, mylist, processes):
        """
        Run a function over a list using multiprocessing with a tqdm progress bar.

        Parameters
        ----------
        func : callable
            Function to apply to each element of `mylist`.
        mylist : list
            Items to process.
        processes : int
            Number of worker processes.

        Returns
        -------
        list
            Collected results in the same order as produced by `imap`.
        """
        with Pool(processes) as pool:
            results = list(tqdm(pool.imap(func, mylist), total=len(mylist)))
        return results

def getOutcome():
    """
    Load patient outcome table and split patients into short vs long survival groups.

    This function supports both English and legacy Chinese column names via
    `_normalize_outcome_dataframe`.

    Returns
    -------
    (list[str], list[str])
        `shortSur`, `longSur` patient IDs that exist under `/home_2//outcome/<patient>/`.
    """
    dat = pd.read_excel('/home_2//80011045_outcome_20201126.xlsx')
    dat = _normalize_outcome_dataframe(dat)
    patient_id_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_id'])
    patient_group_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_group'])

    dat[patient_id_col].fillna(method='ffill', inplace=True)
    patient_group = dat[patient_group_col].astype(str).str.strip()
    short_mask = patient_group.isin(RAPID_PROGRESSION_VALUES)
    long_mask = patient_group.isin(LONG_TERM_SURVIVAL_VALUES)

    shortSur = [str(int(i)) for i in dat[short_mask][patient_id_col].unique()]
    longSur = [str(int(i)) for i in dat[long_mask][patient_id_col].unique()]
    shortSur = [i for i in shortSur if os.path.isdir('/home_2//outcome/{}'.format(i))]
    longSur = [i for i in longSur if os.path.isdir('/home_2//outcome/{}'.format(i))]
    return shortSur, longSur

def getSurvival():
    """
    Load survival metadata table indexed by patient ID.

    Returns
    -------
    pandas.DataFrame
        Indexed by `patient_id` with columns: survival, os, recurrence, TTR.
    """
    dat = pd.read_excel('/home_2//outcome_survival.xlsx')
    dat = _normalize_outcome_dataframe(dat)
    patient_id_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_id'])
    dat = dat[[patient_id_col, 'survival', 'os', 'recurrence', 'TTR']]
    dat = dat.rename(columns={patient_id_col: 'patient_id'})
    dat.index = dat['patient_id']
    return dat

def getSurvivalTumorSize():
    """
    Load tumor size (max diameter) table indexed by patient ID.
    """
    dat = pd.read_excel('/home_2//outcome_survival.xlsx')
    dat = _normalize_outcome_dataframe(dat)
    patient_id_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_id'])
    max_tumor_diameter_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['max_tumor_diameter'])
    dat = dat[[patient_id_col, max_tumor_diameter_col]]
    dat = dat.rename(columns={patient_id_col: 'patient_id', max_tumor_diameter_col: 'max_tumor_diameter'})
    dat.index = dat['patient_id']
    return dat

def getTCGAOutcome():
    """
    Load TCGA clinical outcome table for cases that have local SNV/Indel directories.

    Returns
    -------
    pandas.DataFrame
        Filtered clinical table containing case_submitter_id and days_to_death.
    """
    dat = pd.read_csv('/home_2//ICC_TCGA/clinical.tsv', sep='\t')
    dat = dat[['case_submitter_id', 'days_to_death']]
    dat = dat[dat['days_to_death'] != "'--"]
    dat.drop_duplicates(subset='case_submitter_id', inplace=True)
    maf = pd.read_csv('/home_2//maftools/TCGA/allSample.maf', sep='\t')
    patients = maf['Tumor_Sample_Barcode'].unique()
    a = dat['case_submitter_id'].apply(lambda x: True if os.path.isdir('/home_2//ICC_TCGA/SNVIndel/{}'.format(x)) else False).values
    dat = dat.loc[a, :]
    return dat

def getDriveGene(info ='driver q-value', value = 0.1): 
    """
    Load driver-gene predictions and return genes passing a q-value cutoff.

    Parameters
    ----------
    info : str
        Column name used for sorting/cutoff (e.g., 'driver q-value').
    value : float
        Maximum allowed value in `info`.
    """
    filein = '/home_2//driverGene/merge/pretrained_output/results/r_random_forest_prediction.txt'
    dat = pd.read_csv(filein, sep='\t')
    dat.sort_values(info, ascending = True, inplace=True)
    dat = dat[dat[info] <= value]
    return sorted(dat['gene'].tolist())

def getDriveInto(Nums1=5, Nums2=2):
    """
    Filter IntOGen driver genes by minimum sample and cohort counts.
    """
    dat = pd.read_csv('/home_2//IntOGen-DriverGenes_CH.tsv', sep='\t')
    dat = dat[(dat['Samples'] >= Nums1) & (dat['Cohorts'] >= Nums2)]
    return sorted(dat['Symbol'].tolist())

def getDriver():
    """
    Load Cancer Gene Census and return genes annotated as somatic drivers.
    """
    dat = pd.read_csv('/home//database/cancer_gene_census.csv', sep=',')
    dat = dat[dat['Somatic'] == 'yes']
    dat = dat[dat['Role in Cancer'].isin(['TSG', 'oncogene', 'oncogene, fusion', 'TSG, fusion', 'oncogene, TSG', 'oncogene, TSG, fusion'])]
    return sorted(dat['Gene Symbol'].values.tolist())


def getDriver1():
    """
    Convenience helper: pick the top recurrent driver genes (project heuristic).
    """
    maf = pd.read_csv('/home_2//analysisMetAF5/maftools/sampleLevel/allSample.maf', sep='\t')
    maf1 = maf[maf['IMPACT'].isin(['HIGH', 'MODERATE'])]
    driverGene3 = getDriver()
    most_com = maf1['Hugo_Symbol'].value_counts().head(n=30)
    return [i for i in most_com.index if i in driverGene3]



def getProteinCoding():
    """Return a list/array of protein-coding genes from a local Excel reference."""
    protein = pd.read_excel('/home//database/Genes.xlsx')['Gene_Symbol'].values
    return protein

def getPatientNum():
    """
    Build a mapping from patient name to numeric patient ID.

    Returns
    -------
    dict
        patient_name -> patient_id (int)
    """
    mydict = {}
    dat = pd.read_excel('80011045_outcome_20201126.xlsx')
    dat = _normalize_outcome_dataframe(dat)
    patient_id_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_id'])
    patient_name_col = _resolve_column(dat, OUTCOME_COLUMN_ALIASES['patient_name'])
    dat[patient_id_col].fillna(method='ffill', inplace=True)
    for i, j in zip(dat[patient_name_col], dat[patient_id_col]):
        mydict[i] = int(j)
    return mydict

def getPatientFiles(patient):
    """
    List tissue/sample folders for a given patient, excluding `NC`.

    The list is sorted such that PT* samples come first (project convention).
    """
    files = glob.glob('outcome/{}/WES/*'.format(patient))
    files = [os.path.basename(i) for i in files]; files = [i for i in files if i != 'NC']
    files = sorted(files); files = [i for i in files if i.startswith('PT')] + [i for i in files if not i.startswith('PT')]
    return files

def getDataframe():
    """
    Create a long-form table of (patient, tissue, outcome) for all patients.

    Returns
    -------
    pandas.DataFrame
        Columns: patient, tissue, outcome, Multi
    """
    a =[]; b= []; c = []; d = []
    shortSur, longSur = getOutcome()
    for patient in chain(shortSur, longSur):
        files = getPatientFiles(patient)
        if patient in shortSur:
            a.extend([patient] * len(files)); b.extend(files)
        else:
            c.extend([patient] * len(files)); d.extend(files)
    dat = pd.DataFrame({'patient': a + c, 'tissue': b + d, 'outcome': ['shortSur'] * len(a) + ['longSur'] * len(c) })
    dat['Multi'] = [True if len(getPatientFiles(patient)) >=2 else False for patient in dat['patient']]
    return dat


def getDataframe1():
    """
    Join outcome dataframe with TLS grouping annotations.

    Returns
    -------
    pandas.DataFrame
        `getDataframe()` rows merged with TLS grouping Excel on `sample_id`.
    """
    dat = getDataframe()
    dat['sample_id'] = dat['patient'] + '_' + dat['tissue']
    dat1 = pd.read_excel(_resolve_tls_workbook())
    dat1 = _normalize_tls_dataframe(dat1)
    sample_id_col = _resolve_column(dat1, TLS_COLUMN_ALIASES['sample_id'])
    if sample_id_col != 'sample_id':
        dat1 = dat1.rename(columns={sample_id_col: 'sample_id'})
    dat2 = pd.merge(dat, dat1, on='sample_id')
    return dat2

def gettumorContent(patient, tissue):
    """
    Load tumor purity and ploidy from Sequenza output text files.
    """
    filein1 = 'outcome/{}/WES/{}/sequenza/{}_cellularity.txt'.format(patient, tissue, patient)
    filein2 = 'outcome/{}/WES/{}/sequenza/{}_ploidy.txt'.format(patient, tissue, patient)
    purity = subprocess.getoutput('cat {}'.format(filein1))
    ploidy = subprocess.getoutput('cat {}'.format(filein2))
    return  [np.float32(purity), np.float32(ploidy)]


def f_gettumorContent(isint=False):
    """
    Compute a nested dictionary of tumor purity and ploidy for all patients/tissues.

    Returns
    -------
    dict
        mydict[patient][tissue] = [purity, ploidy]
    """
    allPatient = os.listdir('outcome')
    mydict = defaultdict(dict)
    for patient in allPatient:
        files = getPatientFiles(patient)
        for tissue in files:
            purity, ploidy = gettumorContent(patient, tissue)
            if isint: ploidy = round(ploidy)
            mydict[patient][tissue] = [purity, ploidy]
    return mydict


def getTLSStatus():
    """
    Return hard-coded sample IDs labeled as TLS-negative vs TLS-positive for scRNA data.
    """
    noTLS = ['T213082', 'T220002_3', 'T220008', 'T220010', 'T220421']
    withTLS = ['T213087', 'T220011', 'T220012_4', 'T220442', 'T220452']
    return noTLS, withTLS

    
def getDataframeSC():
    """
    Build a minimal dataframe for single-cell samples with a TLS label.
    """
    noTLS, withTLS = getTLSStatus()
    dat = pd.DataFrame({'patient': noTLS + withTLS, 'TLS': ['noTLS'] * len(noTLS) + ['withTLS'] * len(withTLS)})
    dat['tissue'] = 'PT'
    return dat
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
