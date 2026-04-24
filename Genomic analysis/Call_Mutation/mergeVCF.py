"""
MAF merging and cohort-level table preparation for downstream analysis.

This script aggregates per-tissue MAF files into:
- Patient-level MAFs (one row per unique variant across all tissues for a patient)
- Region-level MAFs (one row per tissue/sample, with optional clonality labels)
- Clinical annotation tables for maftools plotting (clinical.tsv / clinical1.tsv)

It also includes helper routines for:
- Adding clonality labels (`clone_state9.tsv`) to MAF rows
- Merging signatures across outcome/metastasis datasets
- Small expression sanity-checks (anaExp)

Assumptions
- Path layout is hard-coded under `/home_2//outcome` and `/home_2//analysisOutcome`.
- Many intermediate files must already exist (e.g., `output_pass2.maf`, `clone_state9.tsv`).
"""

from itertools import chain
from util import *
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import stats

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rc('font',family='DejaVu Sans Mono')
os.chdir('/home_2//')

def CountMutation5():
    """
    Correlate mutation counts across multi-region samples.

    This is a quick exploratory utility that:
    - Collects variant counts across tissue pairs for patients with >=2 tissues
    - Computes Pearson correlation and draws a robust regression plot
    """
    a = []; b = []
    patients = os.listdir('/home_2//outcome')
    for patient in patients:
        files = getPatientFiles('outcome', patient)
        if len(files) >=2:
            for i in range(len(files) - 1):
                for j in range(i+1, len(files)):
                    maf1 = '/home_2//{}/{}/WES/{}/gatk/lifted_overPASS.maf'.format('outcome', patient, files[i])
                    maf2 = '/home_2//{}/{}/WES/{}/gatk/lifted_overPASS.maf'.format('outcome', patient, files[j])
                    count1 = pd.read_csv(maf1, sep='\t').shape[0]
                    count2 = pd.read_csv(maf2, sep='\t').shape[0]
                    a.append(count1); b.append(count2)
    stats.pearsonr(a, b)
    sns.regplot(a, b, robust=True)

    a = []; b = []
    files = glob.glob('/home_2//outcome/*/WES/PT1-1/gatk/lifted_overPASS.maf')
    for i in range(len(files) - 1):
        for j in range(i+1, len(files)):
            count1 = pd.read_csv(files[i], sep='\t').shape[0]
            count2 = pd.read_csv(files[j], sep='\t').shape[0]
            a.append(count1); b.append(count2)
    stats.pearsonr(a, b)
    sns.regplot(a, b, robust=True)


def RegionLevel():
    """
    Create region-level maftools inputs:
    - `clinical1.tsv` labeled by outcome (shortSur/longSur)
    - `allSample_final.maf` merged from all region-level MAFs
    """
    os.chdir('/home_2//analysisOutcome/maftools/regionLevel')
    dat1 = getDataframe()
    dat1['Tumor_Sample_Barcode'] = dat1['patient'] + dat1['tissue']
    dat1 = dat1[['Tumor_Sample_Barcode', 'outcome']]
    dat1.columns = ['Tumor_Sample_Barcode', 'group']
    dat1.to_csv('clinical1.tsv', sep='\t', header=True, index=False)
    alldat = []
    mafs = glob.glob('/home_2//outcome/*/WES/PT*/gatk/4pyclone_snpindels/output_pass3.maf')
    for i in mafs:
        dat = pd.read_csv(i, sep='\t')
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('allSample_final.maf', sep='\t', header=True, index=False)

def RegionLevel1():
    """
    Add clonality label to each MAF row by joining to `clone_state9.tsv`.

    For each patient/tissue:
    - Read `clone_state9.tsv` into a position->clonality map
    - Add a `clonality` column into `output_pass2.maf`
    - Write `output_pass2Clonality.maf`
    """
    dat1 = getDataframe()
    def fun(x):
        pos = x['Hugo_Symbol'] + '_' + str(x['vcf_pos'])
        if pos in mydict:
            return mydict[pos]
        else:
            return 'NOT'
    for patient, tissue in zip(dat1['patient'], dat1['tissue']):
        mydict = {}
        filein = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/clone_state9.tsv'.format(patient, tissue)
        with open(filein, 'r') as fin:
            for line in fin:
                lines = line.strip().split('\t')
                mydict[lines[0] + '_' +  lines[2]] = lines[3]
        filein = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2.maf'.format(patient, tissue)
        fileout = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2Clonality.maf'.format(patient, tissue)
        dat = pd.read_csv(filein, sep='\t')
        dat['clonality'] = dat.apply(fun, axis=1)
        dat.to_csv(fileout, sep='\t', header=True, index=False)

def f_RegionLevel1():
    """
    Merge all region-level clonality MAFs into a single cohort file and write clinical labels.
    """
    os.chdir('/home_2//analysisOutcome/maftools/regionLevel')
    dat1 = getDataframe1()
    dat1['Tumor_Sample_Barcode'] = dat1['patient'] + dat1['tissue']
    dat1 = dat1[['Tumor_Sample_Barcode', 'Type']]
    dat1.columns = ['Tumor_Sample_Barcode', 'group']
    dat1.to_csv('clinical.tsv', sep='\t', header=True, index=False)
    alldat = []
    dat1 = getDataframe1()
    for patient, tissue in zip(dat1['patient'], dat1['tissue']):
        filein = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2Clonality.maf'.format(patient, tissue)
        dat = pd.read_csv(filein, sep='\t')
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('allSampleClonality.maf', sep='\t', header=True, index=False)



def mergeMAFSignature():
    """
    Merge MAF files across datasets for signature-level analyses.

    Output
    - `/home_2//maftools/merge/allSample_PT.maf`
    """
    alldat = []
    os.chdir('/home_2//maftools/merge')
    mafs1 = ['/home_2//maftools/outcome_merge/allSample.maf']
    mafs2 = glob.glob('/home_2//metastasis/*/WES/PT/gatk/4pyclone_snpindels/output_pass2.maf') + glob.glob('/home_2//metastasis/*/WES/PT1-1/gatk/4pyclone_snpindels/output_pass2.maf')
    for i in chain(mafs1, mafs2):
        dat = pd.read_csv(i, sep='\t')
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('allSample_PT.maf', sep='\t', header=True, index=False)


def fun1(patient):
    """
    Patient-level merge helper:
    Combine all tissues for a patient into one MAF and de-duplicate by (gene, start_position).
    """
    alldat = []
    files = getPatientFiles(patient)
    for tissue in files:
        filein = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/output_pass2.maf'.format(patient, tissue)
        dat = pd.read_csv(filein, sep='\t', dtype=np.str)
        alldat.append(dat)
    alldat = pd.concat(alldat, axis=0)
    alldat['tmp'] = alldat['Hugo_Symbol'] + '_' + alldat['Start_Position']
    alldat.drop_duplicates('tmp', inplace=True)
    alldat.drop(labels='tmp', axis=1, inplace=True)
    alldat['Tumor_Sample_Barcode'] = patient
    return alldat

def PatientLevel():
    """
    Build patient-level maftools files:
    - `short.maf`, `long.maf`, `allSample.maf`
    - `clinical.tsv` labeling each patient as short/long
    """
    os.chdir('/home_2//analysisOutcome/maftools/patientLevel')
    shortSur, longSur = getOutcome()
    alldatshort = []; alldatlong = []
    for patient in shortSur:
        alldat = fun1(patient)
        alldatshort.append(alldat)
    alldatshort = pd.concat(alldatshort, axis=0)
    for patient in longSur:
        alldat = fun1(patient)
        alldatlong.append(alldat)
    alldatlong = pd.concat(alldatlong, axis=0)
    alldatshort.to_csv('short.maf', sep='\t', header=True, index=False)
    alldatlong.to_csv('long.maf', sep='\t', header=True, index=False)
    dat1 = alldatlong.loc[:, ['Tumor_Sample_Barcode']]; dat1['group'] = 'long'
    dat2 = alldatshort.loc[:, ['Tumor_Sample_Barcode']]; dat2['group'] = 'short'
    dat = pd.concat([dat1, dat2], axis=0)
    dat.drop_duplicates(subset='Tumor_Sample_Barcode', keep='first', inplace=True)
    dat.to_csv('clinical.tsv', sep='\t', header=True, index=False)
    alldat = pd.concat([alldatlong, alldatshort], axis=0)
    alldat.to_csv('allSample.maf', sep='\t', header=True, index=False)


def PatientLevelClonality():
    """
    Add clonality labels to patient-level MAFs from precomputed `clone_state9.tsv`.

    Reads:
    - `/home_2//analysisOutcome/cloneStatus/<patient>/clone_state9.tsv`
    Writes:
    - `patientLevel/<group>Clonality.maf` for each group in {allSample, short, long}
    """
    def fun(x):
        pos = x['Hugo_Symbol'] + '_' + str(x['vcf_pos'])
        if pos in mydict[str(x['Tumor_Sample_Barcode'])]:
            return mydict[str(x['Tumor_Sample_Barcode'])][pos]
        else:
            return 'NOT'
    mydict = defaultdict(dict)
    patients = os.listdir('/home_2//outcome')
    for patient in patients:
        filein = '/home_2//analysisOutcome/cloneStatus/{}/clone_state9.tsv'.format(patient)
        with open(filein, 'r') as fin:
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                mydict[patient][lines[0] + '_' + lines[2]] = lines[3]
    
    for i in ['allSample', 'short', 'long']:
        dat = pd.read_csv('/home_2//analysisOutcome/maftools/patientLevel/{}.maf'.format(i), sep='\t')
        dat['clonality'] = dat.apply(fun, axis=1)
        dat = dat[dat['clonality'] != 'NOT']
        dat.to_csv('/home_2//analysisOutcome/maftools/patientLevel/{}Clonality.maf'.format(i), sep='\t', index=False)



def anaExp():
    """
    Small expression sanity-check helper (paired and unpaired t-tests).

    This reads an edgeR-style count matrix and metadata, then performs:
    - paired tumor vs normal comparisons within long/short groups
    - unpaired comparisons between short vs long
    """
    cat = pd.read_csv('/home_2//analysis/edgeR/total.metadata', sep='\t')
    dat = pd.read_csv('/home_2//analysis/edgeR/total.countdata', sep='\t', index_col=2)
    dat.drop(labels=['ENSEMBL', 'ENTREZID'], axis=1, inplace=True)
    dat = dat.round(1)
    long_t = cat[(cat['pat_type'] == 'long') & (cat['tissue_type'] == 'tumor')]['id'].values
    long_n = cat[(cat['pat_type'] == 'long') & (cat['tissue_type'] == 'normal')]['id'].values
    short_t = cat[(cat['pat_type'] == 'short') & (cat['tissue_type'] == 'tumor')]['id'].values
    short_n = cat[(cat['pat_type'] == 'short') & (cat['tissue_type'] == 'normal')]['id'].values
    def fun(gene):
        return  dat.loc[gene,  long_t].values, dat.loc[gene, long_n].values, dat.loc[gene, short_t].values, dat.loc[gene, short_n].values
    long_t_exp, long_n_exp, short_t_exp,short_n_exp = fun('MUC4')
    stats.ttest_rel(long_t_exp, long_n_exp)
    stats.ttest_rel(short_t_exp, short_n_exp)

    stats.ttest_ind(short_t_exp, long_t_exp)
    stats.ttest_ind(short_n_exp, long_n_exp)



if __name__ == '__main__':
    print ('hello, world')
    RegionLevel()
    RegionLevel1()
    mergeMAFSignature()
    PatientLevel()
    PatientLevelClonality()
