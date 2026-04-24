"""
Pathway-level mutation enrichment analysis (Hallmark gene sets).

This script follows the logic described in:
Unraveling tumor–immune heterogeneity in advanced ovarian cancer uncovers immunogenic
effect of chemotherapy

Workflow overview
1) Load Hallmark gene sets (MSigDB-style table stored locally)
2) For each sample, count non-synonymous mutations that fall into each Hallmark gene set
3) Write count tables (all mutations, clone-only, subclone-only)
4) Run a likelihood-ratio test via multinomial logistic regression (MNLogit) to assess
   whether Hallmark counts differ across labels (e.g., TLS group or outcome), optionally
   controlling for patient effects in multi-region samples.

Notes
- This is a project-specific analysis with hard-coded file paths.
- Column naming is normalized to support both English and legacy Chinese names.
"""

from util import *
from itertools import chain
from collections import Counter
import seaborn as sns
import scipy

NonSyn = ['Missense_Mutation', 'Splice_Region', 'Splice_Site', 'Nonsense_Mutation', 
'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']

SUMMARY_COLUMNS = ['sample_id_level1', 'Type', 'Type-50%', 'sample_id', 'patient', 'tissue', 'outcome', 'Multi', 'NonSyn']


def normalize_analysis_columns(dataframe):
    rename_map = {}
    for standard_name, aliases in TLS_COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in dataframe.columns:
                rename_map[alias] = standard_name
                break
    if rename_map:
        dataframe = dataframe.rename(columns=rename_map)
    return dataframe


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def getHallmark():
    """
    Load Hallmark gene sets from a local tab-delimited file.

    Returns
    -------
    dict[str, list[str]]
        Mapping of hallmark name (lower-cased) to list of gene symbols.
    """
    mydict = {}
    dat  = pd.read_csv('/home_2//HallMark.txt', sep='\t')
    for i, j in zip(dat['GeneSet'], dat['Genes']):
        genes = j.strip().split(',')
        mydict[i.lower()] = genes
    return mydict


def f_getHallmark():
    """
    Generate a per-sample Hallmark mutation count table from region-level MAFs.

    Output
    - `/home_2//analysisOutcome/PathwayEnrichment/HallMarkCounts.tsv`
    """
    dat = pd.read_csv('/home_2//analysisOutcome/Analysis1.tsv',sep='\t')
    dat = normalize_analysis_columns(dat)
    mydict = getHallmark()
    for Hallmark in mydict:
        temp = []
        for patient, tissue in zip(dat['patient'], dat['tissue']):
            maf = '/home_2//outcome/{}/WES/{}/gatk/lifted_overPASS.maf'.format(patient, tissue)
            dat1 = pd.read_csv(maf, sep='\t')
            dat2 = dat1[dat1['Variant_Classification'].isin(NonSyn)]
            temp.append(sum(dat2['Hugo_Symbol'].isin(mydict[Hallmark])))
        dat[Hallmark] = temp
    columns = SUMMARY_COLUMNS + list(mydict.keys())
    dat = dat[columns]
    dat.to_csv('/home_2//analysisOutcome/PathwayEnrichment/HallMarkCounts.tsv', sep='\t', index= False, float_format='%.3f')

def f_getHallmark1(clone_state):
    """
    Generate clone-only or subclone-only Hallmark count tables.

    Parameters
    ----------
    clone_state : str
        Either 'clone' or 'subclone' as encoded in `clone_state9.tsv`.

    Outputs
    - HallMarkCloneCounts.tsv (clone)
    - HallMarkSubCloneCounts.tsv (subclone)
    """
    dat = pd.read_csv('/home_2//analysisOutcome/Analysis1.tsv',sep='\t')
    dat = normalize_analysis_columns(dat)
    mydict = getHallmark()
    for Hallmark in mydict:
        temp = []
        for patient, tissue in zip(dat['patient'], dat['tissue']):
            maf = '/home_2//outcome/{}/WES/{}/gatk/4pyclone_snpindels/clone_state9.tsv'.format(patient, tissue)
            dat1 = pd.read_csv(maf, sep='\t')
            dat2 = dat1[dat1['Variant_Classification'].isin(NonSyn)]
            dat2 = dat2[dat2['clone_state'] == clone_state]
            temp.append(sum(dat2['Hugo_Symbol'].isin(mydict[Hallmark])))
        dat[Hallmark] = temp
    columns = SUMMARY_COLUMNS + list(mydict.keys())
    dat = dat[columns]
    if clone_state == 'clone':
        dat.to_csv('/home_2//analysisOutcome/PathwayEnrichment/HallMarkCloneCounts.tsv', sep='\t', index= False, float_format='%.3f')
    else:
        dat.to_csv('/home_2//analysisOutcome/PathwayEnrichment/HallMarkSubCloneCounts.tsv', sep='\t', index= False, float_format='%.3f')



def doAnalysis(label, Hallmark, multi=True):
    """
    Likelihood-ratio test for association between a Hallmark count and a categorical label.

    The full model includes the Hallmark count plus covariates:
    - `NonSyn` (overall non-synonymous mutation count)
    - `patient1` (encoded patient ID) when `multi=True`

    The partial model excludes the Hallmark count. The LRT p-value is returned.
    """
    try:
        from sklearn.preprocessing import LabelEncoder
        import statsmodels.api as sm
        dat = pd.read_csv(filein,sep='\t')
        if multi:
            dat = dat[dat['Multi']]
            a = ['NonSyn', 'patient1']
        else:
            a = ['NonSyn']
        le = LabelEncoder()
        y = le.fit_transform(dat[label])
        dat['patient1'] = le.fit_transform(dat['patient'])
        dat1 = dat[[Hallmark] + a]   ### full model
        dat1 = sm.add_constant(dat1)
        dat2 = dat[a]             ### partial model
        dat2 = sm.add_constant(dat2)
        full_model = sm.MNLogit(y, dat1).fit().llf
        partial_model = sm.MNLogit(y,dat2).fit().llf
        LR_statistic = -2*(partial_model - full_model)
        p_val = scipy.stats.chi2.sf(LR_statistic, len(a))
    except:
        p_val = 1
    return p_val

def f_doAnalysis(label, multi):
    """
    Run enrichment analysis across all Hallmarks and write an output TSV.

    Also prints significant Hallmarks (FDR <= 0.05) with group means.
    """
    dat = pd.read_csv(filein,sep='\t')
    mydict  =  getHallmark()
    Hallmarks = list(mydict.keys())
    p_vals = []
    for Hallmark in Hallmarks:
        p_val = doAnalysis(label, Hallmark, multi)
        p_vals.append(p_val)
    FDR = p_adjust_bh(p_vals)
    for Hallmark, fdr in zip(Hallmarks, FDR):
        if fdr <=0.05:
            a = dat.groupby(label)[Hallmark].mean()
            if label == 'outcome':
                print ('fdr', '\t', round(fdr, 3), '\t' ,Hallmark, '\t', round(a['longSur'], 3), '\t', round(a['shortSur'], 3))
            else:
                print ('fdr', '\t', round(fdr, 3), '\t' ,Hallmark, '\t', round(a['Mature TLS'], 3), '\t', round(a['Immature TLS'], 3), '\t', round(a['No TLS'], 3))
            

    with open(fileout, 'w') as fout:
        if label == 'outcome':
            fout.write('{}\n'.format('\t'.join(['pvalue', 'HallMark', 'longSur', 'shortSur'])))
        else:
            fout.write('{}\n'.format('\t'.join(['pvalue', 'HallMark', 'Mature TLS', 'Immature TLS', 'No TLS'])))
        
        for Hallmark, p_val in zip(Hallmarks, p_vals):
            if p_val <=1:
                a = dat.groupby(label)[Hallmark].mean()
                if label == 'outcome':
                    a1 = str(round(a['longSur'], 3)); a2 = str(round(a['shortSur'], 3));
                    fout.write('{}\n'.format('\t'.join([str(round(p_val, 3)), Hallmark, a1, a2])))
                else:
                    a1 = str(round(a['Mature TLS'], 3));a2 = str(round(a['Immature TLS'], 3)); a3 = str(round(a['No TLS'], 3))
                    fout.write('{}\n'.format('\t'.join([str(round(p_val, 3)), Hallmark, a1, a2, a3])))


   
    
label = 'Type'  ### outcome, Type
clonestat = 'SubClone'  #### '', 'Clone',  'SubClone'
filein = '/NFS_home/NFS_home_3///analysisOutcome/PathwayEnrichment/HallMark{}Counts.tsv'.format(clonestat)

if label == 'outcome':
    fileout = '/home_2//analysisOutcome/PathwayEnrichment/outcome{}_pathwayEnrichment.tsv'.format(clonestat)
else:
    fileout = '/home_2//analysisOutcome/PathwayEnrichment/TLS{}_pathwayEnrichment.tsv'.format(clonestat)

if __name__ == '__main__':
    f_getHallmark()
    f_getHallmark1('clone')
    f_getHallmark1('subclone')
    f_doAnalysis(label, True)
