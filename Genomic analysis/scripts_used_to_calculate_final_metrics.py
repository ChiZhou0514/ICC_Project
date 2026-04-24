"""
Metric calculation utilities for cohort-level analyses.

This script aggregates multiple per-sample signals into a single analysis table, e.g.:
- CNV burden (counts and genome instability indices)
- Mutation burden (all vs non-synonymous vs driver mutations)
- Clonality decomposition (clone vs subclone counts, clonality ratio)
- Diversity metrics (Shannon diversity over cluster sizes)
- Tumor purity/ploidy and MSI score

The functions here are typically used to create summary tables and generate figures.
Many inputs are hard-coded to a specific HPC directory structure.
"""

from util import *
import scipy, math
import seaborn as sns
from scipy.spatial import distance_matrix
from skbio.diversity.alpha import shannon
from statannot import add_stat_annotation
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns



### scripts used to calculate metrics



driverGene1 = getDriveGene(); driverGene2 = getDriveInto(); driverGene3 = getDriver()
driverGene = driverGene1 + driverGene2

polyclonal = [25, 28, 33, 38, 58, 63, 8, 82, 9]
NoTreeSample = [15, 27, 44, 37, 51, 34]
NonSyn = ['Missense_Mutation', 'Splice_Region', 'Splice_Site', 'Nonsense_Mutation', 
'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation']
Syn = ['Intron', 'Silent', "3'UTR", "5'UTR", "3'Flank", "5'Flank"]

def getHLA():
    """
    Load HLA LOH annotations and normalize sample barcodes.

    Returns
    -------
    pandas.DataFrame
        HLALOH table with `TumorSampleBarcode` derived from `sample_name`.
    """
    filein = '/home_2//analysisMetAF5/HLALOH.txt'
    dat = pd.read_csv(filein, sep='\t', dtype=object)
    dat['TumorSampleBarcode'] = dat['sample_name'].apply(lambda x: x.replace('_', ''))
    dat.drop(labels=['sample_name'], inplace=True, axis=1)
    return dat

def getCNVcounts(patient, file, cat):
    """
    Count the number of unique genes affected by CNV events (loss or gain).

    Parameters
    ----------
    patient : str
        Patient identifier.
    file : str
        Tissue/sample identifier.
    cat : str
        'Loss' or 'Amp' determining whether to count deletions or amplifications.
    """
    putiry, ploidy = tumorContentDict[patient][file]
    a = []
    filein = '/home_2//outcome/{}/WES/{}/cnv/tumor.callsegmetrics.cns'.format(patient, file)
    dat = pd.read_csv(filein, sep='\t')
    if cat == 'Loss': dat1 = dat[dat['cn'] < ploidy]
    else: dat1 = dat[dat['cn'] > ploidy]
    for i in dat1['gene']:
        a.extend(i.split(','))
    return len(set(a))


def getCNVlen(patient, file, cat):
    """
    Estimate CNV genome fraction (wGII-like) for a given sample.

    Parameters
    ----------
    cat : str
        'Loss', 'Amp', or 'total' controlling which segments are used.
    """
    def fun():
        filein = '/home//database/NCBI/hg38/bwa/Allhs_ref_GRCh38.p7.fa.fai'
        dat = pd.read_csv(filein, sep='\t', header=None)
        mydict = dat[[0,1]].set_index(0).to_dict()[1]
        return mydict
    mydict = fun()
    putiry, ploidy = tumorContentDict[patient][file]
    filein = '/home_2//outcome/{}/WES/{}/cnv/tumor.callsegmetrics.cns'.format(patient, file)
    dat = pd.read_csv(filein, sep='\t')
    if cat == 'Loss':
        dat = dat[ (dat['cn'] < ploidy) & (~dat['chromosome'].isin(['chrX', 'chrY', 'chrM']))]
    elif cat == 'Amp':
        dat = dat[ (dat['cn'] > ploidy) & (~dat['chromosome'].isin(['chrX', 'chrY', 'chrM']))]
    else:
        dat = dat[ (dat['cn'] != ploidy) & (~dat['chromosome'].isin(['chrX', 'chrY', 'chrM']))]
    if dat.shape[0] == 0: return 0
    dat['len'] = dat['end']  - dat['start']
    mydict1 = dat.groupby('chromosome')['len'].sum().to_dict()
    percentage = np.mean([round((mydict1[i] / mydict[i]), 3) for i in mydict1])
    return percentage

def getMutation1(patient, file):
    """
    Return mutation counts for one region:
    - total variants
    - non-synonymous variants
    - non-synonymous driver-gene variants
    """
    maf = '/home_2//outcome/{}/WES/{}/gatk/4pycloneSnpindelsAF5/output_pass.maf'.format(patient, file)
    dat = pd.read_csv(maf, sep='\t')
    dat1 = dat[dat['Variant_Classification'].isin(NonSyn)]
    dat2 = dat1[dat1['Hugo_Symbol'].isin(driverGene)]
    return dat.shape[0], dat1.shape[0], dat2.shape[0]

def getMutation2(patient, file):
    """
    Return clone/subclone mutation counts based on a precomputed clone-state table.
    """
    filein = '/home_2//metastasis/{}/WES/{}/gatk/4pycloneSnpindelsAF5/cloneState9.tsv'.format(patient, file)
    dat = pd.read_csv(filein, sep='\t')
    dat_clone = dat[dat['clone_state'] == 'clone']
    dat_subclone = dat[dat['clone_state'] == 'subclone']
    dat_cloneNonSyn = dat[(dat['clone_state'] == 'clone') & (dat['Variant_Classification'].isin(NonSyn))]
    dat_subcloneNonSyn = dat[(dat['clone_state'] == 'subclone') & (dat['Variant_Classification'].isin(NonSyn))]
    dat1 = dat_clone[(dat_clone['Hugo_Symbol'].isin(driverGene)) & (dat_clone['Variant_Classification'].isin(NonSyn))]
    dat2 = dat_subclone[(dat_subclone['Hugo_Symbol'].isin(driverGene)) & (dat_subclone['Variant_Classification'].isin(NonSyn))]
    return dat_clone.shape[0], dat_subclone.shape[0], dat_cloneNonSyn.shape[0],  dat_subcloneNonSyn.shape[0],  dat1.shape[0], dat2.shape[0]


def getShannoDiv1(patient, file, cat):
    """
    Compute Shannon diversity over cluster sizes within a clonevol-style table.

    Parameters
    ----------
    cat : str
        'subclone', 'clone', or 'all' controlling which clusters are included.
    """
    filein = '/home_2//analysisMetAF5/pycloneMajorCopy/{}/clonevol.txt'.format(patient)
    dat = pd.read_csv(filein, sep='\t')
    dat = dat[['ID', 'cluster', 'gene', 'Driver', file]]
    a = dat.groupby('cluster')[file].mean()
    a = a[a>=min_ccf]
    if 1 not in a.index: return 0
    clone_ccf = a[1]
    a = a / clone_ccf
    if cat == 'subclone': a = a[a < cutoff]
    elif cat == 'clone': a = a[a >= cutoff]
    else: pass
    mydict = dat['cluster'].value_counts().to_dict()
    mylist = [mydict[i] for i in a.index]
    if len(mylist)  <= 0: return 0
    return shannon(mylist)

def getMsiScore(patient, file):
    """Load MSI score from `msi.tsv` generated by msisensor-pro."""
    filein = '/home_2//metastasis/{}/WES/{}/msi/msi.tsv'.format(patient, file)
    dat = pd.read_csv(filein, sep='\t')
    return dat.iloc[0, 2]

def getCloneStatus(patient, file, mydict):
    """
    Lookup helper for external clone-status dictionaries.
    """
    patient = int(patient)
    if patient in mydict:
        if file in mydict[patient]:
            return mydict[patient][file]
    return 'None'

def doAnalysis(dat):
    """
    Enrich an input dataframe with multiple derived metrics in-place.

    Parameters
    ----------
    dat : pandas.DataFrame
        Must contain `patient` and `tissue` columns.
    """
    dat['LossCount'] = [getCNVcounts(patient, file, 'Loss') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['AmpCount'] = [getCNVcounts(patient, file, 'Amp') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['LosswGII'] = [getCNVlen(patient, file, 'Loss') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['AmpwGII'] = [getCNVlen(patient, file, 'Amp') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['wGII'] = [getCNVlen(patient, file, 'total') for patient, file in zip(dat['patient'], dat['tissue'])]

    tmp = [getMutation1(patient, file) for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['AllMutation'] = [i[0] for i in tmp]
    dat['NonSyn'] = [i[1] for i in tmp]
    dat['Driver'] = [i[2] for i in tmp]
    
    tmp = [getMutation2(patient, file) for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['clone'] = [i[0] for i in tmp]
    dat['subclone'] = [i[1] for i in tmp]
    dat['clonality'] = dat['clone'] / (dat['clone'] + dat['subclone'])
    dat['cloneNonSyn'] = [i[2] for i in tmp]
    dat['subcloneNonSyn'] = [i[3] for i in tmp]
    dat['cloneDriver'] = [i[4] for i in tmp]
    dat['subcloneDriver'] = [i[5] for i in tmp]

    dat['subcloneShanno'] = [getShannoDiv1(patient, file, 'subclone') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['cloneShanno'] = [getShannoDiv1(patient, file, 'clone') for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['Shanno'] = [getShannoDiv1(patient, file, 'all') for patient, file in zip(dat['patient'], dat['tissue'])]

    ### purity and ploidy
    mydict1 = f_gettumorContent(isint=False)
    mydict2 = ff_gettumorContent()
    dat['purity'] =  [mydict1[patient][file][0] for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['ploidy'] =  [mydict1[patient][file][1] for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['purityAdjust'] =  [mydict2[patient][file][0] for patient, file in zip(dat['patient'], dat['tissue'])]
    dat['ploidyAdjust'] =  [mydict2[patient][file][1] for patient, file in zip(dat['patient'], dat['tissue'])]
    ### msi
    dat['msi'] =  [getMsiScore(patient, file) for patient, file in zip(dat['patient'], dat['tissue'])]

def vis_certain_index_subgroup(indexName, dat1, axs):
    """
    Draw a boxplot for one metric with statistical annotations on a given axis.
    """
    sns.boxplot(x="typeOfmet6",
            y = indexName,
            #hue = 'Outcome',
            data = dat1,
            order=['PTLNM', 'LNM','PTHM', 'HM'],
            flierprops = {'marker':'o',
                          'markerfacecolor':'lightgrey',
                          'color':'lightgrey',
                         },
            #hue = 'Type',
            width=0.4, 
            #palette=sns.color_palette(['cornflowerblue','coral','forestgreen']),
            ax = axs,
           )

    add_stat_annotation(axs, data=dat1, x='typeOfmet6', y = indexName, 
                            order=['PTLNM', 'LNM','PTHM', 'HM'],
                            box_pairs=[('PTLNM','LNM'), ('PTHM','HM'), ('PTLNM', 'PTHM'), ('LNM', 'HM')],
                            test='Mann-Whitney', text_format='star', loc='inside', verbose=0) #t-test_ind/Mann-Whitney

    axs.set_xlabel("")

def f_vis_certain_index_subgroup(excludeNoTreeSample=False):
    """
    Generate a multi-panel PDF comparing multiple metrics across metastasis subgroups.
    """
    os.chdir('/home_2//analysisMetAF5/maftools/sampleLevel')
    filein = 'sampleLevel.tsv'
    dat = pd.read_csv(filein, sep='\t')
    if excludeNoTreeSample: dat = dat[~dat['patient'].isin(NoTreeSample)]

    indexs = ['LossCount', 'AmpCount', 'LosswGII', 'AmpwGII', 'wGII',
              'AllMutation', 'NonSyn', 'Driver', 'clone', 'subclone', 'cloneNonSyn',
 'subcloneNonSyn', 'cloneDriver', 'subcloneDriver', 'clonality', 'Shanno', 'cloneShanno', 'subcloneShanno',  
 'purity', 'ploidy', 'purityAdjust', 'ploidyAdjust', 'msi']
    
    fig, axslist = plt.subplots(math.ceil(len(indexs) / 4), 4, figsize=(20, math.ceil(len(indexs) / 4) * 5))
    for axs, indexName in zip(np.array(axslist).flatten(), indexs):
        vis_certain_index_subgroup(indexName, dat, axs)
    if excludeNoTreeSample:
        fileout = 'compareIndex_excludeNoTreeSample.pdf'
    else:
        fileout = 'compareIndex.pdf'
    fig.savefig(fileout, bbox_inches='tight')



def fun1LNM(metric, paired = True, excludeNoTreeSample=False):
    """
    Compare PT vs LNM for a metric either paired (within patient) or unpaired.
    """
    dat = pd.read_csv('/home_2//analysisMetAF5/sampleLevel.tsv', sep='\t')
    dat = dat[dat['typeOfmet5'] != 'HM']
    patients = dat[dat['typeOfmet5'] == 'LNM']['patient'].unique()
    dat = dat[dat['patient'].isin(patients)]
    if excludeNoTreeSample: dat = dat[~dat['patient'].isin(polyclonal)]
    if paired:
        a = dat[dat['typeOfmet5'] == 'PT'].groupby('patient')[metric].mean().values
        b = dat[dat['typeOfmet5'] == 'LNM'].groupby('patient')[metric].mean().values
        print (scipy.stats.ttest_rel(a, b)); print (scipy.stats.wilcoxon(a, b))
    else:
        a = dat[dat['typeOfmet5'] == 'PT'][metric].values
        b = dat[dat['typeOfmet5'] == 'LNM'][metric].values
        print (scipy.stats.ttest_ind(a, b)); print (scipy.stats.mannwhitneyu(a, b))
    sns.boxplot(['PT'] * len(a) +  ['LNM'] * len(b), np.append(a, b), meanline=False, showmeans=True, order = ['PT','LNM'])
    return dat, a, b

def fun1M(metric, paired=False, excludeNoTreeSample=False):
    """
    Compare PT vs M (metastasis) for a metric either paired or unpaired.
    """
    dat = pd.read_csv('/home_2//analysisMetAF5/sampleLevel.tsv', sep='\t')
    dat = dat[dat['typeOfmet5'] != 'LNM']
    patients = dat[dat['typeOfmet5'] == 'M']['patient'].unique()
    dat = dat[dat['patient'].isin(patients)]
    if excludeNoTreeSample: dat = dat[~dat['patient'].isin(polyclonal)]
    if paired:
        a = dat[dat['typeOfmet5'] == 'PT'].groupby('patient')[metric].mean().values
        b = dat[dat['typeOfmet5'] == 'M'].groupby('patient')[metric].mean().values
        print (scipy.stats.ttest_rel(a, b)); print (scipy.stats.wilcoxon(a, b))
    else:
        a = dat[dat['typeOfmet5'] == 'PT'][metric].values
        b = dat[dat['typeOfmet5'] == 'M'][metric].values
        print (scipy.stats.ttest_ind(a, b)); print (scipy.stats.mannwhitneyu(a, b))
    sns.boxplot(['PT'] * len(a) +  ['M'] * len(b), np.append(a, b), meanline=False, showmeans=True, order = ['PT','M'])
    return dat, a, b

def ffun1MET():
    dat, a, b = fun1LNM('Driver',  paired= True, excludeNoTreeSample=False)
    dat, a, b = fun1LNM('clonality', paired= True, excludeNoTreeSample=True)

dat = getDataframe(); dat['TumorSampleBarcode'] = dat['patient'] + dat['tissue']
dat['typeOfmet6'] = dat['outcome'].apply(lambda x: 'PTLTS' if x == 'longSur' else 'PTHM')
dat['group'] = 'outcome'
dat = dat[['patient', 'tissue', 'typeOfmet6', 'group']]

cutoff = 0.9
min_ccf = 10
import pickle


with open('/home_2//analysisMetAF5/outcome_tumorContentDict', 'rb') as fout:
    outcome_tumorContentDict = pickle.load(file=fout)    


if __name__ == '__main__':
    doAnalysis(dat)
