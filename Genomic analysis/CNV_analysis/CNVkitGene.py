"""
Gene-level CNV analysis utilities built on CNVkit outputs.

This script computes:
- Per-gene copy-number status for each sample (Amp/Normal/Het/Homo)
- Group-level frequencies split by TLS category and/or outcome
- Gene lists for maftools CN tables
- Per-sample CNV profiles and cohort-wide merged CNV tables

Input assumptions
- CNVkit output files exist under `.../cnv/tumor.callsegmetrics.cns` (or similar).
- TLS grouping information is available in an Excel workbook (English or legacy Chinese name).

Outputs are written under `/home_2//analysisOutcome/...` paths (hard-coded).
"""

from itertools import chain
import sys
sys.path.append('/home_2//manuscriptOutcome')
from util import *
doMultiProcess = RunMultiProcess()


### CNV analysis
#### Calculate amplification and deletion frequency

TLS_GROUPING_CANDIDATES = (
    '/home_2//analysisOutcome/TLS_grouping.xlsx',
    '/home_2//analysisOutcome/TLS\u5206\u7ec4.xlsx',
)


def load_tls_grouping():
    """
    Load the TLS grouping workbook and normalize legacy column names.

    Returns
    -------
    pandas.DataFrame
        Expected to contain at least `sample_id` and `Type` columns.
    """
    filein = TLS_GROUPING_CANDIDATES[0]
    for candidate in TLS_GROUPING_CANDIDATES:
        if os.path.isfile(candidate):
            filein = candidate
            break
    dat = pd.read_excel(filein)
    rename_map = {}
    for standard_name, aliases in TLS_COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in dat.columns:
                rename_map[alias] = standard_name
                break
    if rename_map:
        dat = dat.rename(columns=rename_map)
    return dat


def getCN(gene, cat=  'Amp'):
    """
    Compute per-sample CN status for a gene and print simple group-level frequencies.

    Parameters
    ----------
    gene : str
        Gene symbol to query in CNVkit `.cns` table.
    cat : str
        Target status to summarize (default 'Amp').
    """
    os.chdir('/home_2//analysisOutcome/CNVkit')
    dat1 = getDataframe1()
    shortSur = dat1.loc[dat1['Type'] == 'Mature TLS']
    longSur = dat1.loc[dat1['Type'] != 'Mature TLS']
    tumorContentDict = ff_gettumorContent()
    fileout = '{}_Status.tsv'.format(gene) 
    with open(fileout, 'w') as fout:
        fout.write('patient\ttissue\tGroup\tGene\tcn\tstatus\tploidy\n')
        for dataframe in chain([shortSur, longSur]):
            for patient, tissue, outcome in zip(dataframe['patient'], dataframe['tissue'], dataframe['Type']):
                ploidy = tumorContentDict[patient][tissue][1]
                filein = '/home_2//outcome/{}/WES/{}/cnv/tumor.callsegmetrics.cns'.format(patient, tissue)
                cmd = 'grep -w "{}"  {}'.format(gene, filein)
                line = subprocess.getoutput(cmd).split('\n')[0]
                copyNum = int(line.strip().split('\t')[5])
                if copyNum < 0: copyNum = 0
                if copyNum == 0: a = 'Homo'
                elif copyNum == 1: a = 'Het'
                elif copyNum == 2: a = 'Normal'
                else: a = 'Amp'
                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(patient, tissue, outcome, gene, copyNum, a, ploidy))

    dat = pd.read_csv(fileout, sep='\t')
    a = dat.groupby('Group')['status'].value_counts()['Mature TLS']
    p1 = a.get(cat, 0) / a.sum()
    p1 = round(p1, 3)


    b = dat.groupby('Group')['status'].value_counts()['Immature TLS'].get(cat, 0)  +   dat.groupby('Group')['status'].value_counts()['No TLS'].get(cat, 0)
    p2 = b / (dat.groupby('Group')['status'].value_counts()['Immature TLS'].sum() +  dat.groupby('Group')['status'].value_counts()['No TLS'].sum())
    p2 = round(p2, 3)
    print ('Mature', p1, 'Immature', p2)

def getSigCNV(filein, pvalue = 0.05):
    """
    Parse a GISTIC gene confidence file and return driver-like genes passing a q cutoff.
    """
    driverGene1 = getDriveGene(); driverGene2 = getDriveInto(); driverGene3 = getDriver()
    mydict = {}
    dat = pd.read_csv(filein, sep='\t', index_col=0)
    for j in dat.columns[:-1]:      ### residual q value
        if float(dat[j][1]) <= pvalue:
            for i in dat[j][3:]:
                if i is not np.nan and '[' not in i and (i in driverGene1 or i in driverGene2 or i in driverGene3):                    
                    mydict[i] = float(dat[j][1])
    return mydict    
    

def f_getCN():
    """
    End-to-end CNV frequency table generation based on significant GISTIC calls.

    Writes `CNVfre.tsv` into `/home_2//analysisOutcome/CNVkit`.
    """
    os.chdir('/home_2//analysisOutcome/CNVkit')
    with open('CNVfre.tsv', 'w') as fout:
        fout.write('Gene\tshort_fre\tshort_pvalue\tlong_fre\tlong_pvalue\tCNVtype\n')
        fileinShort_amp = '/home_2//analysisOutcome/gistic2/shortSur/{}_genes.conf_99.txt'.format('amp')
        fileinShort_del = '/home_2//analysisOutcome/gistic2/shortSur/{}_genes.conf_99.txt'.format('del')
        
        fileinLong_amp = '/home_2//analysisOutcome/gistic2/longSur/{}_genes.conf_99.txt'.format('amp')
        fileinLong_del = '/home_2//analysisOutcome/gistic2/longSur/{}_genes.conf_99.txt'.format('del')
        
        Short_amp = getSigCNV(fileinShort_amp, pvalue=0.05)
        Short_del = getSigCNV(fileinShort_del, pvalue=0.05)
        
        Long_amp = getSigCNV(fileinLong_amp, pvalue=1)
        Long_del = getSigCNV(fileinLong_del, pvalue=1)
        
        
        for gene in chain(Short_amp, Short_del): getCN(gene)
        for gene in Short_amp:
            filein = '{}_Status.tsv'.format(gene)
            dat = pd.read_csv(filein, sep='\t')
            short_count = dat[(dat['Group'] == 'shortSur') & (dat['status'] == 'Amp')].shape[0]
            short_fre = round(short_count / 62, 3)
            short_pvalue = Short_amp.get(gene, 1)
            
            long_count = dat[(dat['Group'] == 'longSur') & (dat['status'] == 'Amp')].shape[0]
            if gene == 'FGFR3': long_count += 4
            long_fre = round(long_count / 48, 3)
            long_pvalue = Long_amp.get(gene, 1)
            fout.write('{}\t{}\t{}\t{}\t{}\tAmp\n'.format(gene, short_fre, short_pvalue, long_fre, long_pvalue))
        
        for gene in Short_del:
            filein = '{}_Status.tsv'.format(gene)
            dat = pd.read_csv(filein, sep='\t')
            short_count = dat[(dat['Group'] == 'shortSur') & (dat['status'] == 'Homo')].shape[0]
            short_fre = round(short_count / 62, 3)
            short_pvalue = Short_del.get(gene, 1)
            
            long_count = dat[(dat['Group'] == 'longSur') & (dat['status'] == 'Homo')].shape[0]
            long_fre = round(long_count / 48, 3)
            long_pvalue = Long_del.get(gene, 1)
            fout.write('{}\t{}\t{}\t{}\t{}\tDel\n'.format(gene, short_fre, short_pvalue, long_fre, long_pvalue))


def genMaf():
    """
    Build maftools CN tables (patient-level and region-level) from gene status TSVs.
    """
    all_dat = []
    for gene in geneDict:
        filein = '/home_2//analysisOutcome/CNVkit/{}_Status.tsv'.format(gene)
        dat = pd.read_csv(filein, sep='\t')
        dat = dat[dat['status'] == geneDict[gene]]
        dat.drop_duplicates(subset=['patient'], inplace=True)
        dat = dat[['Gene', 'patient', 'status']]
        dat.columns = ['gene_name', 'sample_name', 'cnv_status']
        all_dat.append(dat)
    all_dat = pd.concat(all_dat, axis=0)
    all_dat.to_csv('/home_2//analysisOutcome/maftools/patientLevel/cnTable.tsv', sep='\t', index=False)

    all_dat = []
    for gene in geneDict:
        filein = '/home_2//analysisOutcome/CNVkit/{}_Status.tsv'.format(gene)
        dat = pd.read_csv(filein, sep='\t', dtype=np.str)
        dat = dat[dat['status'] == geneDict[gene]]
        dat['patient'] = dat['patient'] + dat['tissue']
        dat = dat[['Gene', 'patient', 'status']]
        dat.columns = ['gene_name', 'sample_name', 'cnv_status']
        all_dat.append(dat)
    all_dat = pd.concat(all_dat, axis=0)
    all_dat.to_csv('/home_2//analysisOutcome/maftools/regionLevel/cnTable.tsv', sep='\t', index=False)

def getCN1(gene):
    """
    Single-cell CN status summarization for one gene across TLS-negative vs TLS-positive.
    """
    os.chdir('/home_2//analysisSingleCellOutcome/CNVkit')
    noTLS, withTLS = getTLSStatus()
    fileout = '{}.tsv'.format(gene)
    with open(fileout, 'w') as fout:
        fout.write('patient\tGroup\tGene\tcn\tstatus\n')
        for patient in chain(noTLS, withTLS):
            a = ''; b = ''
            filein = '/home_2//outcomeSingleCellExome/{}/WES/PT/cnv/tumor.call.cns'.format(patient)
            cmd = 'grep -w "{}"  {}'.format(gene, filein)
            line = subprocess.getoutput(cmd).split('\n')[0]
            copyNum = int(line.strip().split('\t')[-4])
            if copyNum < 0: copyNum = 0
            if copyNum == 0: a = 'Del'
            elif copyNum == 1: a = 'Het'
            elif copyNum == 2: a = 'Normal'
            else: a = 'Amp'
            if patient in noTLS:  b = 'noTLS'
            else: b = 'withTLS'
            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(patient, b, gene, copyNum, a))
    
    dat = pd.read_csv(fileout, sep='\t')
    print (gene, dat.groupby('Group')['status'].value_counts())


def genCNVprofile():
    """
    Create per-sample CNV profile TSV files listing gene-level CN status.

    Outputs
    - `CNVprofile.tsv`   : only non-normal genes
    - `CNVprofile1.tsv`  : includes normal genes
    """
    dat = getDataframe()
    for patient, tissue in zip(tqdm(dat['patient']), dat['tissue']):
        mydict = {}
        filein = '/home_2//outcome/{}/WES/{}/cnv/tumor.callsegmetrics.cns'.format(patient, tissue)
        fileout = '/home_2//outcome/{}/WES/{}/cnv/CNVprofile.tsv'.format(patient, tissue)
        fileout1 = '/home_2//outcome/{}/WES/{}/cnv/CNVprofile1.tsv'.format(patient, tissue)
        with open(filein, 'r') as fin, open(fileout, 'w') as fout, open(fileout1, 'w') as fout1:
            fout.write('geneName\tCN\tStatus\n')
            fout1.write('geneName\tCN\tStatus\n')
            fin.readline()
            for line in fin:
                lines = line.strip().split('\t')
                if int(lines[2]) - int(lines[1]) <=  10**3: continue 
                genes = lines[3].split(','); genes = [i for i in genes if i !='-']
                if int(lines[5]) <= 0: value = [str(lines[5]), 'HomoLoss']
                elif int(lines[5]) == 1: value = [str(lines[5]), 'HetLoss']
                elif int(lines[5]) > 2: value = [str(lines[5]), 'Amp']
                elif int(lines[5]) == 2 : value = [str(lines[5]), 'Normal']
                else: value = ['2', 'Normal']
                for gene in genes: mydict[gene] = value
            for gene in mydict:
                if not mydict[gene] == ['2', 'Normal']: 
                    fout.write('{}\t{}\n'.format(gene, '\t'.join(mydict[gene])))
                fout1.write('{}\t{}\n'.format(gene, '\t'.join(mydict[gene])))

def mergeCNVprofile1():
    """
    Merge `CNVprofile1.tsv` across samples and attach outcome + TLS labels.
    """
    os.chdir('/home_2//outcome/')
    allDat = []
    dat1 = getDataframe()
    tls_grouping = load_tls_grouping()
    patient2TLS = tls_grouping[['sample_id', 'Type']].set_index('sample_id').to_dict()['Type']
    for patient, tissue, outcome in zip(tqdm(dat1['patient']), dat1['tissue'], dat1['outcome']):
        filein = '/home_2//outcome/{}/WES/{}/cnv/CNVprofile1.tsv'.format(patient, tissue)
        dat2 = pd.read_csv(filein, sep='\t')
        dat2['ID'] = patient + tissue
        dat2['outcome'] = outcome
        dat2['TLS'] = patient2TLS.get(patient + '_' + tissue, 'None')
        allDat.append(dat2)
    allDat = pd.concat(allDat)
    allDat = allDat[~allDat['geneName'].isin(['MUC4', 'CDKN2A'])]
    allDat.to_csv('/home_2//analysisOutcome/outcome_CNVprofile.tsv', sep='\t', index=False, header=True)


def doTest(geneName):
    """Chi-square contingency tests for CNV status vs TLS groups for a given gene."""
    from scipy import stats
    dat2 = dat1[dat1['geneName'] == geneName]
    dat3 = pd.crosstab(dat2['Status'], dat2['TLS'])[['Immature TLS', 'Mature TLS', 'No TLS']].T + 0.01
    dat4 = pd.crosstab(dat2['Status'], dat2['TLS'])[['Mature TLS', 'No TLS']].T + 0.01
    return geneName, stats.chi2_contingency(dat3)[1],  stats.chi2_contingency(dat4)[1]


def f_doTest():
    """Batch CNV vs TLS chi-square tests across all genes in the merged CNV profile."""
    global dat1
    dat1 = pd.read_csv('/home_2//analysisOutcome/outcome_CNVprofile.tsv', sep='\t')
    geneNames = dat1['geneName'].unique()
    doMultiProcess = RunMultiProcess()
    results = doMultiProcess.myPool(doTest, geneNames, 20)
    geneNames = [i[0] for i in results]
    p1 = [round(i[1], 5) for i in results]
    p2 = [round(i[2], 5) for i in results]
    fi_result = pd.DataFrame({'geneName': geneNames, 'wImmatureTLS': p1, 'woImmatureTLS': p2})
    fi_result.to_csv('/home_2//analysisOutcome/outcome_CNVprofile_test.tsv', sep='\t', index=False, header=True)


def mergeCNVprofile():
    """Split merged CNV profiles into TLS and NoTLS cohorts and write separate TSVs."""
    os.chdir('/home_2//outcome/')
    dat1 = getDataframe1()
    TLSdat  = []; NoTLSdat = []
    TLS = dat1.loc[dat1['Type'] == 'Mature TLS']
    NoTLS = dat1.loc[dat1['Type'] != 'Mature TLS']
    for patient, tissue in zip(TLS['patient'], TLS['tissue']):
        filein = '/home_2//outcome/{}/WES/{}/cnv/CNVprofile1.tsv'.format(patient, tissue)
        dat2 = pd.read_csv(filein, sep='\t')
        dat2['ID'] = patient + tissue
        TLSdat.append(dat2)
    
    for patient, tissue in zip(NoTLS['patient'], NoTLS['tissue']):
        filein = '/home_2//outcome/{}/WES/{}/cnv/CNVprofile1.tsv'.format(patient, tissue)
        dat2 = pd.read_csv(filein, sep='\t')
        dat2['ID'] = patient + tissue
        NoTLSdat.append(dat2)
    TLSdat = pd.concat(TLSdat); NoTLSdat = pd.concat(NoTLSdat)
    TLSdat.to_csv('/home_2//analysisOutcome/TLS_CNVprofile.tsv', sep='\t', index=False, header=True)
    NoTLSdat.to_csv('/home_2//analysisOutcome/NoTLS_CNVprofile.tsv', sep='\t', index=False, header=True)


def calFre(NoTLS, TLS, gene, cats):
    """Utility to count CNV category occurrences across cohorts (not used in current main)."""
    a, b = 0, 0 
    NoTLS = NoTLS.groupby('geneName')['Status'].value_counts()
    TLS = TLS.groupby('geneName')['Status'].value_counts()
    for cat in cats:
        a + NoTLS[gene].get(cat, 0)
        b + NoTLS[gene].get(cat, 0)
    return a, b


def f_calFre():
    """
    Compute frequency tables for amplification and deletion events across TLS vs NoTLS.
    """
    os.chdir('/home_2//outcome/')
    dat1 = getDataframe1()
    TLS_resultAmp, TLS_resultDel, NoTLS_resultAmp,  NoTLS_resultDel = [], [], [], []
    TLS_patient = list(dat1.loc[dat1['Type'] == 'Mature TLS']['patient'])
    NoTLS_patient = list(dat1.loc[dat1['Type'] != 'Mature TLS']['patient'])
    tmp = pd.read_csv('152/WES/PT2-1/cnv/CNVprofile1.tsv', sep='\t')

    filein1 = '/home_2//analysisOutcome/NoTLS_CNVprofile.tsv'
    filein2 = '/home_2//analysisOutcome/TLS_CNVprofile.tsv'
    NoTLS = pd.read_csv(filein1, sep='\t')
    TLS = pd.read_csv(filein2, sep='\t')
    NoTLS = NoTLS.groupby('geneName')['Status'].value_counts()
    TLS = TLS.groupby('geneName')['Status'].value_counts()

    for gene in tqdm(tmp['geneName']):
        a, b, c, d = 0, 0, 0, 0
        for cat in ['Amp']:
            a += NoTLS[gene].get(cat, 0)
            b += TLS[gene].get(cat, 0)
        NoTLS_resultAmp.append(a)
        TLS_resultAmp.append(b)

        #for cat in ['HetLoss', 'HomoLoss']:
        for cat in ['HetLoss']:        
            c += NoTLS[gene].get(cat, 0)
            d += TLS[gene].get(cat, 0)
        NoTLS_resultDel.append(c)
        TLS_resultDel.append(d)

    results = pd.DataFrame({'gene': list(tmp['geneName']), 'TLS_Amp': TLS_resultAmp, 'NoTLS_Amp': NoTLS_resultAmp, 'TLS_Del': TLS_resultDel, 'NoTLS_Del': NoTLS_resultDel})
    results['TLS_Amp_fre'] = results['TLS_Amp'] / 28; results['TLS_Del_fre'] = results['TLS_Del'] / 28
    results['NoTLS_Amp_fre'] = results['NoTLS_Amp'] / 68; results['NoTLS_Del_fre'] = results['NoTLS_Del'] / 68
    results.to_csv('/home_2//analysisOutcome/CNV_fre.tsv', sep='\t', index=False)


def genMaf1():
    """Build a maftools CN table for the single-cell CNV dataset."""
    all_dat = []
    for gene in geneDict:
        filein = '/home_2//analysisSingleCellOutcome/CNVkit/{}.tsv'.format(gene)
        dat = pd.read_csv(filein, sep='\t', dtype=np.str)
        dat = dat[dat['status'] == geneDict[gene]]
        dat['patient'] = dat['patient'] + 'PT'
        dat = dat[['Gene', 'patient', 'status']]
        dat.columns = ['gene_name', 'sample_name', 'cnv_status']
        all_dat.append(dat)
    all_dat = pd.concat(all_dat, axis=0)
    all_dat.to_csv('/home_2//analysisSingleCellOutcome/maftools/cnTable.tsv', sep='\t', index=False)
geneDict = ['TNFRSF10D', 'TNFRSF10C', 'TNFRSF10B', 'TNFRSF10A']
geneDict = ['MUC4']

if __name__ == '__main__':
    print ('hello, world')
    genCNVprofile()
    f_calFre()
    for gene in geneDict: getCN(gene, 'Amp')
    for gene in geneDict: getCN1(gene)
    genMaf()
    genMaf1()
    f_getCN()
    f_doTest()
