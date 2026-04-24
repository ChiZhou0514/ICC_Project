"""
CNVkit wrapper for calling copy-number variants (CNVs) in WES data.

This script orchestrates CNVkit commands to:
- build target/antitarget bins (preProcess)
- create a normal reference panel (creatNormalRef)
- call CNVs for each tumor sample (CNVkit)
- optionally run ABSOLUTE on derived segment files (Absolute)
- export inputs for GISTIC (merge4gistic)
- generate per-sample CNV heatmaps (CNVHeatMap / ff_CNVHeatMap)

Assumptions
- Designed for an HPC environment with hard-coded paths under `/home_2//`.
- Uses purity/ploidy from Sequenza (`ff_gettumorContent` in util).
- Uses `f_subprocess` to log command output into `log.txt`.
"""

from util import *
import seaborn as sns


doMultiProcess = RunMultiProcess()
os.chdir('/home_2/')
def f_subprocess(cmd):
    """
    Run command segments joined by '&&' and append output to `log.txt`.
    """
    with open('log.txt', 'a') as fout:
        p = subprocess.Popen('&&'.join(cmd), shell=True, stderr=subprocess.STDOUT, stdout=fout)
        p.wait()

def preProcess():
    """
    Create CNVkit target and anti-target bins.

    Outputs (in the current CNVkit work directory)
    - `targets.bed` plus derived antitarget/target bin sizes via `cnvkit.py autobin`.
    """
    os.chdir('/home_2//CNVkitFile/outcome')
    cmd = 'cnvkit.py  target  {bed} --annotate  {refFlat}  --split -o targets.bed'.format(**doMultiProcess.__dict__) #### generate targets.bed
    subprocess.call(cmd, shell=True)
    ### automatically generate anti-target and target bins with cnvkit autobin
    cmd = 'cnvkit.py  autobin  /home_2//outcomeSingleCellExome/*/WES/*/gatk/deduped.bam  -m  hybrid  -g  {mappable_file}  --annotate {refFlat}  -t  targets.bed'.format(**doMultiProcess.__dict__)
    subprocess.call(cmd, shell=True)

def creatNormalRef():
    """Create a CNVkit reference profile from normal (NC) coverage files."""
    os.chdir('/home_2//outcomeSingleCellExome')
    cmd = 'cnvkit.py reference */WES/NC/cnv/*coverage.cnn -f {fasta} -o /home_2//CNVkitFile/outcomeSingleCell/Reference.cnn'.format(**doMultiProcess.__dict__)
    subprocess.call(cmd, shell=True)

 
def CNVkit(dirName):
    """
    Call CNVs for one tumor tissue directory using CNVkit.

    Parameters
    ----------
    dirName : str
        Tissue directory containing `gatk/deduped.bam` and CNV output folders.
    """
    os.chdir(dirName)
    targets_bed = '/home_2//CNVkitFile/outcomeMetastasis/targets.target.bed'
    antitargets_bed = '/home_2//CNVkitFile/outcomeMetastasis/targets.antitarget.bed'
    ref_cnn = '/home_2//CNVkitFile/outcomeMetastasis/Reference.cnn'
    patient = os.path.basename(os.path.dirname(os.path.dirname(dirName)))
    tissue = os.path.basename(dirName)
    caseid = patient + tissue
    print (caseid)
    if not os.path.isdir('cnv1'): os.makedirs('cnv1')

    cmd = []; cmd1 = cmd2 = cmd3 = cmd4 = ''
    cmd1 = 'cnvkit.py coverage  gatk/deduped.bam  {} -o cnv1/targetcoverage.cnn -p 15 -q 20'.format(targets_bed)
    cmd2 = 'cnvkit.py coverage  gatk/deduped.bam  {} -o cnv1/antitargetcoverage.cnn -p 15 -q 20'.format(antitargets_bed)
    #f_subprocess([cmd1, cmd2])
    cmd3 = r'cnvkit.py fix cnv1/targetcoverage.cnn  cnv1/antitargetcoverage.cnn  {} -o cnv1/tumor.cnr -i {}'.format(ref_cnn, caseid)
    cmd4 = r'cnvkit.py segment cnv1/tumor.cnr -o  cnv1/tumor.cns -p 30 --drop-low-coverage  --drop-outliers 10'
    f_subprocess([cmd3, cmd4])
    
    purity, ploidy = tumorContentDict[patient][tissue]
    cmd5 = r'cnvkit.py call  cnv1/tumor.cns  -m  clonal  --min-variant-depth 50  -v ../NC/germline/germlinePASS.vcf  -o  cnv1/tumor.callGermline.cns  --purity {:.2f}  --ploidy {:.0f}  --drop-low-coverage -i {} -n {}NC '.format(purity, ploidy, caseid, patient)
    cmd6 = r'/home//anaconda3/bin/cnvkit.py  segmetrics    cnv1/tumor.cnr  -s cnv1/tumor.cns  --ci   --pi  --sem  -o  cnv1/tumor.segmetrics.cns'
    f_subprocess([cmd6])

    ### call
    cmd7 = r'/home//anaconda3/bin/cnvkit.py call  cnv/tumor.cns  -m  clonal  --min-variant-depth 50  -o  cnv/tumor.call.cns  --purity {:.2f}  --ploidy {:.0f}  --drop-low-coverage'.format(purity, ploidy) 
    cmd8 = r'/home//anaconda3/bin/cnvkit.py call  cnv1/tumor.segmetrics.cns  --filter  ci  -m  clonal  --min-variant-depth 50  -o  cnv1/tumor.callsegmetrics.cns  --purity {:.2f}  --ploidy {:.0f}  --drop-low-coverage'.format(purity, ploidy) ### segmetrics has limited impact on downstream results
    cmd9 = r'cnvkit.py call  cnv/tumor.cns   --min-variant-depth 50  -o  cnv/tumor.callthreshold.cns  --drop-low-coverage'.format(purity, ploidy)  ### keep an uncorrected seg_mean because clonal mode strongly adjusts seg_mean by purity and ploidy
    subprocess.call(cmd8, shell=True)

    cmd1 = r'cnvkit.py  export seg  cnv/tumor.call.cns  -o  cnv/tumor.call.seg'
    cmd2 = r'cnvkit.py  export seg  cnv1/tumor.callsegmetrics.cns  -o  cnv1/tumor.callsegmetrics.seg'
    cmd3 = r'cnvkit.py  export seg  cnv/tumor.callthreshold.cns  -o  cnv/tumor.callthreshold.seg'
    subprocess.call(cmd2, shell=True)
    def fun1(filein, fileout):
        dat = pd.read_csv(filein, sep='\t')
        dat.columns = ['Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean']
        dat.to_csv(fileout, sep='\t', header=True, index=False)

    #fun1('cnv/tumor.call.seg', 'cnv/plotCBSsegments.call.seg')
    #fun1('cnv/tumor.callthreshold.seg', 'cnv/plotCBSsegments.callthreshold.seg')
    #fun1('cnv/tumor.callsegmetrics.seg', 'cnv/plotCBSsegments.callsegmetrics.seg')

    # dat = pd.read_csv('cnv/tumor.call.seg', sep='\t')
    # dat.columns = ['ID', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean']
    # dat['Chromosome'] = dat['Chromosome'].apply(lambda x: x.strip('chr'))
    # dat = dat[dat['End'] - dat['Start'] >= 500000]
    # dat.to_csv('cnv/tumor2absolute.seg', index=False, sep='\t')


def f_CNVkit():
    """Batch CNV calling for all tumor tissues under `/home_2//outcome/*/WES/*`."""
    dirNames = glob.glob('/home_2//outcome/*/WES/*')
    dirNames = [i for i in dirNames if os.path.basename(i) != 'NC']
    doMultiProcess.myPool(CNVkit, dirNames, processes=10)


def Absolute(dirName):
    """Run ABSOLUTE on a tissue directory (removes existing `absolute/` results first)."""
    os.chdir(dirName)
    if os.path.isdir('absolute'): shutil.rmtree('absolute/')
    cmd = '/usr/local/bin/Rscript   /home_2//manuscriptOutcome/absolute.r'
    f_subprocess([cmd])

def f_Absolute():
    """Batch ABSOLUTE runs for all non-NC samples in the single-cell exome dataset."""
    #dirNames = glob.glob('/home_2//outcome/*/WES/*')
    dirNames = glob.glob('/home_2//outcomeSingleCellExome/*/WES/*')
    dirNames = [i for i in dirNames if os.path.basename(i) != 'NC']
    doMultiProcess.myPool(Absolute, dirNames, processes=5)


###  /home_2//analysisOutcome/gistic2
#cmd = 'cp  /home//software/GISTIC2/gistic2 /home//software/GISTIC2/gp_gistic2_from_seg .'
#subprocess.call(cmd, shell=True)
def merge4gistic(cat = 'longSur'):
    """
    Merge CNVkit segment files into a single `.seg` file for GISTIC.

    Parameters
    ----------
    cat : str
        'shortSur', 'longSur', or 'all' controlling which samples are merged.
    """
    alldat = []
    dat = getDataframe()
    if cat != 'all':
        dat = dat[dat['outcome'] == cat]
    for patient, tissue in zip(dat['patient'], dat['tissue']):
        filein = '/home_2//outcome/{}/WES/{}/cnv1/tumor.callsegmetrics.seg'.format(patient, tissue)
        dat1 = pd.read_csv(filein, sep='\t')
        dat1['ID'] = patient + tissue
        alldat.append(dat1)
    alldat = pd.concat(alldat, axis=0)
    alldat.to_csv('/home_2//analysisMetAF5/gistic3/sampleLevel_ta0.3/{}.seg'.format(cat), sep='\t', index=False, header=False)
    cmd = 'export LD_LIBRARY_PATH="/home//software/MCR/v83/runtime/glnxa64:/home//software/MCR/v83/bin/glnxa64:/home//software/MCR/v83/sys/os/glnxa64:${LD_LIBRARY_PATH}"'
    print (cmd)
    refgenes = '/home//software/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat'
    cmd = './gistic2  -b {} -seg {}.seg  -refgene {}  -genegistic 1   -broad 0  -twoside 1 ' \
        ' -armpeel 1 -savegene 1 -gcm mean  -js 4  -conf 0.99    -ta 0.3 -td 0.3 -savedata 0'.format(cat, cat, refgenes)
    print (cmd)

def getSigCNV(filein, pvalue = 0.001):
    """
    Parse a GISTIC confidence file and return significant genes at a q-value cutoff.
    """
    mydict = {}
    dat = pd.read_csv(filein, sep='\t', index_col=0)
    for j in dat.columns[:-1]:      ### residual q value
        if float(dat[j][1]) <= pvalue:
            for i in dat[j][3:]:
                if i is not np.nan and '[' not in i:
                    mydict[i] = j
    return mydict


def genGenes():
    """Generate a gene list from significant GISTIC amplifications and deletions."""
    mydict1 = getSigCNV('/home_2//gistic2/outcome/tumorVSlong/del_genes.conf_99.txt')
    mydict2 = getSigCNV('/home_2//gistic2/outcome/tumorVSlong/amp_genes.conf_99.txt')
    with open('/home_2//gistic2/outcome/gene.tsv', 'w') as fout:
        fout.write('SYMBOL\n')
        for i in chain(mydict1, mydict2):
            fout.write('{}\n'.format(i))



def CNVHeatMap():
    """Create per-sample symlinks and prepare inputs for CNVkit heatmap plotting."""
    os.chdir('/home_2//analysisOutcome/CNVHeatMap')
    patients = os.listdir('/home_2//outcome')
    for patient in patients:
        files = getPatientFiles(patient)
        for file in files:
            cmd = 'ln -s  /home_2//outcome/{}/WES/{}/cnv/tumor.cns  {}{}.cns'.format(patient, file, patient, file)
            subprocess.call(cmd, shell=True)

def f_CNVHeatMap(patient):
    """Plot a CNVkit heatmap PDF for one patient across all tissues."""
    os.chdir('/home_2//analysisOutcome/CNVHeatMap')
    files = getPatientFiles(patient)
    files = [patient + file + '.cns' for file in files]
    file = '  '.join(files)
    cmd = 'cnvkit.py  heatmap  {}  -d  -o  {}.pdf'.format(file, patient)
    subprocess.call(cmd, shell=True)

def ff_CNVHeatMap():
    """Batch CNV heatmaps across all patients."""
    patients = os.listdir('/home_2//outcome')
    doMultiProcess.myPool(f_CNVHeatMap, patients, 10)

"""
cnvkit.py scatter -s  tumor.cn{s,r} -c chr6  -g HLA-A
"""

doMultiProcess = RunMultiProcess()
tumorContentDict = ff_gettumorContent()
if __name__ == "__main__":
    print ('hello, world!')
    #preProcess()
    #creatNormalRef()
    #f_CNVkit()
    #f_Absolute()
    #CNVHeatMap()
    #ff_CNVHeatMap()
    merge4gistic(cat = 'longSur')
    merge4gistic(cat = 'shortSur')
    #merge4gistic(cat = 'all')

    
