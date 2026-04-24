"""
Prepare REVOLVER cohort input from per-patient clonevol outputs.

REVOLVER is an approach for detecting recurrent evolutionary trajectories across a cohort.
This script converts a per-patient `clonevol.txt` table into the REVOLVER "cohort" input
format (one line per variant with cluster assignments and cancer cell fractions).

Inputs
- `/home_2//analysis/pyclone_outcome_majorCopy/<patient>/clonevol.txt`
  Expected columns include: `cluster`, `gene`, `Driver`, and one column per tissue/timepoint.

Outputs
- `/home_2//analysis/revolver_outcome/input.tsv`
  Columns: Misc, patientID, variantID, cluster, is.driver, is.clonal, CCF

Notes
- Paths are hard-coded and assume an HPC file system layout.
- CCF values are normalized to [0, 1] and zeroed below a 0.1 cutoff.
"""

from util import *
import os, re, sys, glob, subprocess
import numpy as np, pandas as pd
from multiprocessing import Pool
from itertools import chain
from collections import defaultdict
def preFile(patient):
    """
    Build a REVOLVER-compatible table for a single patient.

    Parameters
    ----------
    patient : str
        Patient identifier (directory name under `/home_2//outcome`).

    Returns
    -------
    pandas.DataFrame | None
        A per-variant table for this patient. Returns None if no variants are present.
    """
    filein = '/home_2//analysis/pyclone_outcome_majorCopy/{}/clonevol.txt'.format(patient)
    dat = pd.read_csv(filein, sep='\t')
    #dat = dat[dat['Driver']]
    if dat.shape[0] == 0: return
    # Tissue/timepoint columns used to compute per-variant CCF strings.
    files = getPatientFiles('outcome', patient)
    if len(files) == 1: pass
    dat.sort_values(by='cluster', ascending=True, inplace=True)
    dat.drop_duplicates(subset='gene', inplace=True)
    Misc = []; patientID = []; variantID = []; cluster = []; isdriver = []; isclonal = []; CCF = []
    for i in dat.iterrows():
        lines = i[1]
        Misc.append('Nothing'); patientID.append(patient); variantID.append(lines['gene'])
        cluster.append(lines['cluster']); isdriver.append(str(lines['Driver']).upper())
        # Convention: cluster 1 is treated as clonal (founding) for REVOLVER.
        if lines['cluster'] == 1: isclonal.append('TRUE')
        else: isclonal.append('FALSE')

        # Convert clonevol percentage (0-100) to a fraction (0-1) and apply a 0.1 cutoff.
        CCFs = [round(lines[file] / 100, 3) if lines[file] / 100 >= 0.1 else 0 for file in files]
        CCFs = ['{}:{}'.format(file, tmp) for file, tmp in zip(files, CCFs)]
        CCF.append(';'.join(CCFs))
    dat = pd.DataFrame({'Misc': Misc, 'patientID': patientID, 'variantID': variantID, 
    'cluster': cluster, 'is.driver': isdriver, 'is.clonal': isclonal, 'CCF': CCF})
    return dat


def f_preFile():
    """
    Batch conversion for all patients and write REVOLVER cohort input to disk.
    """
    os.chdir('/home_2//analysis/revolver_outcome')
    patients = os.listdir('/home_2//outcome')
    all_dat = [preFile(patient) for patient in patients]
    all_dat = pd.concat(all_dat, axis=0)
    all_dat.to_csv('input.tsv', sep='\t', index=False)


def analysisResult():
    """
    Quick sanity-check utility:
    Merge REVOLVER output categories with patient outcomes and print contingency counts.
    """
    filein = '/home_2//analysis/revolver_outcome/output.tsv'
    dat = pd.read_csv(filein, sep='\t', header=None, dtype=np.object)
    dat.columns = ['patient', 'cat']
    dat1 = getDataframe(); dat1 = dat1.astype(np.object)
    dat1.drop(labels=['tissue'], inplace=True, axis=1)
    dat1.drop_duplicates(subset=['patient'], inplace=True)
    alldat = pd.merge(left=dat, right=dat1)
    print (alldat.groupby('cat')['outcome'].value_counts())


if __name__ == '__main__':
    f_preFile()
