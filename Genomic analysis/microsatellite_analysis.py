"""
Microsatellite instability (MSI) calling wrapper.

This script is a thin orchestrator around `msisensor-pro` to generate per-sample MSI
metrics (e.g., `msi.tsv`) for all WES samples under an `outcome/<patient>/WES/<tissue>`
directory structure.

High-level workflow
1) Enumerate patients under `/home_2//outcome`
2) For each patient, enumerate tissue/sample folders via `getPatientFiles(...)`
3) For each tissue:
   - Create `/home_2//outcome/<patient>/WES/<tissue>/msi/`
   - Run `msisensor-pro msi` with a fixed reference site list

Notes / assumptions
- Paths are hard-coded for an HPC/NFS environment.
- Normal BAM is assumed at `../../NC/gatk/deduped.bam`, tumor BAM at `../gatk/deduped.bam`.
- Execution uses `f_subprocess` from CNVkit utilities which logs to `log.txt`.
"""

from CNVkit import f_subprocess
from util import *


doMultiProcess = RunMultiProcess()
def doMsi(patient):
    """
    Run msisensor-pro for all tissues belonging to a patient.

    Parameters
    ----------
    patient : str
        Patient identifier (directory name under `/home_2//outcome/`).
    """
    # Tissue/sample identifiers for this patient (excluding "NC" normal folder).
    files = getPatientFiles('outcome', patient)
    for file in files:
        # Working directory where `msisensor-pro` output is written.
        dirName = '/home_2//outcome/{}/WES/{}/msi'.format(patient, file)
        if not os.path.isdir(dirName): os.makedirs(dirName)
        os.chdir(dirName)

        # msisensor-pro command:
        # -d: microsatellite reference site file
        # -n: normal BAM
        # -t: tumor BAM
        # -o: output TSV
        # -c/-b: tool-specific filtering thresholds
        cmd = '/usr/software/anaconda3/bin/msisensor-pro   msi -d /home//software/msisensor/reference.site ' \
            '-n ../../NC/gatk/deduped.bam   -t  ../gatk/deduped.bam  -o msi.tsv  -c  50 -b 5 '
        f_subprocess([cmd])



def f_doMsi():
    """
    Entry point: run MSI for all patients under `/home_2//outcome`.
    """
    patients = os.listdir('/home_2//outcome')
    doMultiProcess.myPool(doMsi, patients, processes=8)


if __name__ == '__main__':
    f_doMsi()
