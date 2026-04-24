library(revolver)
require(tidyverse)
## REVOLVER cohort analysis.
##
## This script loads the cohort input table (`input.tsv`) prepared by `REVOLVER.py`,
## runs REVOLVER model fitting and clustering, and produces summary plots.
##
## Inputs (working directory)
## - `input.tsv` : cohort table with variant clusters and CCF strings
##
## Outputs
## - `output.tsv` : cluster labels (written by write.table below)
## - Multiple PDF plots (clusters, trajectories, dendrograms, etc.)
setwd("/home_2//analysisOutcome/revolver/")
dat = read.table('input.tsv', sep = '\t', header = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
dat$cluster = as.character(dat$cluster)
dat$patientID = as.character(dat$patientID)
inputdat = revolver_cohort(dat, MIN.CLUSTER.SIZE = 0,  annotation = "CHOL")
revolver_check_cohort(inputdat)
non_recurrent = Stats_drivers(inputdat) %>% filter(N_tot == 1) %>% pull(variantID)
inputdat = remove_drivers(inputdat, non_recurrent)
inputdat = compute_clone_trees(inputdat)
inputdat = revolver_fit(inputdat, parallel = TRUE, n = 10, max.iterations = 10, initial.solution = NA)
inputdat = revolver_cluster(inputdat, split.method = 'cutreeHybrid', min.group.size = 3)
plot_clusters(inputdat, cutoff_drivers = 0, cutoff_trajectories = 1)
plot_trajectories_per_cluster(inputdat, min_counts = 1)
plot_drivers_graph(inputdat)
plot_dendrogram(inputdat)
plot_drivers_clonality(inputdat)
plot_patient_trees(inputdat, '97')
#Phylo(inputdat, '97', rank = NULL)
write.table(inputdat$cluster$fits$labels, file='output.tsv', sep='\t', quote = FALSE, col.names = NA, row.names = TRUE)

pdf('plot_clusters.pdf', onefile = FALSE, width = 12,height = 10)















