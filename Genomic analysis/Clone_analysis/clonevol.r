library('clonevol')

## Clonevol visualization of tumor evolution.
##
## This script reads a `clonevol.txt` table (prepared from PyClone outputs) and generates:
## - cluster flow plots
## - inferred clonal models (consensus trees) using bootstrap testing
## - annotated clonal evolution plots highlighting driver mutations (if available)
##
## Inputs
## - `clonevol.txt` : per-variant table with cluster assignments and CCF columns per sample
##
## Notes
## - The working directory is set to one patient folder; change it before running.



setwd("/home_2/wzt/analysisOutcome/pyclone_majorCopy/131")
dat = read.table('clonevol.txt', sep = '\t', header = 1, row.names = 1, stringsAsFactors = FALSE, quote = '', check.names = TRUE)
dat <- dat[order(dat$cluster),]
sampleNames = colnames(dat)[4:length(colnames(dat))]
plot.cluster.flow(dat, vaf.col.names = sampleNames, min.cluster.vaf = 0.05, low.vaf.no.line = FALSE)

y = infer.clonal.models(variants = dat,
                        cluster.col.name = 'cluster',
                        ccf.col.names = sampleNames, 
                        #cancer.initiation.model='polyclonal',
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        min.cluster.vaf = 0.05,
                        sum.p = 0.01,
                        alpha = 0.01)

if (TRUE){
  if (nrow(dat[dat$Driver == 'True', ]) >= 1){
  y <- transfer.events.to.consensus.trees(y, dat[dat$Driver == 'True', ], cluster.col.name = 'cluster', event.col.name = 'gene')}
y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
plot.clonal.models(y,
                  # box plot parameters
                  box.plot = TRUE,
                  fancy.boxplot = TRUE,
                  fancy.variant.boxplot.highlight = 'is.driver',
                  fancy.variant.boxplot.highlight.shape = 21,
                  fancy.variant.boxplot.highlight.fill.color = 'red',
                  fancy.variant.boxplot.highlight.color = 'black',
                  fancy.variant.boxplot.highlight.note.col.name = 'gene',
                  fancy.variant.boxplot.highlight.note.color = 'blue',
                  fancy.variant.boxplot.highlight.note.size = 2,
                  fancy.variant.boxplot.jitter.alpha = 1,
                  fancy.variant.boxplot.jitter.center.color = 'grey50',
                  fancy.variant.boxplot.base_size = 12,
                  fancy.variant.boxplot.plot.margin = 1,
                  fancy.variant.boxplot.vaf.suffix = '.VAF',
                  # bell plot parameters
                  clone.shape = 'bell',
                  bell.event = TRUE,
                  bell.event.label.color = 'blue',
                  bell.event.label.angle = 60,
                  clone.time.step.scale = 1,
                  bell.curve.step = 2,
                  # node-based consensus tree parameters
                  merged.tree.plot = TRUE,
                  tree.node.label.split.character = NULL,
                  tree.node.shape = 'circle',
                  tree.node.size = 30,
                  tree.node.text.size = 0.5,
                  merged.tree.node.size.scale = 1.25,
                  merged.tree.node.text.size.scale = 2.5,
                  merged.tree.cell.frac.ci = FALSE,
                  # branch-based consensus tree parameters
                  merged.tree.clone.as.branch = TRUE,
                  mtcab.event.sep.char = ',',
                  mtcab.branch.text.size = 1,
                  mtcab.branch.width = 0.75,
                  mtcab.node.size = 3,
                  mtcab.node.label.size = 1,
                  mtcab.node.text.size = 1.5,
                  # cellular population parameters
                  cell.plot = TRUE,
                  num.cells = 100,
                  cell.border.size = 0.25,
                  cell.border.color = 'black',
                  clone.grouping = 'horizontal',
                  #meta-parameters
                  scale.monoclonal.cell.frac = TRUE,
                  show.score = FALSE,
                  cell.frac.ci = TRUE,
                  disable.cell.frac = FALSE,
                  # output figure parameters
                  out.dir = '.',
                  out.format = 'pdf',
                  overwrite.output = TRUE,
                  width = 12,
                  height = 8,
                  # vector of width scales for each panel from left to right
                  panel.widths = c(3,4,2,4,4)
)
}

















                      
                      
