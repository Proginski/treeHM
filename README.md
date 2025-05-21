# treeHM.R

This script plots a circular phylogenetic tree with heatmap annotation.

## Usage

**Command line:**
```sh
Rscript treeHM.R -t tree.nwk -d data.tsv -o output.pdf
```
- `tree.nwk`: Newick tree file
- `data.tsv`: TSV file, first column = leaves (matching tree tip labels), next columns = heatmap values
- `output.pdf`: Output PDF file (default: `treeHM.pdf`)

**As a function in R:**
```r
source("treeHM.R")
tree <- ape::read.tree("tree.nwk")
hm_data <- read.table("data.tsv", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
pdf("output.pdf", paper="A4")
print(create_hm_tree(tree, hm_data))
dev.off()
```

## Requirements

- R packages: `ape`, `ggtree`, `ggplot2`, `ggnewscale`, `optparse`

## Features

- Warns about missing leaves and missing values.
- Supports both discrete and continuous heatmap columns.
- Output is always A4-sized PDF.