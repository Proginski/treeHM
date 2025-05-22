# treeHM.R

This script plots a circular phylogenetic tree with heatmap annotation.

## Usage

**Command line:**
```sh
Rscript treeHM.R --tree data/test/test_tree.nwk --data data/test/test_heatmap_values.tsv --output data/test/test_output.pdf
```
- `"-t", "--tree"`: Newick tree file
- `"-d", "--data"`: TSV file, first column = leaves (matching tree tip labels), next columns = heatmap values
- `"-o", "--output"`: Output PDF file (default: `treeHM.pdf`)

**As a function in R:**
```r
source("treeHM.R")
tree <- ape::read.tree("tree.nwk")
hm_data <- read.table("data.tsv", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
pdf("output.pdf", paper="A4")
print(create_hm_tree(tree, hm_data))
dev.off()
```

## Example Output

![Example Output](data/images/example1.png)

## Requirements

- R packages: `ape`, `ggtree`, `ggplot2`, `ggnewscale`, `optparse`

## Features

- Warns about missing leaves and missing values.
- Supports both discrete and continuous heatmap columns.
- Output is always A4-sized PDF.
