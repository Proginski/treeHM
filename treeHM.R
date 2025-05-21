#' treeHM.R
#'
#' Plot a circular phylogenetic tree with heatmap annotation.
#'
#' Authors: Ambre Baumann, Paul Roginski
#' 
#' Usage (command line):
#'   Rscript treeHM.R -t tree.nwk -d data.tsv -o output.pdf
#'
#' Or source the file and use create_hm_tree(tree, hm_data) in your R script.
#'
#' - tree.nwk: Newick tree file
#' - data.tsv: TSV file, first column = leaves (matching tree tip labels), next columns = heatmap values
#' - output.pdf: Output PDF file (default: treeHM.pdf)
#'
#' Requires: ape, ggtree, ggplot2, ggnewscale, optparse

library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(optparse)

#' Create a circular tree with heatmap annotation
#'
#' @param tree      A phylo object (from ape::read.tree)
#' @param hm_data   A data.frame with rownames as leaves, columns as heatmap values
#' @param width_per_col Width of each heatmap column (default 0.04)
#' @return          A ggplot object
create_hm_tree <- function(tree, hm_data, width_per_col=0.04) {
    # Warn for missing leaves
    missing_in_tree <- setdiff(rownames(hm_data), tree$tip.label)
    missing_in_data <- setdiff(tree$tip.label, rownames(hm_data))
    if (length(missing_in_tree) > 0) {
        warning("The following leaves are in the data file but not in the tree: ", paste(missing_in_tree, collapse=", "))
    }
    if (length(missing_in_data) > 0) {
        warning("The following leaves are in the tree but not in the data file: ", paste(missing_in_data, collapse=", "))
    }
    # Keep only tips present in both tree and data
    common_tips <- intersect(tree$tip.label, rownames(hm_data))
    tree <- keep.tip(tree, common_tips)
    hm_data <- hm_data[common_tips, , drop=FALSE]

    # Warn for missing values in the data
    na_idx <- which(is.na(hm_data), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
        for (i in seq_len(nrow(na_idx))) {
            warning(sprintf(
                "Missing value for leaf '%s' in column '%s'",
                rownames(hm_data)[na_idx[i, 1]],
                colnames(hm_data)[na_idx[i, 2]]
            ))
        }
    }

    # Build base tree
    circ <- ggtree(tree, layout="circular", size=0.4)

    # Prepare palettes for discrete columns
    palette_discrete <- lapply(seq_len(ncol(hm_data)), function(i) {
        if (is.factor(hm_data[[i]]) || is.character(hm_data[[i]])) {
            vals <- unique(hm_data[[i]])
            setNames(rainbow(length(vals)), vals)
        } else {
            NULL
        }
    })

    # Build heatmap layers (one per column)
    layers <- lapply(seq_len(ncol(hm_data)), function(i) hm_data[, i, drop=FALSE])

    # Plot
    nb_layers <- length(layers)
    palette_continous <- c("Purples", "Greens", "Oranges", "BuPu", "GnBu")
    nb_discrete = 0
    nb_continuous = 0
    circ_with_heatmap <- circ
    for (i in 1:nb_layers) {
        if (i > 1) circ_with_heatmap <- circ_with_heatmap + new_scale_fill()
        if (is.numeric(hm_data[[i]])) {
            nb_continuous = nb_continuous + 1
            circ_with_heatmap <- gheatmap(
                circ_with_heatmap, layers[[i]], colnames=FALSE, color=NA,
                width=width_per_col * ncol(layers[[i]]),
                offset = ifelse(i == 1, 0, width_per_col * sum(sapply(layers[1:(i-1)], ncol)) * 2)
            ) +
            scale_fill_distiller(palette = palette_continous[nb_continuous], direction=1) +
            labs(fill = colnames(layers[[i]])[1])
        } else {
            nb_discrete = nb_discrete + 1
            circ_with_heatmap <- gheatmap(
                circ_with_heatmap, layers[[i]], colnames=FALSE, color=NA,
                width=width_per_col * ncol(layers[[i]]),
                offset = ifelse(i == 1, 0, width_per_col * sum(sapply(layers[1:(i-1)], ncol)) * 2)
            ) +
            scale_fill_manual(values = palette_discrete[[i]]) +
            labs(fill = colnames(layers[[i]])[1])
        }
    }
    circ_with_heatmap + theme(
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        legend.key.height = unit(0.8, "cm")
    )
}

# If run as a script, parse command line and plot
if (sys.nframe() == 0) {
    option_list <- list(
      make_option(c("-t", "--tree"), type="character", help="Tree file in Newick format"),
      make_option(c("-d", "--data"), type="character", help="TSV file: first column = leaves, next columns = heatmap values"),
      make_option(c("-o", "--output"), type="character", default="treeHM.pdf", help="Output PDF file [default: %default]")
    )
    opt <- parse_args(OptionParser(option_list=option_list))

    # Read tree and data
    tree <- read.tree(opt$tree)
    hm_data <- read.table(opt$data, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

    # Plot and save to A4 PDF
    pdf(opt$output, paper="A4")
    print(create_hm_tree(tree, hm_data))
    dev.off()
}