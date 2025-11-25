#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

# Parse options
    option_list <- list(
                        make_option(c('-i', '--inFiles', type='character'),
                                        default = NULL, help = 'Comma-separated list of HTseq-count summarized .tsv files',
                                        metavar = 'character'),
                        make_option(c('-s', '--sampleNames', type='character'),
                                        default = NULL,
                                        help = 'Commaseparated list of sample names',
                                        metavar = 'character'))

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

# Checks
    if (is.null(opt$inFiles))
        stop('You have to provide HTseq-count files')

    if (is.null(opt$sampleNames))
        stop('Unspecified sample names')


# Parse params
    files <- opt$inFiles |> strsplit(split = ',') |> unlist()
    names <- opt$sampleNames |> strsplit(split = ',') |> unlist()

# Load data
    data <- files |>
            read_tsv(id = 'path', col_names = c('gene_name', 'gene_id', 'gene_biotype', 'read_count')) |>
            mutate(sample = str_split_fixed(basename(path), pattern = '_htseq', n = 2)[ ,1]) |>
            dplyr::rename(sample_id = sample) |>
            select(-path) |>
            filter(!str_starts(gene_name, '__')) |>
            filter(gene_biotype == 'protein_coding') |>
            arrange(sample_id) |>
            mutate(sample_id = factor(sample_id))


# Pivot to wide
    data_wide <- data |>
        pivot_wider(names_from = 'sample_id', values_from = 'read_count') |>
        column_to_rownames(var = 'gene_name') |>
        select(-gene_id, -gene_biotype)

# Prepare data for DESeq2
    countdata <- data_wide |>
                    as.matrix()

    coldata <- data.frame(sample_id = names)

# Run normalization with DESeq2
    library(DESeq2)

    dds <- DESeqDataSetFromMatrix(
            countData = countdata,
            colData = coldata,
            design = ~ 1)

    dds <- estimateSizeFactors(dds)

# Normalize factors to mean
    factors <- dds |>
        colData() |>
        data.frame() |>
        mutate(factor = mean(sizeFactor) / sizeFactor) |>
        dplyr::select(factor) |>
        rownames_to_column(var = 'sampleId')

# Write out to file
    write_tsv(factors, 'norm_factors.tsv', col_names = FALSE)

