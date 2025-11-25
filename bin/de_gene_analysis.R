#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(cowplot)
library(ggtext)
library(glue)
library(GGally)


# Parse options
    option_list <- list(
                        make_option(c('-i', '--inFiles', type='character'),
                                        default = NULL, help = 'Comma-separated list of HTseq-count summarized .tsv files',
                                        metavar = 'character'),
                        make_option(c('-o', '--outDir', type='character'),
                                        default = NULL,
                                        help = 'Output directory for figures',
                                        metavar = 'character'),
                        make_option(c('-a', '--annotation', type='character'),
                                        default = NULL,
                                        help = 'Absolute path to .csv file with information about samples. First three columns must be: SampleId, read1, read2; followed by batches and conditions.',
                                        metavar = 'character'),
                        make_option(c('-c', '--conditions', type='character'),
                                        default = NULL,
                                        help = 'Comma-separated list of experimental factors used for DE analysis',
                                        metavar = 'character'),
                        make_option(c('-b', '--batches', type='character'),
                                        default = NULL,
                                        help = 'Comma-separated list of factors which effects should be removed (optional)',
                                        metavar = 'character'),
                        make_option(c('-r', '--reflevels', type='character'),
                                        default = NULL,
                                        help = 'Comma-separated list of factors values that should be treated as reference levels. More than one value per factor is allowed. If not provided, reference levels will be determined alphabetically. (optional)',
                                        metavar = 'character'),
                        make_option(c('-l', '--log2fc', type='double'),
                                        default = 0.5,
                                        help = 'Threshold for log2 fold-change of gene expression (optional)',
                                        metavar = 'number'),
                        make_option(c('-p', '--pval', type='double'),
                                        default = 0.05,
                                        help = 'Threshold for p-value (optional)',
                                        metavar = 'number'))

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

# Checks
    if (is.null(opt$inFiles))
        stop('You have to provide HTseq-count files')

    if (is.null(opt$outDir))
        stop('Unspecified output folder')

    if (is.null(opt$annotation))
        stop('You have to provide samples annotation in .csv format')

    if (is.null(opt$conditions))
        stop('You have to provide name of at least one experimental condition')


# Set working dir
    setwd(opt$outDir)

# Set params
    lfcT <- opt$log2fc          # log2 fold-change of gene expression cut-off
    alpha <- opt$pval           # p-adj value cut-off
    ref_levels <- opt$reflevels |> strsplit(split = ',') |> unlist()    # reference conditions for DE comparisons

# Parse factors
    conditions <- opt$conditions |> strsplit(split = ',') |> unlist()
    
    if (!is.null(opt$batches)) {
        batches <- opt$batches |> strsplit(split = ',') |> unlist() 
    } else {
        batches <- NULL
    }
    

# Load data
    files <- opt$inFiles |> strsplit(split = ',') |> unlist()

    data <- files |>
            read_tsv(id = 'path', col_names = c('gene_name', 'gene_id', 'gene_biotype', 'read_count')) |>
            mutate(sample = str_split_fixed(basename(path), pattern = '_htseq', n = 2)[ ,1]) |>
            dplyr::rename(sample_id = sample) |>
            select(-path) |>
            filter(!str_starts(gene_name, '__'))

            if (is.numeric(data$sample_id)) {
                data <- data |>
                    arrange(as.numeric(sample_id)) |>
                    mutate(sample_id = fct_inseq(sample_id))
            } else {
                data <- data |>
                    arrange(sample_id) |>
                    mutate(sample_id = factor(sample_id))


            }


## Pivot to wide
    data_wide <- data |>
        pivot_wider(names_from = 'sample_id', values_from = 'read_count') |>
        column_to_rownames(var = 'gene_name') |>
        select(-gene_id, -gene_biotype)


# Labels
    sample_info <- read_csv(opt$annotation, col_types = 'fccfffff') |>
        select(-read1, -read2)


    annotation <- data_wide |>
        colnames() |>
        data.frame(sample_id = _) |>
        left_join(sample_info, by = join_by( sample_id == sampleId )) |> # this ensures that annotation is in the same order as data
        unite('group', all_of(conditions), sep = '-', remove = FALSE) |> # combine conditions for easier modeling with DEseq2
        mutate(group = factor(group))
    
# Define model for DE analysis
    formula_de <- paste('~', paste(c(batches, 'group'), collapse = ' + ')) |> as.formula()


# Build pairs of contrasts for DE comparison
    # For comparisons we don't want permutations but combinations of groups to avoid comparing the same groups twice. e.g. 'T1-Y vs. T1-N' and 'T1-N vs. T1-Y' is the same.
    # In fact we don't want all combinations, only combinations within specific conditions. e.g. We want 'T1-Y vs T1-N' or 'T0-Y vs T0-N'.
    #       ,but we don't want e.g. 'T1-Y vs. T0-N'  or  'T1-N vs. T1-Y'
    if (length(conditions) == 1) {
        comparisons <-
            annotation |>
            mutate(group = as.character(group)) |>
            pull(group) |>
            unique() |>
            combn(m = 2) |>
            t() |>
            data.frame() |>
            dplyr::rename(c('comp1' = 1, 'comp2' = 2))
    } else {
        comparisons <- list()
        for (condition in conditions) {
            comparisons[[condition]] <-
                annotation |>
                select(any_of(conditions), group) |>
                mutate(group = as.character(group)) |>
                unique() |>
                group_by(.data[[condition]]) |>
                group_map(~ data.frame(t(combn(.$group, m = 2)))) |>
                bind_rows()
        }
        comparisons <- comparisons |> 
                        bind_rows() |>
                        dplyr::rename(c('comp1' = 1, 'comp2' = 2))
    }

# Change order of groups within comparisons to make sure that reference levels are always second in order
    if (!is.null(ref_levels)) {

        for (i in 1:nrow(comparisons)) {
            for (ref in ref_levels) {
                change <- str_detect(comparisons[i, ], ref) == c(TRUE, FALSE) # Determine if reference level is in the first position instead of second
                
                if (sum(change) == 2) { # If yes, switch the two values
                    comp_temp <- comparisons$comp1[i]
                    comparisons$comp1[i] <- comparisons$comp2[i]
                    comparisons$comp2[i] <- comp_temp
                    break # Stop at first identified reference level to prevent further switching in case multiple levels were provided.
                
                }
            }
        }
    }

# Helper functions
## Create better scientific log10 notation
    scientific_10 = function(x) {
        ifelse(
            x==0, '0',
            parse(text = sub('e[+]?', ' %*% 10^', scales::scientific_format()(x)))
        )
    }

## PCA plot
    # works only for 2 interest groups
    plot_pca <- function(vsd, intgroup, subtit) {

        pcaData <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE) 

        percentVar <- round(100 * attr(pcaData, 'percentVar'))
        
        if (length(intgroup) == 1) {
            p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[intgroup]]))
        } else {
            p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[intgroup[1]]], shape = .data[[intgroup[2]]]))
        }

        p <- p +
            geom_point(size = 3) +
            geom_text_repel(aes(label = name), color = 'black', size = 2.5) +
            xlab(paste0('PC1: ', percentVar[1], '% variance')) +
            ylab(paste0('PC2: ', percentVar[2], '% variance')) +
            coord_fixed() +
            theme_half_open(font_size = 10) +
            background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
            labs(title = 'PCA (vst)',
                subtitle = subtit)
    }

## MA plot
    ma_plot <- function(deseq_lfc, lfcT, alpha, subtit) {
        ma_data <- deseq_lfc |> 
            data.frame() |> 
            mutate(sig = case_when(
                            padj < alpha ~ 'sig',
                            .default = 'nosig')) |>
            arrange(sig)
        
        all <- ma_data |> filter(baseMean != 0) |> nrow()
        up <- ma_data |> filter(sig == 'sig', log2FoldChange > 0) |> nrow()
        down <- ma_data |> filter(sig == 'sig', log2FoldChange < 0) |> nrow()
        
        color_leg <- paste0('total:  ', all, '\n', 
                            'up:     ', up, ' (', round(up/all * 100, digits = 1), ' %)\n',
                            'down: ', down, ' (', round(down/all * 100, digits = 1), ' %)\n',
                            '\n\n\nsignificance')
        
        # Y axis scale to 99.9% of the distribution
        log2FC <- ma_data |> pull(log2FoldChange)
        lower_bound <- quantile(log2FC, 0.001, na.rm = TRUE)  # Lower bound for outliers
        upper_bound <- quantile(log2FC, 0.999, na.rm = TRUE)  # Upper bound for outliers
        
        y_limit <- max(1.5, abs(lower_bound), abs(upper_bound))
        
        outliers <- ma_data |>
            filter(abs(log2FoldChange) > y_limit) |>
            mutate(log2FoldChange = sign(log2FoldChange) * y_limit * 1.05)
            
        
        ma_data |>
            filter(abs(log2FoldChange) <= y_limit) |>
            ggplot() + 
                geom_hline(yintercept = 0, color = 'gray40', linewidth = 1) +
                geom_hline(yintercept = lfcT, color = 'pink', linewidth = 1) +
                geom_hline(yintercept = -lfcT, color = 'pink', linewidth = 1) +
                geom_point(aes(y = log2FoldChange, x = baseMean, color = sig), show.legend = TRUE) +
                geom_point(data = outliers |> filter(log2FoldChange > 0),
                           aes(y = log2FoldChange, x = baseMean, color = sig), shape = 2, show.legend = FALSE) +
                geom_point(data = outliers  |> filter(log2FoldChange < 0),
                           aes(y = log2FoldChange, x = baseMean, color = sig), shape = 6, show.legend = FALSE) +
                scale_x_log10(labels = scales::label_log()) +
                coord_cartesian(ylim = c(-y_limit, y_limit), expand = TRUE) +
                scale_color_manual(limits = c('nosig', 'sig'),
                                   values = c('gray50', '#00BFC4')) +
                labs(title = 'MA plot',
                     subtitle = subtit,
                     x = 'mean of normalized counts (log<sub>10</sub>)',
                     y = 'log<sub>2</sub> (fold change)',
                     caption = glue('FDR < {alpha}'),
                     color = color_leg) +
                theme_half_open(font_size = 10) +
                background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
                theme(plot.subtitle = element_markdown(),
                      axis.title.x = element_markdown(),
                      axis.title.y = element_markdown())
                
                
    }

## Volcano plot
    volcano_plot <- function(deseq_lfc, lfcT, alpha, subtit) {
        volc_data <- deseq_lfc |> 
            data.frame() |> 
            rownames_to_column(var = 'gene_name') |>
            mutate(padj = replace_na(padj, 1)) |>
            mutate(sig = case_when(
                            padj < alpha & abs(log2FoldChange) > lfcT ~ 'Pval + logFC',
                            padj >= alpha & abs(log2FoldChange) > lfcT ~ 'logFC',
                            padj < alpha & abs(log2FoldChange) <= lfcT ~ 'Pval',
                            padj >= alpha & abs(log2FoldChange) <= lfcT ~ 'NS',
                            .default = 'NS') )

        all <- volc_data |> filter(baseMean != 0) |> nrow()
        up <- volc_data |> filter(sig == 'Pval + logFC', log2FoldChange > 0) |> nrow()
        down <- volc_data |> filter(sig == 'Pval + logFC', log2FoldChange < 0) |> nrow()
        color_leg <- paste0('total:  ', all, '\n', 
                            'up:     ', up, ' (', round(up/all * 100, digits = 1), ' %)\n',
                            'down: ', down, ' (', round(down/all * 100, digits = 1), ' %)\n',
                            '\n\n\nsignificance')
        
        # Y axis scale 99.9% of the distribution
        log2FC <- volc_data |> pull(log2FoldChange)
        lower_bound <- quantile(log2FC, 0.001, na.rm = TRUE)  # Lower bound for outliers
        upper_bound <- quantile(log2FC, 0.999, na.rm = TRUE)  # Upper bound for outliers
        
        x_limit <- max(1.5, abs(lower_bound), abs(upper_bound))
        
        outliers <- volc_data |>
            filter(abs(log2FoldChange) > x_limit) |>
            mutate(log2FoldChange = sign(log2FoldChange) * x_limit * 1.05)
        

        volc_data |>
            filter(abs(log2FoldChange) <= x_limit) |>
            ggplot() + 
                geom_hline(yintercept = -log10(alpha), color = 'gray60', linewidth = 0.5, linetype = 'longdash') +
                geom_vline(xintercept = lfcT, color = 'gray60', linewidth = 0.5, linetype = 'longdash') + 
                geom_vline(xintercept = -lfcT, color = 'gray60', linewidth = 0.5, linetype = 'longdash') +
                geom_point(aes(y = -log10(padj), x = log2FoldChange, color = sig), show.legend = TRUE) +
                geom_point(data = outliers |> filter(log2FoldChange > 0),
                           aes(y = -log10(padj), x = log2FoldChange, color = sig), shape = 2, show.legend = FALSE) +
                geom_point(data = outliers |> filter(log2FoldChange < 0),
                           aes(y = -log10(padj), x = log2FoldChange,  color = sig), shape = 2, show.legend = FALSE) +
                geom_text_repel(aes(y = -log10(padj), x = log2FoldChange, label = gene_name), 
                                data = ~ filter(.x, sig == 'Pval + logFC'), color = 'black', size = 2.5) +
                geom_text_repel(aes(y = -log10(padj), x = log2FoldChange, label = gene_name), 
                                data = outliers |> filter(sig == 'Pval + logFC'), color = 'black', size = 2.5) +
                scale_color_manual(limits = c('Pval + logFC', 'logFC', 'Pval', 'NS'),
                                   values = c('#F8766D', '#00BF7D', '#00B6EB', 'gray50')) +
                coord_cartesian(xlim = c(-x_limit, x_limit), expand = TRUE) +
                labs(title = 'Volcano plot',
                     subtitle = subtit,
                     x = 'log<sub>2</sub> (fold change)',
                     y = '-log<sub>10</sub> (FDR)',
                     color = color_leg,
                     caption = glue('|log<sub>2</sub>FC| > {lfcT} ; FDR < {alpha}')) +
                theme_half_open(font_size = 10) +
                background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
                theme(plot.subtitle = element_markdown(),
                      plot.caption = element_markdown(),
                      axis.title.x = element_markdown(),
                      axis.title.y = element_markdown())
    }

## Gene expression plot
    gene_expression_plot <- function(deseq_lfc, deseq_dds, lfcT, alpha, condition, genesToPlot = 12, subtit) {
        top_genes <- deseq_lfc |>
            data.frame() |>
            rownames_to_column(var = 'gene_name') |>
            mutate(sig = case_when(
                            padj < alpha & abs(log2FoldChange) > lfcT ~ 'Pval + logFC',
                            padj >= alpha & abs(log2FoldChange) > lfcT ~ 'logFC',
                            padj < alpha & abs(log2FoldChange) <= lfcT ~ 'Pval',
                            padj >= alpha & abs(log2FoldChange) <= lfcT ~ 'NS',
                            .default = 'NS') ) |>
            filter(sig == 'Pval + logFC') |>
            arrange(desc(abs(log2FoldChange))) |>
            pull(gene_name) |>
            head(n = genesToPlot)

        if (length(top_genes) > 0 ) {
            p <- NULL

            for (i in seq_along(top_genes)) {
                plot_data <- plotCounts(deseq_dds, gene = top_genes[i], intgroup = condition,
                                        returnData = TRUE)
                
                if (length(condition) == 1) {
                    p[[i]] <- plot_data |>
                            ggplot() +
                            geom_beeswarm(aes(x = .data[[condition]], y = count, color = .data[[condition]]), cex = 4) +
                            geom_point(data = plot_data |>
                                                group_by(.data[[condition]]) |>
                                                summarize(mean_count = mean(count), .groups = 'drop'),
                                       aes(x = .data[[condition]], y = mean_count),
                                       shape = 95, size = 10)
                } else {
                    p[[i]] <- plot_data |>
                            ggplot() +
                            geom_beeswarm(aes(x = .data[[condition[1]]], y = count, color = .data[[condition[1]]]), cex = 4) +
                            geom_point(data = plot_data |>
                                                group_by(.data[[condition[1]]], .data[[condition[2]]]) |>
                                                summarize(mean_count = mean(count), .groups = 'drop'),
                                       aes(x = .data[[condition[1]]], y = mean_count),
                                       shape = 95, size = 10) +
                            facet_wrap(vars(.data[[condition[2]]]))
                }
                p[[i]] <- p[[i]] +
                            scale_y_log10() +
                            labs(title = top_genes[i],
                                 y = 'normalized counts') +
                            theme_half_open(font_size = 10) +
                            background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
                            theme(legend.position = 'none',
                                  axis.title.y = element_markdown())
            }


            wrap_plots(p) +
                plot_layout(axis_titles = 'collect') +
                plot_annotation(title = 'Differentially Expressed genes',
                                subtitle = glue('Top {genesToPlot} genes; {subtit}'),
                                caption = glue('|log<sub>2</sub>FC| > {lfcT} ; FDR < {alpha}'),
                                theme = theme(plot.caption = element_markdown()))
        } else {
            ggplot() + geom_blank()
        }
    }


# GO term enrichment analysis plot
    GO_dot_plot <- function(deseq_lfc, updown, lfcT, alpha, subtit, fname) {

        if (updown == 'up') {
                significant_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(log2FoldChange > lfcT, padj < alpha) |>
                                rownames_to_column(var = 'gene_symbol')
        } else if (updown == 'down') {
                significant_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(log2FoldChange < -lfcT, padj < alpha) |>
                                rownames_to_column(var = 'gene_symbol')
        }

        background_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(baseMean != 0) |>
                                rownames_to_column(var = 'gene_symbol')


        significant_genes_map <- bitr(geneID = significant_genes |> pull(gene_symbol),
                                    fromType = 'SYMBOL',
                                    toType = 'ENTREZID',
                                    OrgDb = 'org.Hs.eg.db')

        background_genes_map <- bitr(geneID = background_genes |> pull(gene_symbol),
                                    fromType = 'SYMBOL',
                                    toType = 'ENTREZID',
                                    OrgDb = 'org.Hs.eg.db')

        # GO term enrichment
        # MF: Molecular Function : molecular activities of gene products
        # CC: Cellular Component : where gene products are active
        # BP: Biological Process : pathways and larger processes made up of the activities of multiple gene products

        go_types <- tribble(
            ~type, ~tit, ~subtit,
            'BP',  'Biological Process',    'pathways of activities of multiple gene products',
            'MF',  'Molecular Function',    'molecular activities of gene products',
            'CC',  'Cellular Component',    'where gene products are active'
        )

        p <- NULL
        for (i in 1:nrow(go_types)) {
            type     <- go_types[i, ]$type
            t_tit    <- go_types[i, ]$tit
            t_subtit <- go_types[i, ]$subtit

            go_an <- enrichGO(gene          = significant_genes_map$ENTREZID,
                              universe      = background_genes_map$ENTREZID,
                              OrgDb         = org.Hs.eg.db,
                              ont           = type,
                              pAdjustMethod = 'BH',
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

            if (!is.null(go_an)) {
                go_an@result |> write_csv(glue('{fname}.GO_{type}_results.{updown}.csv'))

                if (go_an@result |> filter(p.adjust < 0.01) |> nrow() > 0) {
                    p[[type]] <- dotplot(go_an, showCategory = 20, font.size = 10, title = t_tit) + labs(subtitle = t_subtit)
                }
            }
        }

        if (!is.null(p)) {
            if (updown == 'up') {
                ud_title <- 'upregulated' 
                ud_caption <- glue('log<sub>2</sub>FC > {lfcT}')
            } else if (updown == 'down') {
                ud_title <- 'downregulated'
                ud_caption <- glue('log<sub>2</sub>FC < -{lfcT}')
            }

            wrap_plots(p) +
                plot_annotation(title = glue('GO term analysis: {ud_title} genes'),
                                subtitle = subtit,
                                caption = glue('(Shown are max. top 20 terms); {ud_caption} ; FDR < {alpha} ; GO-FDR < 0.01'),
                                theme = theme(plot.caption = element_markdown()))
        } else {
            ggplot() + geom_blank() + labs(title = 'GO term analysis:  No statistically significant enriched sets found', subtitle = subtit)
        }
    }

# Molecular signatures analysis
    mol_signatures_plot <- function(deseq_lfc, updown, lfcT, alpha, subtit, fname) {

        if (updown == 'up') {
                significant_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(log2FoldChange > lfcT, padj < alpha) |>
                                rownames_to_column(var = 'gene_symbol')
        } else if (updown == 'down') {
                significant_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(log2FoldChange < -lfcT, padj < alpha) |>
                                rownames_to_column(var = 'gene_symbol')
        }

        background_genes <- deseq_lfc |>
                                data.frame() |>
                                filter(baseMean != 0) |>
                                rownames_to_column(var = 'gene_symbol')

        # https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
        # H: hallmark gene sets
        # C1: positional gene sets
        # C2: curated gene sets
        # C3: motif gene sets
        # C4: computational gene sets
        # C5: GO gene sets
        # C6: oncogenic signatures
        # C7: immunologic signatures

        mol_sig <- tribble(
            ~mark, ~desc,
            'H',   'Hallmark genes',
            'C2',  'Curated online genes sets',
            'C5',  'Ontology genes',
            'C6',  'Oncogenic signature genes',
            'C7',  'Immunologic signature genes',
            'C8',  'Cell type signature genes'
        )

        p <- NULL

        for (i in 1:nrow(mol_sig)) {
            mark <- mol_sig[i,]$mark
            desc <- mol_sig[i,]$desc
            
            hallmark_genes <- msigdbr(species = 'Homo sapiens', category = mark) |>
                dplyr::select(gs_name, gene_symbol)


            em <- enricher(significant_genes$gene_symbol,
                           TERM2GENE = hallmark_genes, 
                           universe = background_genes$gene_symbol)

            if (!is.null(em)) {
                em@result |> write_csv(glue('{fname}.mol_sig_{mark}_results.{updown}.csv'))

                if (em@result |> filter(p.adjust < 0.05) |> nrow() > 0) {
                    p[[mark]] <- barplot(em, 
                                         title = glue('{mark}: {desc}'), 
                                         showCategory = 20, 
                                         font.size = 7)
                }
            }
        }

        if (!is.null(p)) {
            if (updown == 'up') {
                ud_title <- 'upregulated' 
                ud_caption <- glue('log<sub>2</sub>FC > {lfcT}')
            } else if (updown == 'down') {
                ud_title <- 'downregulated'
                ud_caption <- glue('log<sub>2</sub>FC < -{lfcT}')
            }

            wrap_plots(p) +
                plot_layout(axis_titles = 'collect') +
                plot_annotation(title = glue('Molecular signatures enrichment: {ud_title} genes'),
                                subtitle = subtit,
                                caption = glue('{ud_caption} ; FDR < {alpha} ; GO-FDR < 0.05'),
                                theme = theme(plot.caption = element_markdown()))
        } else {
            ggplot() + geom_blank() + labs(title = 'Molecular signatures enrichment:  No statistically significant enriched sets found', subtitle = subtit)
        }
    }

# GSEA
    gsea_analysis <- function(deseq_lfc, category = 'C5') {
        # category: H, C1-C8

        # Rank by shrunken log2 fold change ( alternatives are: p-value or sign(log2FoldChange) * -log10(pvalue))
        ranking <- deseq_lfc |>
                        data.frame() |>
                        filter(baseMean != 0) |>
                        rownames_to_column(var = 'gene_symbol') |>
                        mutate(rank = log2FoldChange) |>
                        arrange(desc(rank))


        gene_list <- ranking$rank
        names(gene_list) <- ranking$gene_symbol

        hallmark_genes <- msigdbr(species = 'Homo sapiens', category = category) |>
                            dplyr::select(gs_name, gene_symbol)


        gsea_results <- GSEA(gene_list, TERM2GENE = hallmark_genes, seed = 42)

        gsea_results
    }

    gsea_plot <- function(gsea_results, nplots = 30, category = 'C5', subtit, fname) {
        p <- NULL

        if (!is.null(gsea_results)) {
            gsea_results@result |> write_csv(glue('{fname}.GSEA_results.csv'))

            if (nrow(gsea_results@result) > 0) {

                # Plot up to top N categories
                for (i in seq_len(gsea_results@result |> head(n = nplots) |> nrow())) {
                    res <- gsea_results@result[i, ]

                    # Parse category names
                    id_split <- str_split_fixed(res$ID, '_', n =2)
                    if (nchar(id_split[, 1]) == 4) {
                        categ <- str_sub(id_split[, 1], 3, 4)
                    } else {
                        categ <- id_split[, 1]
                    }
                    tit <- id_split[, 2] |> str_replace_all('_', ' ') |> str_to_title()
                    
                    # Plot
                    p[[i]] <- gseaplot(gsea_results, 
                                       geneSetID = res$ID, 
                                       by = 'runningScore', 
                                       title = glue('{tit} ({categ})') |> str_wrap(width = 35) |> paste(collapse = '\n')) +
                                theme(axis.text.x = element_text(size = 10),
                                      axis.text.y = element_text(size = 10),
                                      axis.title = element_text(size = 10),
                                      plot.title = element_text(size = 10))
                }
            }
        }

        if (!is.null(p)) {
            wrap_plots(p) +
                plot_layout(axis_titles = 'collect') +
                plot_annotation(title = 'GSEA',
                                subtitle = subtit,
                                caption = glue('{category}; top {nplots} categories'))
        } else {
            ggplot() + geom_blank() + labs(title = 'GSEA:   No statistically significant enriched sets found', subtitle = subtit)
        }
    }

    gsea_ma_subplots <- function(gsea_results, deseq_lfc, lfcT, alpha, subtit, fname) {
        if (!is.null(gsea_results)) {
            if (nrow(gsea_results@result) > 0) {

                # Plot MA plot per category
                for (i in seq_len(gsea_results@result |> nrow())) {
                    res <- gsea_results@result[i, ]

                    # Parse category names
                    id_split <- str_split_fixed(res$ID, '_', n =2)
                    if (nchar(id_split[, 1]) == 4) {
                        categ <- str_sub(id_split[, 1], 3, 4)
                    } else {
                        categ <- id_split[, 1]
                    }
                    tit <- id_split[, 2] |> str_replace_all('_', ' ') |> str_to_title()
                    
                    # Get get leading edge genes per category
                    leading_edge_genes <- gsea_results@result[i,]$core_enrichment |> str_split(pattern = '/') |> unlist()

                    p <- deseq_lfc |>
                            data.frame() |>
                            rownames_to_column(var = 'gene_symbol') |>
                            filter(gene_symbol %in% leading_edge_genes) |>
                            ma_plot(lfcT = lfcT, alpha = alpha, subtit = subtit) +
                                labs(title = glue('{tit} ({categ})'))
                    

                    ggsave(glue('{fname}.GSEA.MAplot-{i}.png'), p, bg = 'white', width = 7, height = 4)

                }
            }
        }
    }



# Descriptive stats
## Number of assigned reads (sequencing depth)
### All
    data |>
        group_by(sample_id) |>
        summarize(assigned_reads = sum(read_count)) |>
        ggplot() +
            geom_bar(aes(x = sample_id, y = assigned_reads), stat = 'identity') +
            labs(title = 'Number of reads assigend to genes',
                 x = 'sample ID',
                 y = 'number of reads') +
            theme_half_open(font_size = 10) +
            background_grid(major = 'y') +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) + 
            scale_y_continuous(labels = scientific_10)

    ggsave('seq_depth_all.png', bg = 'white', width = 7, height = 4)

### Per catagory
    for (condition in conditions) {
        data |>
            group_by(sample_id) |>
            summarize(assigned_reads = sum(read_count)) |>
            left_join(annotation) |>
            ggplot() +
                geom_violin(aes(x = .data[[condition]], y = assigned_reads, fill = .data[[condition]])) +
                geom_boxplot(aes(x = .data[[condition]], y = assigned_reads), width = 0.1) +
                labs(title = 'Sequencing depth',
                     subtitle = 'Number of reads assigend to genes',
                     y = 'number of reads') +
                theme_half_open(font_size = 10) +
                background_grid(major = 'y') +
                scale_y_continuous(labels = scientific_10)

        ggsave(glue('seq_depth_{condition}.png'), bg = 'white', width = 7, height = 4)
    }

## Number of detected genes
    for (condition in conditions) {
        data |>
            filter(read_count > 0) |>
            group_by(sample_id) |>
            summarize(number_of_detected_genes = n()) |>
            left_join(annotation) |>
            ggplot() +
                geom_violin(aes(x = .data[[condition]], y = number_of_detected_genes, fill = .data[[condition]])) +
                geom_boxplot(aes(x = .data[[condition]], y = number_of_detected_genes), width = 0.1) +
                labs(title = 'Detected genes (>0 reads)',
                    y = 'number of genes') +
                theme_half_open(font_size = 10) +
                background_grid(major = 'y')

        ggsave(glue('detected_genes_{condition}.png'), bg = 'white', width = 7, height = 4)
    }

## Dependency No. reads and No. detected genes
### All
    data |>
        group_by(sample_id) |>
        summarize(total_reads = sum(read_count),
                  gene_count = sum(read_count > 0)) |>
        ggplot(aes(x = total_reads,
                   y = gene_count)) +
        geom_point() +
            geom_smooth(method = 'lm') +
            geom_text_repel(aes(label = sample_id), color = 'black', size = 2.5) +
            theme_half_open(font_size = 10) +
            scale_x_log10(labels = scientific_10) +
            scale_y_continuous(labels = scientific_10) +
            background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
            labs(title = 'Relationship between read number and detected genes',
                x = 'gene assigned reads (log<sub>10</sub>)',
                y = 'number of detected genes') +
            theme(axis.title.x = element_markdown())

    ggsave('depth_vs_detectedGenes_all.png', bg = 'white', width = 7, height = 4)

### Per category
    for (condition in conditions) {
        data |>
            group_by(sample_id) |>
            summarize(total_reads = sum(read_count),
                    gene_count = sum(read_count > 0)) |>
            left_join(annotation) |>
            ggplot(aes(x = total_reads,
                       y = gene_count,
                       color = .data[[condition]])) +
                geom_point() +
                geom_smooth(method = 'lm') +
                geom_text_repel(aes(label = sample_id), color = 'black', size = 2.5) +
                scale_x_log10(labels = scientific_10) +
                scale_y_continuous(labels = scientific_10) +
                theme_half_open(font_size = 10) +
                background_grid(major = 'xy', minor = 'xy', size.major = 0.3, size.minor = 0.1) +
                labs(title = 'Relationship between read number and detected genes',
                    x = 'gene assigned reads (log<sub>10</sub>)',
                    y = 'number of detected genes') +
                theme(axis.title.x = element_markdown())

        ggsave(glue('depth_vs_detectedGenes_{condition}.png'), bg = 'white', width = 7, height = 4)
    }


## Correlation
### Pearson
    data_wide |>
        ggcorr(limits = c(0,1), 
            midpoint = 0.5,
            label = FALSE,
            label_round = 2,
            label_color = 'gray30',
            label_size = 4,
            nudge_x = -5,
            size = 2) +
        labs(title = 'Pearson correlation')

    ggsave('correlation_pearson.png', bg = 'white', width = 7, height = 4)

### Spearman
    data_wide |>
        ggcorr(limits = c(0,1), 
            midpoint = 0.5,
            method = c('pairwise', 'spearman'),
            label = FALSE,
            label_round = 2,
            label_color = 'gray30',
            label_size = 4,
            nudge_x = -5,
            size = 2) +
        labs(title = 'Spearman correlation')

    ggsave('correlation_spearman.png', bg = 'white', width = 7, height = 4)


# PCA
## Normalize data with DESeq2
    library(DESeq2)

    countdata <- data_wide |>
                    as.matrix()
            
    data_dds <- 
        DESeqDataSetFromMatrix(
            countData = countdata,
            colData = annotation,
            design = formula_de)
        

    data_vsd <- vst(data_dds, blind = FALSE)    # nsub = 100, VST (variance stabilizing transformation)

## Plot PCA
    plot_pca(data_vsd, conditions, 'Before batch effects removal')
    ggsave('PCA_before_conditions.png', bg = 'white', width = 7, height = 4)

    if (!is.null(batches)) {
        plot_pca(data_vsd, batches, 'Before batch effects removal')
        ggsave('PCA_before_batches.png', bg = 'white', width = 7, height = 4)
    }

## Remove batch effects with Limma
    # TODO: Limma supports removal only 2 batch effects. For more batch effects variables they would need to be combined into one group variable
    if (!is.null(batches)) {
        library(limma)

        data_noBatch_vsd <- data_vsd
        mat <- assay(data_noBatch_vsd)
        mm <- model.matrix(~ group, colData(data_noBatch_vsd))

        if (length(batches) == 1) {
            mat <- limma::removeBatchEffect(mat, batch = data_noBatch_vsd[[batches]], design = mm)
        } else {
            mat <- limma::removeBatchEffect(mat, batch = data_noBatch_vsd[[batches[1]]], batch2 = data_noBatch_vsd[[batches[2]]], design = mm)
        }
        
        assay(data_noBatch_vsd) <- mat

        detach('package:limma', unload=TRUE)
    }

## Plot PCA
    if (!is.null(batches)) {
        plot_pca(data_noBatch_vsd, batches, 'After batch effects removal')
        ggsave('PCA_after_batches.png', bg = 'white', width = 7, height = 4)

        plot_pca(data_noBatch_vsd, conditions, 'After batch effects removal')
        ggsave('PCA_after_conditions.png', bg = 'white', width = 7, height = 4)
    }

# Hierarchical clustering
    if (!is.null(batches)) { 
        heat_vsd <- data_noBatch_vsd
    } else {
        heat_vsd <- data_vsd
    }
        
## Heatmap
    png(filename = 'clustering_heatmap.png', units = 'in', width = 5, height = 5, bg = 'white', res = 300)
    
    heat_vsd |>
        assay() |>
        t() |>
        dist() |>
        as.matrix() |> 
        heatmap(symm = TRUE)

    dev.off()

## Pheatmap
    sampleDists <- 
        heat_vsd |>
        assay() |>
        t() |>
        dist()

    sampleDistMatrix <- 
        sampleDists |>
        as.matrix()

    png(filename = 'clustering_pheatmap.png', units = 'in', width = 8, height = 7, bg = 'white', res = 100)
    pheatmap(sampleDistMatrix,
            clustering_distance_rows = sampleDists,
            clustering_distance_cols = sampleDists,
            annotation_col = colData(heat_vsd) |> data.frame() |> select(all_of(conditions)),
            border_color = NA)
    dev.off()

## Pheatmap of most variable genes
    heat_vsd |>
        assay() |>
        as.data.frame() |>
        rowwise() |>
        mutate(variance = var(c_across(everything()))) |>
        arrange(desc(variance)) |>
        head(n = 50) |>
        dplyr::select(-variance) -> mat

    png(filename = 'clustering_top50var_genes.png', units = 'in', width = 8, height = 7, bg = 'white', res = 100)
    pheatmap(mat, 
        annotation_col = colData(heat_vsd) |> data.frame() |> select(all_of(conditions)), 
        show_rownames = FALSE,
        scale = 'row',
        border_color = NA,
        main = 'Top 50 variable genes')
    dev.off()

# DE analysis
    library(ashr)

    data_de_dds <- DESeq(data_dds)


    for (i in 1:nrow(comparisons)) {
        comp1 <- comparisons[i, ]$comp1
        comp2 <- comparisons[i, ]$comp2

        # Change which comparison is going to be reference level and recalculate test
        data_de_dds$group <- relevel(data_de_dds$group, comp2)
        data_de_dds <- nbinomWaldTest(data_de_dds)

        # Calculate log2 fold change with DEseq2
        deseq_lfc <- lfcShrink(data_de_dds, contrast = c('group', comp1, comp2), type = 'ashr', lfcThreshold = 0)

        # Save results table
        deseq_lfc |>
            as.data.frame() |>
            rownames_to_column(var = 'gene_name') |>
            arrange(padj) |>
            write_csv(glue('{comp1}_vs_{comp2}.DE_results.csv'))

        # How many significantly changed genes do we have
        sig_genes_count <- deseq_lfc |>
                            data.frame() |>
                            filter(padj < alpha & abs(log2FoldChange) > lfcT) |>
                            nrow()

        # MA plot
        ma_plot(deseq_lfc, lfcT, alpha, subtit = glue('Comparing {comp1} vs. {comp2}'))
        ggsave(glue('{comp1}_vs_{comp2}.DE_MAplot.png'), bg = 'white', width = 7, height = 4)

        # Volcano plot
        volcano_plot(deseq_lfc, lfcT, alpha, subtit = glue('Comparing {comp1} vs. {comp2}'))
        ggsave(glue('{comp1}_vs_{comp2}.DE_volcano.png'), bg = 'white', width = 7, height = 4)


        # TODO: note that this works only for human because of 'org.Hs.eg.db' library
        library(msigdbr)
        library(fgsea)
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(enrichplot)

        if (sig_genes_count > 0) {

            # Plot expression of top DE genes
            gene_expression_plot(deseq_lfc, data_de_dds, lfcT, alpha, conditions, genesToPlot = 12, subtit = glue('{comp1} vs. {comp2}; ({condition})'))
            ggsave(glue('{comp1}_vs_{comp2}.DE_genes_{conditions[1]}.png'), bg = 'white', width = 10, height = 6)
            
            if (length(conditions) > 1) {
                # If there are two condition, create alternative plot whete second condition is mapped as color and first condition as facet 
                gene_expression_plot(deseq_lfc, data_de_dds, lfcT, alpha, c(conditions[2], conditions[1]), genesToPlot = 12, subtit = glue('{comp1} vs. {comp2}; ({condition})'))
                ggsave(glue('{comp1}_vs_{comp2}.DE_genes_{conditions[2]}.png'), bg = 'white', width = 10, height = 6)
            }

            # GO term enrichment analysis
            GO_dot_plot(deseq_lfc, updown = 'up', lfcT, alpha, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
            ggsave(glue('{comp1}_vs_{comp2}.GO_dotplot.up.png'), bg = 'white', width = 20, height = 10)
            GO_dot_plot(deseq_lfc, updown = 'down', lfcT, alpha, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
            ggsave(glue('{comp1}_vs_{comp2}.GO_dotplot.down.png'), bg = 'white', width = 20, height = 10)

            # Molecular signatures enrichment
            mol_signatures_plot(deseq_lfc, updown = 'up', lfcT, alpha, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
            ggsave(glue('{comp1}_vs_{comp2}.mol_sig.up.png'), bg = 'white', width = 18, height = 12)
            mol_signatures_plot(deseq_lfc, updown = 'down', lfcT, alpha, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
            ggsave(glue('{comp1}_vs_{comp2}.mol_sig.down.png'), bg = 'white', width = 18, height = 12)

        }

        # GSEA (Gene set enrichment analysis)
        gsea_results <- gsea_analysis(deseq_lfc, category = 'C5')
        
        gsea_plot(gsea_results, nplots = 30, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
        ggsave(glue('{comp1}_vs_{comp2}.GSEA.png'), bg = 'white', width = 18, height = 12)

        gsea_ma_subplots(gsea_results, deseq_lfc, lfcT, alpha, subtit = glue('{comp1} vs. {comp2}'), fname = glue('{comp1}_vs_{comp2}'))
    }
    







