# rna-seq-pipeline
Differential expression analysis of RNA-seq data

Usage:

    nextflow run RNA-seq_pipeline.nf -profile cluster --samples sample_config.csv [params]

### Parameters
    --samples                      - path to sample configuration file in csv format. Header: sampleId,read1,read2,condition[,condition2,..][,batch1,..]
    --outdir                       - output directory for results
    --hisat2_index                 - path to hisat2 genome index files
    --gene_annotation              - path to gene annotation in gff3 format
    --annotation_gene_name         - gff3 field name for common gene name
    --annotation_gene_id           - gff3 field name for gene id
    --annotation_gene_biotype>     - gff3 field name for gene biotype
    --check_biotype                - summarize gene biotypes of annotated reads (default: true)
    --reads_strand                 - Reads mapping orientation towards genes [forward, reverse, no] (default: reverse)
    --conditions                   - Comma-separated list of treatment condition names from csv file header e.g. mutation_type
    --batches                      - Comma-separated list of names of sample batches from csv file header e.g. libraryprep_batch
    --reflevels                    - Values of treatmet that are considered refrence. One per experimental condition. e.g. heathy,T0
    --log2fc                       - Significance threshold for log2 fold change of gene expression (default: 0.5)
    --pval                         - Threshold for FDR p-value (default: 0.05)


