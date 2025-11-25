#!/usr/bin/env nextflow

// Script parameters
params.samples = "./sample_config.csv"
params.outdir = "results"
params.hisat2_index = "/data/genome/chm13v2"
params.gene_annotation = "/data/genome/chm13.draft_v2.0.gene_annotation.gff3"
params.annotation_gene_name = "gene_name"         // GFF field for common gene name
params.annotation_gene_id = "source_gene"         // GFF field for gene id
params.annotation_gene_biotype = "gene_biotype"   // GFF field for biotype of the gene
params.check_biotype = true                       // Should check for reads aligning to e.g. rRNA or mtRNA be run
params.reads_strand = "reverse"                   // forward | reverse | no
params.conditions = "timePoint,longCovid"
params.batches = null                             //"seqBatch,prepBatch"
params.reflevels = "T0,N"
params.log2fc = "0.5"
params.pval = "0.05"

// Print config summary
log.info """\
    RNA-seq   pipeline
 ===================================
 input reads                 : ${params.samples}
 outdir                      : ${params.outdir}
 hisat2 index                : ${params.hisat2_index}
 gene annotation             : ${params.gene_annotation}
 """


// Define processes
process TRIM_ADAPT {
    /*
    Trim adapter seqeunce. 
    Filter reads that are less than 10 nt long.
    */
    label 'small'
    tag "$sample"
    conda 'conda-forge::python=3.9.13 conda-forge::python-isal=0.11.1 bioconda::cutadapt=4.1'

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("${sample}_R*.t.fastq")
    path "*.log"

    script:
    """
    cutadapt \
        -a AGATCGGAAGAGC \
        -A AGATCGGAAGAGC \
        --minimum-length 10 \
        --cores=$task.cpus \
        -o ${sample}_R1.t.fastq \
        -p ${sample}_R2.t.fastq \
        ${files[0]} ${files[1]} \
        >> ${sample}.cutadapt.log 2>&1
    """
}


process ALIGN {
    /*
    Align reads to chm13v2 human genome using hisat2
    */
    label 'small'
    tag "$sample"
    conda 'bioconda::hisat2'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.sam")
    path "*.log"

    script:
    """
    hisat2 \
        -p ${task.cpus} \
        -x ${params.hisat2_index} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -S ${sample}.sam \
        >> ${sample}.hisat2.log 2>&1
    """
}


process SORT_SAM {
    /*
    Sort .sam files and convert to .bam format
    */
    label 'small'
    tag "$sample"
    conda 'bioconda::samtools=1.15.1'
    publishDir "${params.outdir}/bam", mode:'copy'

    input:
    tuple val(sample), path(sam_file)

    output:
    tuple val(sample), path("*.bam")

    script:
    """
    # -q 2
    samtools view -q 2 -@ ${task.cpus} -h ${sam_file} \
        | samtools sort -@ ${task.cpus} -o ${sample}.bam -
    """
}


process HTSEQ_GENES {
    /*
    Count reads towards all genes (exon overlap only).
    */
    label 'single'
    tag "$sample"
    conda 'bioconda::htseq'
    publishDir "${params.outdir}/htseq", mode:'copy'

    input:
    tuple val(sample), path(bam_file)

    output:
    val(sample)
    path("*.genes.tsv")

    script:
    if ( params.reads_strand == "reverse" ) {
        strand = "reverse"
    }
    else if (params.reads_strand == "forward") {
        strand = "yes"
    }
    else {
        strand = "no"
    }

    def gene_name = "${params.annotation_gene_name}"
    def gene_id = "${params.annotation_gene_id}"
    def gene_biotype = "${params.annotation_gene_biotype}"

    """
    htseq-count \
        -f bam \
        -t exon \
        -i ${gene_name} \
        --additional-attr ${gene_id} \
        --additional-attr ${gene_biotype} \
        -m intersection-strict \
        -r pos \
        --stranded ${strand} \
        -c ${sample}_htseq.genes.tsv \
        ${bam_file} \
        ${params.gene_annotation}
    """
}


process HTSEQ_BIOTYPES {
    /*
    Count reads towards all genes (exon overlap only).
    */
    label 'single'
    tag "$sample"
    conda 'bioconda::htseq'

    input:
    tuple val(sample), path(bam_file)

    output:
    tuple val(sample), path("*.biotypes.tsv")

    script:
    if ( params.reads_strand == "reverse" ) {
        strand = "reverse"
    }
    else if (params.reads_strand == "forward") {
        strand = "yes"
    }
    else {
        strand = "no"
    }

    def gene_name = "${params.annotation_gene_name}"
    def gene_id = "${params.annotation_gene_id}"
    def gene_biotype = "${params.annotation_gene_biotype}"

    """
    htseq-count \
        -f bam \
        -t exon \
        -i ${gene_id} \
        --additional-attr ${gene_name} \
        --additional-attr ${gene_biotype} \
        -m intersection-strict \
        --nonunique random \
        --minaqual 0 \
        -r pos \
        --stranded ${strand} \
        -c ${sample}_htseq.biotypes.tsv \
        ${bam_file} \
        ${params.gene_annotation}
    """
}


process SUMMARIZE_BIOTYPES {
    /*
    Create summary of how many reads were assigned to each gene biotype.
    */
    label 'single'
    tag "$sample"
    publishDir "${params.outdir}/biotypes", mode:'copy'

    input:
    tuple val(sample), path(htseq_file)

    output:
     path "*.biotypes.tsv"

    shell:
    '''
    awk -v OFS="\t" \
        '$1 !~ /__/ {count[$3] += $4}
         END {for (key in count) {
                print key, count[key]}}' \
        !{htseq_file} > !{sample}.biotypes.tsv
    
    '''
}


process CALC_NORM {
    /*
    Calculate normalization factors for all samples
    */
    label 'single'
    tag "all"
    conda 'conda-forge::r-base conda-forge::r-optparse conda-forge::r-tidyverse bioconda::bioconductor-deseq2'
    publishDir "${params.outdir}/norm_factors", mode:'copy'

    input:
    val(names)
    path(files)

    output:
    path("norm_factors.tsv")

    script:
    """
    norm_calc.R --inFiles ${files.join(',')} \
                --sampleNames ${names.join(',')}
    """
}


process DE_GENE_ANALYSIS {
    /*
    Perform differential gene expression analysis
    */
    label 'single'
    tag "all"
    conda 'conda-forge::r-base conda-forge::r-optparse conda-forge::r-tidyverse conda-forge::r-patchwork conda-forge::r-ggrepel conda-forge::r-ggbeeswarm conda-forge::r-ggtext conda-forge::r-pheatmap bioconda::bioconductor-deseq2 bioconda::bioconductor-limma conda-forge::r-ashr conda-forge::r-msigdbr bioconda::bioconductor-fgsea bioconda::bioconductor-clusterprofiler bioconda::bioconductor-org.hs.eg.db bioconda::bioconductor-enrichplot conda-forge::r-cowplot conda-forge::r-ggally'
    publishDir "${params.outdir}/DE_stats", mode:'copy'

    input:
    path(files)

    output:
    path("*.csv")
    path("*.png")

    script:
    def conds = "${params.conditions}"
    def batches = params.batches ? "--batches ${params.batches}" : ""
    def sample_annot = "${params.samples}"
    def reference_levels = "${params.reflevels}"
    def log2fc = "${params.log2fc}"
    def pval = "${params.pval}"
    """
    de_gene_analysis.R --inFiles ${files.join(',')} \
                       --outDir \$(pwd) \
                       --annotation ${sample_annot} \
                       --conditions ${conds} \
                       ${batches} \
                       --reflevels ${reference_levels} \
                       --log2fc ${log2fc} \
                       --pval ${pval}
    """
}



process MAKE_CHROM_SIZES {
    /*
    Get chromosome sizes from alignment samfile
    */
    label 'single'
    tag "$sample"

    input:
    tuple val(sample), path(sam_file)

    output:
    path "chm13v2.chrom.sizes"

    shell:
    '''
    # Extract chromosome lengths
    head -n 1000 !{sam_file} \
        | awk -v OFS="\t" '$1 ~ /^@SQ/ {split($2, chr, ":")
                                        split($3, size, ":")
                                        print chr[2], size[2]}' > chm13v2.chrom.sizes
    '''
}

process CREATE_BIGWIG_TRACKS {
    /*
    Create .bigwig tracks of aligned data
    */
    label 'single'
    tag "$sample"
    conda 'bioconda::bedtools bioconda::ucsc-bedgraphtobigwig'
    publishDir "${params.outdir}/bigwig", mode:'copy'

    input:
    tuple val(sample), path(bam_file)
    path(chrom_sizes)
    path(norm_factors)

    output:
    path "*.bigWig"

    script:
    """
    # Get scaling factor
    normVal=\$(awk -v sam=${sample} '\$1 == sam {print \$2}' ${norm_factors})

    # Calculate read coverage
    bedtools genomecov -bg \
                       -split \
                       -ibam ${bam_file} \
        | LC_COLLATE=C sort -k1,1 -k2,2n \
        | awk -v OFS='\t' -v nomr=\${normVal} \
            '{\$4 = nomr*\$4
              print \$0}'  \
        > ${sample}.bdg

    # Convert read coverage to .bigWig
    bedGraphToBigWig \
            ${sample}.bdg \
            ${chrom_sizes} \
            ${sample}.bigWig
    """
}


process CREATE_STRANDED_TRACKS {
    /*
    Create coverage tracks for each strand
    */
    tag "$sample"
    conda 'bioconda::star bioconda::ucsc-bedgraphtobigwig conda-forge::parallel'
    publishDir "${params.outdir}/bigwig-stranded", mode:'copy'
    cpus 2
    memory '20 GB'

    input:
    tuple val(sample), path(bam_file)
    path(chrom_sizes)
    path(norm_factors)

    output:
    path "*.bigWig"

    script:
    if ( params.reads_strand == "reverse" ) {
        first  = "str2"
        second = "str1"
    } else {
        first  = "str1"
        second = "str2"
    }
    """
    # Get scaling factor
    normVal=\$(awk -v sam=${sample} '\$1 == sam {print \$2}' ${norm_factors})

    # Calculate read coverage
    STAR \
        --runMode inputAlignmentsFromBAM \
            --inputBAMfile ${bam_file} \
            --outWigType bedGraph \
            --outWigNorm None \
            --outWigStrand Stranded \
            --outFileNamePrefix ./${sample}_

    # Sort and scale
    awk -v OFS='\t' -v nomr=\${normVal} \
                    '{\$4 = nomr*\$4
                      print \$0}'  \
                    ${sample}_Signal.Unique.${first}.out.bg \
                   | LC_COLLATE=C sort -k1,1 -k2,2n \
                   > ${sample}.pos.bdg

    awk -v OFS='\t' -v nomr=\${normVal} \
                    '{\$4 = "-"nomr*\$4
                      print \$0}'  \
                    ${sample}_Signal.Unique.${second}.out.bg \
                   | LC_COLLATE=C sort -k1,1 -k2,2n \
                   > ${sample}.neg.bdg

    # Convert read coverage to .bigWig
    parallel -j 2 "bedGraphToBigWig \
                    ${sample}.{1}.bdg \
                    ${chrom_sizes} \
                    ${sample}.{1}.bigWig" \
        ::: pos neg
   """
}


process FASTQC {
    tag "$sample"
    conda 'bioconda::fastqc'
    publishDir "${params.outdir}/fastqc", mode:'copy'
    cpus 5
    memory '300 GB'

    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}.fastqc.log"

    script:
    """
    mkdir ${sample}.fastqc.log
    fastqc --outdir ${sample}.fastqc.log -f fastq -t $task.cpus -q ${reads[0]}
    """
}


process MULTIQC {
    label 'single'
    tag "Generating final stats"
    conda 'bioconda::multiqc'
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc --config ${baseDir}/multiqc_config.yml .
    """
}


workflow {
    /* ch_files = channel.fromFilePairs( params.samples, size: 2, checkIfExists: true) */

    // Trim and align reads
    ch_files = Channel.fromPath(params.samples, checkIfExists: true) \
                | splitCsv(header:true, sep:',') \
                | map { row -> tuple(row.sampleId, [file(row.read1), file(row.read2)])}

    (ch_trimmed, ch_trimmed_log) = TRIM_ADAPT( ch_files )
    (ch_aligned, ch_aligned_log) = ALIGN( ch_trimmed )

    // Filter, sort, deduplicate
    ch_sorted_bam = SORT_SAM( ch_aligned )
    //(ch_dedup_bam, ch_dedup_log) = REMOVE_DUP( ch_sorted_bam )

    // Quality control
    ch_fastqc = FASTQC( ch_trimmed )

    // Count reads with htseq-count
    (ch_sample_name, ch_htseq_gene) = HTSEQ_GENES( ch_sorted_bam )

    // Biotype quality check
    if ( params.check_biotype ) {
        ch_htseq_biotype = HTSEQ_BIOTYPES( ch_sorted_bam )
        ch_biotype_summary = SUMMARIZE_BIOTYPES( ch_htseq_biotype )
    } else {
        ch_biotype_summary = Channel.empty()
    }

    // Run differential gene analysis
    DE_GENE_ANALYSIS( ch_htseq_gene.collect() )

    // Calculate normalization factors and chromosome sizes
    ch_norm_fact = CALC_NORM( ch_sample_name.collect(), ch_htseq_gene.collect() )
    ch_chrom_sizes = MAKE_CHROM_SIZES( ch_aligned.first() )

    // Create tracks from total RNA
    CREATE_BIGWIG_TRACKS( ch_sorted_bam, ch_chrom_sizes, ch_norm_fact )
    CREATE_STRANDED_TRACKS( ch_sorted_bam, ch_chrom_sizes, ch_norm_fact )

    // Create summary
    MULTIQC(ch_fastqc.mix(ch_trimmed_log, ch_aligned_log, ch_htseq_gene, ch_biotype_summary).collect(flat: true))

}

workflow.onComplete {
    log.info ( workflow.success 
        ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" 
        : "Oops .. something went wrong" )

    def gitCommit =  workflow.commitId ? workflow.commitId : "git -C ${workflow.projectDir} rev-parse --short HEAD".execute().text.trim()
    def gitRepo = workflow.repository ? workflow.repository : "git -C ${workflow.projectDir} config --get remote.origin.url".execute().text.trim()
    def gitRevision = workflow.revision ? workflow.revision : "git -C ${workflow.projectDir} symbolic-ref --short HEAD".execute().text.trim()

    def meta = file("${params.outdir}/run_info.txt")
    meta.text = """
    ==============
    Pipeline info
    ==============
    Nextflow version: ${workflow.nextflow.version} - ${workflow.nextflow.build}
    Pipeline repository: ${gitRepo}
    Pipeline revision: ${gitRevision}
    Commit ID: ${gitCommit}
    Run name: ${workflow.runName}
    Date: ${workflow.complete}

    ===============
    Run parameters
    ===============
    config file: ${params.samples}
    outdir: ${params.outdir}
    hisat2 index: ${params.hisat2_index}
    genome annotation: ${params.gene_annotation}
        field - gene name: ${params.annotation_gene_name} 
        field - gene id: ${params.annotation_gene_id}
        field - gene biotype: ${params.annotation_gene_biotype}
    library strandness: ${params.reads_strand}
    DESeq2:
        experimental conditions: ${params.conditions}
        batch effects: ${params.batches}
        reference conditions levels: ${params.reflevels}
        log2 fold-change cutoff: ${params.log2fc}
        p-value cutoff: ${params.pval}
    run biotype check: ${params.check_biotype}
    """

    // Make a copy of the config file for a future reference
    def source = file(params.samples).toFile()
    def dest   = new File("${params.outdir}/samples_config.csv")

    dest.bytes = source.bytes
}

// nextflow run RNA-seq_pipeline.nf -with-report -with-trace -with-timeline -with-dag dag.png

// nextflow run RNA-seq_pipeline.nf -profile cluster --samples sample_config.csv

