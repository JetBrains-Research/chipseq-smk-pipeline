from pipeline_util import *

localrules: all_filter_aligned_reads_results, bam_filtered_multiqc

######## Step: Alignment QC ##################
rule all_filter_aligned_reads_results:
    input:
         bams=expand("bams/{sample}.bam", sample=fastq_aligned_names(FASTQ_PATHS)),
         multiqc_bam='multiqc/bam_filtered/multiqc.html',
         deduplicated=expand('deduplicated/{sample}.bam', sample=fastq_aligned_names(FASTQ_PATHS)),

# Filter aligned reads with good mapping score: (remove not aligned and bad mapped)
rule filter_sort_bam_single:
    input: '{anywhere}/{sample}.bam.raw'
    output: '{anywhere}/{sample}.bam'
    log: 'logs/bam_filtered/{anywhere}/{sample}.log'

    conda: '../envs/bio.env.yaml'
    shell:
        'samtools view -bh -q30 {input} > {output}.filtered 2> {log} &&'
        ' samtools sort {output}.filtered -o {output}  &> {log} &&'
        ' rm {output}.filtered  &>> {log}'

# Filtered reads qc
rule bam_filtered_stats:
    input: 'bams/{sample}.bam'
    output: 'qc/bam_filtered/samtools_stats/{sample}.txt'
    log: 'logs/bam_filtered/samtools_stats/{sample}.log'

    wrapper: '0.36.0/bio/samtools/stats'

rule bam_filtered_multiqc:
    input:
        expand(
            'qc/bam_filtered/samtools_stats/{sample}.txt',
            sample=fastq_aligned_names(FASTQ_PATHS)
        )
    output: 'multiqc/bam_filtered/multiqc.html'
    log: 'multiqc/bam_filtered/multiqc.log'

    wrapper: '0.36.0/bio/multiqc'


rule bam_deduplication:
    input: "bams/{sample}.bam"
    output:
          bam="deduplicated/{sample}.bam",
          metrics="qc/picard/{sample}.txt"
    log: "logs/bam_deduplicated/{sample}.log"

    threads: 4
    resources:
        threads = 4,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    params: "REMOVE_DUPLICATES=True"
    wrapper: "0.31.1/bio/picard/markduplicates"