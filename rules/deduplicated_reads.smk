from pipeline_util import *

localrules: all_deduplicated_reads_results

######## Step: Deduplication ##################
rule all_deduplicated_reads_results:
    input:
        deduplicated=expand('deduplicated/{sample}.bam',sample=fastq_aligned_names(config,FASTQ_PATHS)),

rule bam_deduplication:
    input: "bams/{sample}.bam"
    output:
        bam="deduplicated/{sample}.bam",
        metrics="qc/picard/{sample}.txt"
    log: "logs/bam_deduplicated/{sample}.log"

    threads: 2
    resources:
        threads=2,
        mem=8,mem_ram=4,
        time=60 * 120
    params: "REMOVE_DUPLICATES=True"
    wrapper: "0.31.1/bio/picard/markduplicates"