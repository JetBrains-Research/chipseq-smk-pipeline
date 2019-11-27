from pipeline_util import *

localrules: all_raw_qc_results, multiqc_fastq

####### Step: RAW Reads QC ##################
rule all_raw_qc_results:
    input:
         multiqc_fastq='multiqc/fastqc/multiqc.html',

rule fastqc:
    input: f"{config['fastq_dir']}/{{sample}}.{config['fastq_ext']}"
    output:
          html='qc/fastqc/{sample}_fastqc.html',
          zip='qc/fastqc/{sample}_fastqc.zip'
    log: 'logs/fastqc/{sample}.log'

    resources:
        threads = 1,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    wrapper: '0.36.0/bio/fastqc' # https://bitbucket.org/snakemake/snakemake-wrappers/src/0.31.1/bio/fastqc/

rule multiqc_fastq:
    input: expand('qc/fastqc/{sample}_fastqc.zip', sample=fastq_names_wo_ext(FASTQ_PATHS))
    output: 'multiqc/fastqc/multiqc.html'
    log: 'multiqc/fastqc/multiqc.log'

    wrapper: '0.36.0/bio/multiqc'
