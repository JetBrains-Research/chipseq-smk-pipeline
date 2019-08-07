import os
from pipeline_util import fastq_paths, trimmed_fastq_sample_names

ruleorder: trim_paired_fastq > trim_single_fastq
localrules: step2_trim_fastq_results, multiqc_trimmed_fastq

rule step2_trim_fastq_results:
    input:
        multiqc_fastq='multiqc/trimmed/multiqc.html'

rule trim_single_fastq:
    input: f"{config['fastq_dir']}/{{sample}}.{config['fastq_ext']}"
    output:
        "trimmed/{sample}_trimmed.fq.gz",
        f"trimmed/{{sample}}.{config['fastq_ext']}_trimming_report.txt"
    log: 'logs/trimmed/{sample}.log'

    threads: 4
    resources:
        threads=4,
        mem=16, mem_ram=12,
        time=60 * 120
    conda: '../envs/bio.env.yaml'
    # params:
    #     extra = "--cores {threads}" # "--illumina -q 20"
    # wrapper: '0.36.0/bio/trim_galore/se' # only 0.4.3 trim galore version
    conda: '../envs/bio.env.yaml'
    params:
        extra=lambda wildcards, threads: f"--gzip --cores {threads}",  # "--illumina -q 20"
        out_dir="trimmed"
    shell:
        "trim_galore {params.extra} -o {params.out_dir} {input} &> {log}"
        # " && mv trimmed/{wildcards.sample}"
        # f"_trimmed.{trim_galore_file_suffix()}"
        # " {output} &>> {log}"

# MultiQC handles "_trimmed" suffix, so it could combine cutadapt (trimgalore)
# and fastqc reports correctly for single-end reads. Let's rename paired-end
# reads to *_1_trimmed.fq.gz, *_2_trimmed.fq.gz
rule trim_paired_fastq:
    input:
        first=f"{config['fastq_dir']}/{{sample}}_1.{config['fastq_ext']}",
        second=f"{config['fastq_dir']}/{{sample}}_2.{config['fastq_ext']}"
    output:
        fq1 = f"trimmed/{{sample}}_1_trimmed.fq.gz",
        fq1_rep = f"trimmed/{{sample}}_1.{config['fastq_ext']}_trimming_report.txt",
        fq2 = f"trimmed/{{sample}}_2_trimmed.fq.gz",
        fq2_rep = f"trimmed/{{sample}}_2.{config['fastq_ext']}_trimming_report.txt"
    log: 'logs/trimmed/{sample}.log'

    threads: 4
    resources:
        threads=4,
        mem=16, mem_ram=12,
        time=60 * 120
    # wrapper: '0.36.0/bio/trim_galore/pe' only 0.4.3 trim galore version
    params:
        extra=lambda wildcards, threads: f"--gzip --cores {threads}",  # "--illumina -q 20"
        out_dir="trimmed"
    conda: '../envs/bio.env.yaml'
    shell:
        "trim_galore {params.extra} --paired -o {params.out_dir} {input} &> {log} &&"
        " mv trimmed/{wildcards.sample}_1_val_1.fq.gz {output.fq1} &>> {log} &&"
        " mv trimmed/{wildcards.sample}_2_val_2.fq.gz {output.fq2} &>> {log}"

rule trimmed_fastqc:
    input: f"trimmed/{{any}}.fq.gz"
    output:
        html='qc/trimmed/fastqc/{any}_fastqc.html',
        zip='qc/trimmed/fastqc/{any}_fastqc.zip'
    log: 'logs/trimmed/fastqc/{any}.log'
    wildcard_constraints: any="[^/]+"

    resources:
        threads=1,
        mem=8, mem_ram=4,
        time=60 * 120

    # https://bitbucket.org/snakemake/snakemake-wrappers/src/0.31.1/bio/fastqc/
    wrapper: '0.36.0/bio/fastqc'


rule multiqc_trimmed_fastq:
    input:
        expand(
            'qc/trimmed/fastqc/{trimmed_sample}_fastqc.zip',
            trimmed_sample=trimmed_fastq_sample_names(config)
        ),
        expand(
            'trimmed/{fastq_file}_trimming_report.txt',
            fastq_file=[os.path.basename(p) for p in fastq_paths(config)]
        ),

    output: 'multiqc/trimmed/multiqc.html'
    log: 'multiqc/trimmed/multiqc.log'

    wrapper: '0.36.0/bio/multiqc'
