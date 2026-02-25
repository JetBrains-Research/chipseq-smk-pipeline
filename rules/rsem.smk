from pipeline_util import *

rule rsem_all:
    input: expand('rsem-quant/{sample}.genes.results', sample=get_rnaseq_samples())

######## RSEM index ##################
rule rsem_index:
    input:
        fa='fa',
        gtf=f"{config['genome']}.gtf"
    output:
        grp=f'rsem-index/{config["genome"]}.grp',
        ti=f'rsem-index/{config["genome"]}.ti',
        seq=f'rsem-index/{config["genome"]}.seq',
        idx_fa=f'rsem-index/{config["genome"]}.idx.fa',
        n2g=f'rsem-index/{config["genome"]}.n2g.idx.fa'
    log: 'logs/rsem-index.log'
    conda: '../envs/rnaseq.env.yaml'
    threads: 8
    resources:
        mem=32, mem_ram=28,
        time=60 * 180
    params:
        fa_files=lambda wildcards: ' '.join(glob('fa/*.fa.gz')),
        prefix=f'rsem-index/{config["genome"]}'
    shell:
        'mkdir -p rsem-index && '
        'gunzip -c {params.fa_files} > rsem-index/genome.fa && '
        'rsem-prepare-reference '
        '--gtf {input.gtf} '
        '--num-threads {threads} '
        'rsem-index/genome.fa '
        '{params.prefix} &> {log} && '
        'rm rsem-index/genome.fa'

######## RSEM quantification ##################
rule rsem_quantify:
    input:
        bam='star-aligned/{sample}.toTranscriptome.out.bam',
        rsem_idx=f'rsem-index/{config["genome"]}.grp'
    output:
        genes='rsem-quant/{sample}.genes.results',
        isoforms='rsem-quant/{sample}.isoforms.results'
    log: 'logs/rsem-quant/{sample}.log'
    conda: '../envs/rnaseq.env.yaml'
    threads: 8
    resources:
        mem=32, mem_ram=28,
        time=60 * 180
    params:
        prefix='rsem-quant/{sample}',
        rsem_ref=f'rsem-index/{config["genome"]}'
    shell:
        # --paired-end --bam --estimate-rspd --no-bam-output
        'rsem-calculate-expression '
        '--paired-end '
        '--bam '
        '--no-bam-output '
        '--num-threads {threads} '
        '--estimate-rspd '
        '{input.bam} '
        '{params.rsem_ref} '
        '{params.prefix} &> {log}'
