from pipeline_util import *

ruleorder: bowtie2_align_paired > bowtie2_align_single
localrules: all_alignment_results, download_chrom_sizes, download_fa, bam_raw_multiqc

######## Step: Alignment QC ##################
rule all_alignment_results:
    input:
         bams=expand("bams/{sample}.bam", sample=fastq_aligned_names(config)),
         multiqc_bam_raw='multiqc/bam_raw/multiqc.html',

# Indexes:
rule download_chrom_sizes:
    output: f"{config['genome']}.chrom.sizes"
    log: f"logs/{config['genome']}.chrom.sizes.log"

    shell:
        'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes &> {log}'

rule download_fa:
    output: directory('fa')
    log: 'logs/fa.log'

    shell:
        "rsync -avz --partial --exclude='*.txt' " # To include single chromosome use: "rsync -avz --partial  --include='chr15*' --exclude='*' "
        'rsync://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/chromosomes/ {output} &> {log}'

rule bowtie2_index:
    input: 'fa'
    output: directory('bowtie2-index')
    log: 'logs/bam_raw/bowtie2/bowtie2-index.log'

    conda: '../envs/bio.env.yaml'
    params:
        files_list=lambda wildcards: ','.join(glob('fa/*.fa.gz')),
        target='bowtie2-index/{genome}'.format(genome=config['genome'])
    shell: 'mkdir -p {output} && bowtie2-build {params.files_list} {params.target} &> {log}'

# Align
rule bowtie2_align_single:
    input:
        sample=bowtie2_input_paths(config, False),
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam.raw")
    log: "logs/bam_raw/bowtie2/{sample}.log"

    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    params:
        index=lambda wildcards, input: os.path.join(
            str(input.bowtie2_index_path), config['genome']
        ),
        extra=''
    wrapper: "0.36.0/bio/bowtie2/align"

rule bowtie2_align_paired:
    input:
        sample=bowtie2_input_paths(config, True),
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam.raw")
    log: "logs/bam_raw/bowtie2/{sample}.log"

    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    conda: '../envs/bio.env.yaml'
    params:
        index=lambda wildcards, input: os.path.join(
            str(input.bowtie2_index_path), config['genome']
        ),
        extra=config["bowtie2_params"]
    wrapper: "0.36.0/bio/bowtie2/align"

# Aligned bams qc
rule bam_raw_stats:
    input: 'bams/{sample}.bam.raw'
    output: 'qc/bam_raw/samtools_stats/{sample}.txt'
    log: 'logs/bam_raw/samtools_stats/{sample}.log'
    wrapper: '0.36.0/bio/samtools/stats'

rule bam_raw_multiqc:
    input:
        expand(
            'logs/bam_raw/bowtie2/{sample}.log',
            sample=fastq_aligned_names(config)
        ),
        expand(
            'qc/bam_raw/samtools_stats/{sample}.txt',
            sample=fastq_aligned_names(config)
        )
    output: 'multiqc/bam_raw/multiqc.html'
    log: 'multiqc/bam_raw/multiqc.log'
    wrapper: '0.36.0/bio/multiqc'
