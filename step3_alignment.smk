from pipeline_util import *

ruleorder: bowtie2_align_paired > bowtie2_align_single
localrules: step3_alignment_results, download_chrom_sizes, download_fa, bam_raw_multiqc

######## Step: Alignment QC ##################
rule step3_alignment_results:
    input:
         bams=expand("bams/{sample}.bam", sample=fastq_aligned_names(config)),
         multiqc_bam_raw='multiqc/bam_raw/multiqc.html',
         multiqc_bam='multiqc/bam_filtered/multiqc.html',

# Indexes:
rule download_chrom_sizes:
    output: '{}.chrom.sizes'.format(config['genome'])
    shell:
        'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes'

rule download_fa:
    output: directory('fa')
    shell:
        "rsync -avz --partial --exclude='*.txt' " # To include single chromosome use: "rsync -avz --partial  --include='chr15*' --exclude='*' "
        'rsync://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/chromosomes/ {output}'

rule bowtie2_index:
    input: 'fa'
    output: directory('bowtie2-index')
    conda: 'envs/bio.env.yaml'

    params:
        files_list=lambda wildcards: ','.join(glob('fa/*.fa.gz')),
        target='bowtie2-index/{genome}'.format(genome=config['genome'])
    shell: 'mkdir -p {output} && bowtie2-build {params.files_list} {params.target}'

# Align
rule bowtie2_align_single:
    input:
        sample=[("trimmed" if is_trimmed(config) else config['fastq_dir']) + f"/{{sample}}.{config['fastq_ext']}"],
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam.raw")
    log: "logs/bowtie2/{sample}.log"

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

def bowtie2_paired_input_paths(config):
    if is_trimmed(config):
        return [
            f"trimmed/{{sample}}_1_trimmed.{trim_galore_file_suffix(config)}",
            f"trimmed/{{sample}}_2_trimmed.{trim_galore_file_suffix(config)}"
        ]
    else:
        return [
            config['fastq_dir'] + f"/{{sample}}_1.{config['fastq_ext']}",
            config['fastq_dir'] + f"/{{sample}}_2.{config['fastq_ext']}",
        ]

rule bowtie2_align_paired:
    input:
        sample=bowtie2_paired_input_paths(config),
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam.raw")
    log: "logs/bowtie2/{sample}.log"

    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    conda: 'envs/bio.env.yaml'
    params:
        index=lambda wildcards, input: os.path.join(
            str(input.bowtie2_index_path), config['genome']
        ),
        extra="-X 2000 --dovetail"
    wrapper: "0.36.0/bio/bowtie2/align"

# Aligned bams qc
rule bam_raw_stats:
    input: 'bams/{sample}.bam.raw'
    output: 'qc/bam_raw_samtools_stats/{sample}.txt'
    wrapper: '0.36.0/bio/samtools/stats'

rule bam_raw_multiqc:
    input:
        expand(
            'logs/bowtie2/{sample}.log',
            sample=fastq_aligned_names(config)
        ),
        expand(
            'qc/bam_raw_samtools_stats/{sample}.txt',
            sample=fastq_aligned_names(config)
        )
    output: 'multiqc/bam_raw/multiqc.html'
    log: 'multiqc/bam_raw/multiqc.log'
    wrapper: '0.36.0/bio/multiqc'
