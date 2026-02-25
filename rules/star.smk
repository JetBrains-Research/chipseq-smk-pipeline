from pipeline_util import *

localrules: download_gtf

rule star_all:
    input: expand('star-aligned/{sample}.toTranscriptome.out.bam', sample=get_rnaseq_samples())

######## GTF download ##################
rule download_gtf:
    output: f"{config['genome']}.gtf"
    log: f'logs/{config["genome"]}.gtf.log'
    shell:
        'wget -O {output}.gz http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/genes/refGene.gtf.gz &> {log} && '
        'gunzip {output}.gz'


######## STAR index ##################
rule star_index:
    input:
        fa='fa',
        gtf=f"{config['genome']}.gtf"
    output: directory('star-index')
    log: 'logs/star-index.log'
    conda: '../envs/rnaseq.env.yaml'
    threads: 8
    resources:
        threads=8,
        mem=64, mem_ram=60,
        time=60 * 180
    shell:
        'mkdir -p {output} {output}/tmp_fa && '
        'for f in {input.fa}/*.fa.gz; do '
        '  gunzip -c "$f" > {output}/tmp_fa/$(basename "$f" .gz); '
        'done && '
        'FA_FILES=$(find {output}/tmp_fa -name "*.fa" | sort | tr "\\n" " ") && '
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles $FA_FILES '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100 &> {log}; '
        'rm -rf {output}/tmp_fa'

######## STAR alignment ##################
rule star_align:
    input:
        sample = paired_reads_paths(config,True),
        star_index='star-index',
        gtf=f"{config['genome']}.gtf"
    output:
        bam='star-aligned/{sample}.sortedByCoord.out.bam',
        transcriptome_bam='star-aligned/{sample}.toTranscriptome.out.bam',
    log: 'logs/star-align/{sample}.log'
    # Workaround for broken conda STAR on macosx
    # conda: '../envs/rnaseq.env.yaml'
    params:
        r1=lambda wildcards, input: input.sample[0],
        r2=lambda wildcards, input: input.sample[1],
        prefix='star-aligned/{sample}/',
        workdir=WORK_DIR,
    threads: 16
    resources:
        threads=16,
        mem=64, mem_ram=60,
        time=60 * 180
    shell:
        'mkdir -p {params.prefix} && '
        'R1={params.r1} && '
        'R2={params.r2} && '
        'if [[ "$R1" == *.gz ]]; then '
        '  gunzip -c "$R1" > {params.prefix}R1.fastq && R1={params.prefix}R1.fastq; '
        'fi && '
        'if [[ "$R2" == *.gz ]]; then '
        '  gunzip -c "$R2" > {params.prefix}R2.fastq && R2={params.prefix}R2.fastq; '
        'fi && '
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.workdir}/{input.star_index} '
        '--readFilesIn {params.workdir}/"$R1" {params.workdir}/"$R2" '
        '--outFileNamePrefix {params.workdir}/{params.prefix} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outSAMunmapped Within '
        '--outSAMattributes Standard '
        '--quantMode TranscriptomeSAM &> {log} && '
        'mv {params.workdir}/{params.prefix}Aligned.sortedByCoord.out.bam {output.bam} && '
        'mv {params.workdir}/{params.prefix}Aligned.toTranscriptome.out.bam {output.transcriptome_bam} && '
        'rm -rf {params.workdir}/{params.prefix}'