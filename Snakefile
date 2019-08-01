from pipeline_util import *

MACS2_PARAMS = '-q 0.05'
MACS2_PARAMS_ATAC_SEQ = '-q 0.05 -f BAMPE --nomodel --nolambda -B --call-summits'
MACS2_PARAMS_BROAD = '--broad --broad-cutoff 0.05'

print('CONFIG\n{}'.format('\n'.join(['{}: {}'.format(k, v) for k, v in config.items()])))

workdir: config['work_dir']

localrules: download_chrom_sizes, download_fa, multiqc_fastq, download_phantompeakqualtools, multiqc_bowtie2, bam_stats_multiqc, download_span

rule download_chrom_sizes:
    output: '{}.chrom.sizes'.format(config['genome'])
    shell:
         'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes'

rule download_fa:
    output: directory('fa')
    shell: "rsync -avz --partial --exclude='*.txt' "
           'rsync://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/chromosomes/ {output} && '
           'gunzip -f {output}/*.fa.gz'

rule bowtie2_index:
    input: directory('fa')
    output: directory('bowtie2-index')
    params:
          files_list=lambda wildcards: ','.join(glob('fa/*.fa')),
          target='bowtie2-index/{genome}'.format(genome=config['genome'])
    conda: 'envs/bio.env.yaml'
    resources:
        threads = 1,
        mem = 16, mem_ram = 10,
        time = 60 * 120
    shell: 'mkdir -p {output} && bowtie2-build {params.files_list} {params.target}'

rule fastqc:
    input: os.path.join(config['fastq_dir'], '{sample}.fastq')
    output:
          html='qc/fastqc/{sample}_fastqc.html',
          zip='qc/fastqc/{sample}_fastqc.zip'
    log: 'logs/fastqc/{sample}.log'
    resources:
        threads = 1,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    wrapper: '0.31.1/bio/fastqc'

rule multiqc_fastq:
    input: expand('qc/fastqc/{sample}_fastqc.zip', sample=fastq_names(config))
    output: 'multiqc/fastqc/multiqc.html'
    log: 'multiqc/fastqc/multiqc.log'
    wrapper: '0.31.1/bio/multiqc'

rule cleaned_fastqc:
    input: 'cleaned/{sample}.fastq'
    output:
          html='qc/cleaned/fastqc/{sample}_fastqc.html',
          zip='qc/cleaned/fastqc/{sample}_fastqc.zip'
    log: 'logs/cleaned/fastqc/{sample}.log'
    resources:
        threads = 1,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    wrapper: '0.31.1/bio/fastqc'

rule cleaned_multiqc_fastq:
    input: expand('qc/cleaned/fastqc/{sample}_fastqc.zip', sample=fastq_names(config))
    output: 'multiqc/cleaned/fastqc/multiqc.html'
    log: 'multiqc/cleaned/fastqc/multiqc.log'
    wrapper: '0.31.1/bio/multiqc'

rule trim_single_fastq:
    input: os.path.join(config['fastq_dir'], '{sample}.fastq')
    output: 'cleaned/{sample}_se.fastq'
    threads: 4
    log: 'logs/trim/{sample}.log'
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    conda: 'envs/bio.env.yaml'
    shell: 'trim_galore --cores {threads} {input} -o cleaned/ &> {log}; '
           'mv cleaned/{wildcards.sample}_trimmed.fq {output}'

rule trim_paired_fastq:
    input:
         first=os.path.join(config['fastq_dir'], '{sample}_1.fastq'),
         second=os.path.join(config['fastq_dir'], '{sample}_2.fastq')
    output:
         first='cleaned/{sample}_pe_1.fastq',
         second='cleaned/{sample}_pe_2.fastq'
    threads: 4
    log: 'logs/trim/{sample}.log'
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    conda: 'envs/bio.env.yaml'
    shell: 'trim_galore --cores {threads} --paired {input.first} {input.second} -o cleaned/ &> {log}; '
           'mv cleaned/{wildcards.sample}_1_val_1.fq {output.first}; mv cleaned/{wildcards.sample}_2_val_2.fq {output.second}'

rule bowtie2_align_single:
    input:
        sample=["cleaned/{sample}_se.fastq"],
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
          index=lambda wildcards, input : os.path.join(
              config['work_dir'], str(input.bowtie2_index_path), config['genome']
          ),
          extra=''
    wrapper: "0.31.1/bio/bowtie2/align"

rule bowtie2_align_paired:
    input:
        sample=["cleaned/{sample}_pe_1.fastq", "cleaned/{sample}_pe_2.fastq"],
        bowtie2_index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}.bam2.raw")
    log: "logs/bowtie2/{sample}.log"
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    conda: 'envs/bio.env.yaml'
    params:
          index=lambda wildcards, input : os.path.join(
              config['work_dir'], str(input.bowtie2_index_path), config['genome']
          ),
          extra="-X 2000 --dovetail"    
    wrapper: "0.31.1/bio/bowtie2/align"

rule filter_sort_bam_single:
    input: '{anywhere}/{sample}.bam.raw'
    output: '{anywhere}/{sample}.bam'
    conda: 'envs/bio.env.yaml'
    resources:
        threads = 2,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    shell: 'samtools view -bh -q30 {input} > {output}.filtered; samtools sort {output}.filtered -o {output}; rm {output}.filtered'

rule filter_sort_bam_paired:
    input: '{anywhere}/{sample}.bam2.raw'
    output: '{anywhere}/{sample}.bam'
    conda: 'envs/bio.env.yaml'
    resources:
        threads = 2,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    shell: 'samtools view -bh -f2 -q30 {input} > {output}.filtered; samtools sort {output}.filtered -o {output}; rm {output}.filtered'

rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample, [^/]*}.bam.bai'
    wrapper: '0.31.1/bio/samtools/index'

rule remove_unmapped:
    input: "bams/{sample}.bam"
    output: "mapped/{sample}.bam"
    conda: 'envs/bio.env.yaml'
    resources:
        threads = 2,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    shell: 'samtools view -b -F 4 {input} > {output}'

rule download_phantompeakqualtools:
    output: directory('bin/phantompeakqualtools')
    params:
          targz='phantompeakqualtools.tar.gz'
    shell: 'cd bin; '
           'curl --location '
           'https://storage.googleapis.com/google-code-archive-downloads/v2/'
           'code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz '
           '--output {params.targz}; '
           'tar xvf {params.targz}'

# This rule requires spp R package installed
rule bam_qc_phantom:
    input:
         ppqt_dir=rules.download_phantompeakqualtools.output,
         bam='bams/{sample}.bam'
    output: 'qc/phantom/{sample}.phantom.tsv'
    params:
          run_spp=lambda wildcards, input: os.path.join(str(input.ppqt_dir), 'run_spp.R')
    shell: 'Rscript {params.run_spp} -c={input.bam} -savp -out={output} -rf'

rule bam_to_pileup:
    input: 'bams/{sample}.bam'
    output: temp('bams/pileup/{sample}.bed')
    conda: 'envs/bio.env.yaml'
    shell: 'bedtools bamtobed -i {input} > {output}'

rule bam_qc_pbc_nrf:
    input: rules.bam_to_pileup.output
    output: 'qc/pbc_nrf/{sample}.pbc_nrf.tsv'
    params:
          tmp_dir='tmp'
    shell: '''
mkdir -p {params.tmp_dir} &&
(T=$'\\t'
>&2 echo "TotalReadPairs${{T}}DistinctReadPairs${{T}}OneReadPair${{T}}TwoReadPairs${{T}}\
NRF=Distinct/Total${{T}}PBC1=OnePair/Distinct${{T}}PBC2=OnePair/TwoPair"

cat {input} | \
    sort -k1,1 -k3,3n -k2,2n -k6,6 -T {params.tmp_dir} | \
    awk -v OFS='\\t' '{{print $1,$2,$3,$6}}' | uniq -c | \
    awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}}
    ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}}
    END{{
        if (mt!=0){{m0_t=m0/mt}} else {{m0_t=-1.0}};
        if (m0!=0){{m1_0=m1/m0}} else {{m1_0=-1.0}};
        if (m2!=0){{m1_2=m1/m2}} else {{m1_2=-1.0}};
        printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0_t,m1_0,m1_2;
    }}') > {output}
    '''

rule multiqc_bowtie2:
    input: expand('logs/bowtie2/{sample}.log', sample=fastq_aligned_names(config))
    output: 'multiqc/bowtie2/multiqc.html'
    log: 'multiqc/bowtie2/multiqc.log'
    wrapper: '0.31.1/bio/multiqc'

rule bam_stats:
    input: 'bams/{sample}.bam'
    output: 'qc/samtools_stats/{sample}_samtools_stats.txt'
    wrapper: '0.31.1/bio/samtools/stats'

rule bam_stats_multiqc:
    input: expand('qc/samtools_stats/{sample}_samtools_stats.txt', sample=fastq_aligned_names(config))
    output: 'multiqc/samtools_stats/multiqc.html'
    log: 'multiqc/samtools_stats/multiqc.log'
    wrapper: '0.31.1/bio/multiqc'

rule bam2bw:
    input:
         bam='bams/{filename}.bam',
         bai='bams/{filename}.bam.bai'
    output: 'bw/{filename, [^/]*}.bw'
    conda: 'envs/deeptools.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output}'

rule call_peaks_macs2:
    input: unpack(macs2_input_fun(config))
    output: 'macs2/{sample}_{macs2_suffix}_peaks.{type}Peak'
    log: 'logs/macs2/macs2_{macs2_suffix}/{sample}_{macs2_suffix}_macs2_{type}.log'

    conda: 'envs/py27.env.yaml'
    params:
        macs2_params=config['macs2_params'],
        species=macs_species(config['genome']),
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control',
            None) else ""
    shell:
        'macs2 callpeak -t {input.signal} {params.control_arg} --outdir {params.outdir} ' 
        '-n {wildcards.sample}_{wildcards.macs2_suffix} -g {params.species} ' 
        '{params.macs2_params} &> {log}'


rule pileup_bed_effective_genome_fraction:
    input:
        pileup_bed=rules.bam_to_pileup.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        temp(str(rules.bam_to_pileup.output) + ".egf")
    run:
        value = effective_genome_fraction(
            config['genome'], input.chrom_sizes, input.pileup_bed
        )
        shell("echo -n '{value}' > {output}")


# Access to rules required
def sicer_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control_pileup'] = f'bams/pileup/{control_sample}.bed'

    return dict(
        # pileup_bed='bams/pileup/{sample}.bed',
        signal_pileup = f'bams/pileup/{sample}.bed',
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
        effective_genome_fraction=rules.pileup_bed_effective_genome_fraction.output,
    )

rule call_peaks_sicer:
    input: unpack(sicer_input_fun)
    output: 'sicer/{sample}-W{width}-G{gap, \d+}-{any_suffix}'
    log: 'logs/sicer/{sample}-W{width}-G{gap}-{any_suffix}.log'

    conda: 'envs/py27.env.yaml'
    shadow: "shallow"
    params:
        significance=lambda wildcards, input: config['sicer_fdr'] if input.get('control_pileup',
            None) else config['sicer_evalue'],
        signal_pileup_bed_fname=lambda wildcards, input: os.path.basename(input.signal_pileup),
        control_arg=lambda wildcards, input: os.path.basename(input.control_pileup) if input.get('control_pileup', None) else "",
        pileups_dir=lambda wildcards, input: os.path.split(str(input.signal_pileup))[0],
        peaks_file=lambda wildcards, output: os.path.basename(output[0]),
        fragment=config['sicer_fragment'],
        genome=config['genome'],
        script=lambda wildcards, input: "SICER.sh" if input.get('control_pileup', None) else "SICER-rb.sh"
    shell:
        # SICER.sh ["InputDir"] ["bed file"] ["control file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] [“FDR"]
        #
        # SICER-rb.sh ["InputDir"] ["bed file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] ["E-value"]
        'echo "Significance threshold: {params.significance}" > {log} &&'
        ' mkdir -p tmp_sicer &&'
        ' cd tmp_sicer && '
        '  {params.script} ../{params.pileups_dir} {params.signal_pileup_bed_fname} {params.control_arg}'
        '    $(pwd) {params.genome} 1 {wildcards.width}'
        '    {params.fragment} $(cat "../{input.effective_genome_fraction}")'
        '    {wildcards.gap} {params.significance} &>> ../{log} &&'
        ' ls -lah  &>> ../{log} &&'
        ' mv {params.peaks_file} ../{output} &>> ../{log}'


rule download_span:
    output: 'bin/span-0.11.0.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.11.0.4882.jar'

# Access to rules required
def span_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control'] = f'bams/{control_sample}.bam'

    return dict(
        signal=f'bams/{sample}.bam',
        **control_args,
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_span:
    input: unpack(span_input_fun)
    output:
        peaks=f'span/{{sample}}_{{bin}}_{config["span_fdr"]}_{config["span_gap"]}.peak',
        model='span/fit/{sample}_{bin}.span'
    log: f'logs/span/{{sample}}_{{bin}}{config["span_fdr"]}_{config["span_gap"]}.log'

    conda: 'envs/java8.env.yaml'
    threads: 4
    params:
        fdr = config["span_fdr"],
        gap = config["span_gap"],
        span_params=config['span_params'],
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else ""
    shell:
        'java -Xmx3G -jar {input.span} analyze -t {input.signal} --chrom.sizes {'
        'input.chrom_sizes} '
        '--peaks {output.peaks} --model {output.model} --workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {params.fdr} --gap {params.gap} {params.span_params} &> {log}'


rule all:
    input:
         multiqc_fastq='multiqc/fastqc/multiqc.html',

         cleaned_multiqc_fastq='multiqc/cleaned/fastqc/multiqc.html',

         multiqc_bowtie2='multiqc/bowtie2/multiqc.html',

         multiqc_samtools_stats='multiqc/samtools_stats/multiqc.html',

         bws=expand('bw/{sample}.bw',
                    sample=fastq_aligned_names(config)),

         mapped=expand('mapped/{sample}.bam',
                       sample=fastq_aligned_names(config)),

         bam_qc_phantom=expand('qc/phantom/{sample}.phantom.tsv',
                               sample=fastq_aligned_names(config)),

         bam_qc_pbc=expand('qc/pbc_nrf/{sample}.pbc_nrf.tsv',
                           sample=fastq_aligned_names(config)),

         macs2_peaks=expand('macs2/{sample}_{macs2_suffix}_peaks.{type}Peak',
                            sample=fastq_aligned_names(config),
                            type=config['macs2_mode'],
                            macs2_suffix=config['macs2_suffix']),

         sicer_peaks=sicer_all_peaks_input(config),

         span_peaks=expand('span/{sample}_{span_bin}_{span_fdr}_{span_gap}.peak',
                           sample=fastq_aligned_names(config),
                           span_bin=config['span_bin'],
                           span_fdr=config['span_fdr'],
                           span_gap=config['span_gap'])
