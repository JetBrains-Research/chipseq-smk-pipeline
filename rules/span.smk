from pipeline_util import *

localrules: download_span

######## Step: Peak Calling: SPAN ##################
rule all_span:
    input:
        span_peaks=expand(f'span/{{sample}}_{config["span_bin"]}_{config["span_fdr"]}_{config["span_gap"]}.peak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        ) if bool(config['span']) else []


rule download_span:
    output: 'bin/span-1.1.5628.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-1.1.5628.jar'


def span_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = SAMPLE_2_CONTROL_MAP[sample]
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
        peaks=f'span/{{sample}}_{{bin}}_{{fdr}}_{{gap}}.peak'
    log: f'logs/span/{{sample}}_{{bin}}_{{fdr}}_{{gap}}.log'
    conda: '../envs/java.env.yaml'
    threads: 4
    params:
        span_params=config['span_params'],
        span_iterations=config['span_iterations'],
        span_threshold=config['span_threshold'],
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else ""
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
         'java -Xmx{resources.mem_ram}G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
         '{params.control_arg} --peaks {output.peaks} --model span/fit/{wildcards.sample}_{wildcards.bin}.span '
         '--workdir span --iterations {params.span_iterations} --threshold {params.span_threshold} --threads {threads} '
         '--bin {wildcards.bin} --fdr {wildcards.fdr} --gap {wildcards.gap} {params.span_params} &> {log}'
