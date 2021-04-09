from pipeline_util import *

localrules: all_span, all_span_tuned, download_span

######## Step: Peak Calling: SPAN ##################
rule all_span:
    input:
        span_peaks=expand(f'span/{{sample}}_{config["span_bin"]}_{config["span_fdr"]}_{config["span_gap"]}.peak',
            sample=fastq_aligned_names(FASTQ_PATHS)
        )

rule all_span_tuned:
    input:
        span_tuned_peaks=tuned_peaks_input_files(config, FASTQ_PATHS)

rule download_span:
    output: 'bin/span-0.13.5244.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.13.5244.jar'


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
    conda: '../envs/java8.env.yaml'
    threads: 4
    params:
        span_params=config['span_params'],
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else ""
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
         'java -Xmx8G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
         '{params.control_arg} --peaks {output.peaks} --model span/fit/{wildcards.sample}_{wildcards.bin}.span '
         '--workdir span --threads {threads} '
         '--bin {wildcards.bin} --fdr {wildcards.fdr} --gap {wildcards.gap} {params.span_params} &> {log}'


def span_tuned_input_fun(config):
    """Return inner function here to be launched during execution"""
    def inner(wildcards):
        args = dict(
            span=rules.download_span.output,
            chrom_sizes=rules.download_chrom_sizes.output
        )

        sample = wildcards.sample
        anns_file = find_labels_for_sample(sample, config)
        assert anns_file, f"Peaks annotations file is missing for {sample}"
        args['span_markup'] = anns_file
        return args
    return inner

rule call_peaks_span_tuned:
    input:
        unpack(span_tuned_input_fun(config))
    output: 'span/{sample}_{bin}_tuned.peak'
    log: 'logs/span/{sample}_{bin}_tuned.log'

    conda: '../envs/java8.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
         'java -Xmx8G -jar {input.span} analyze --model span/fit/{wildcards.sample}_{wildcards.bin}.span '
         '--workdir span --threads {threads}  --labels {input.span_markup} --peaks {output} &> {log}'