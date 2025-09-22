from pipeline_util import *

localrules: download_span

######## Step: Peak Calling: SPAN ##################
rule all_span_results:
    input:
        span_peaks=expand(f'span/{{sample}}_{config["span_bin"]}_{config["span_fdr"]}.peak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


rule download_span:
    output: 'bin/span.jar'
    shell: 'wget -O {output}  https://download.jetbrains.com/biolabs/span/span-1.6.6499.jar'


def span_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        **control_args,
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_span:
    input: unpack(span_input_fun)
    output:
        peaks='span/{sample}_{bin}_{fdr}.peak'
    log: 'logs/span/{sample}_{bin}_{fdr}.log'
    conda: '../envs/java.env.yaml'
    threads: config['span_threads']
    params:
        span_params=config['span_params'],
        span_fragment=config['span_fragment'],
        span_iterations=config['span_iterations'],
        span_threshold=config['span_threshold'],
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else "",
        additional_arg=lambda wildcards, input: "",
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
         'java -Xmx{resources.mem_ram}G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
         '{params.control_arg} --peaks {output.peaks} --model span/{wildcards.sample}_{wildcards.bin}.span '
         '--workdir span --iterations {params.span_iterations} --threshold {params.span_threshold} '
         '--bin {wildcards.bin} --fragment {params.span_fragment} --fdr {wildcards.fdr} --threads {threads} '
         '{params.span_params} {params.additional_arg} &> {log}'
