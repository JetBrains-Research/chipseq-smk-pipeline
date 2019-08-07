from pipeline_util import *

localrules: step9_span, step9_span_tuned, download_span

######## Step: Peak Calling: SPAN ##################
rule step9_span:
    input:
        span_peaks=expand('span/{sample}_{span_bin}_{span_fdr}_{span_gap}.peak',
            sample=fastq_aligned_names(config),
            span_bin=config['span_bin'],
            span_fdr=config['span_fdr'],
            span_gap=config['span_gap']
        ),

rule step9_span_tuned:
    input:
        span_tuned_peaks=tuned_peaks_input_files(config)

rule download_span:
    output: 'bin/span-0.11.0.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.11.0.4882.jar'


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
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control',
            None) else ""
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        'java -Xmx8G -jar {input.span} analyze -t {input.signal} --chrom.sizes {'
        'input.chrom_sizes} '
        '--peaks {output.peaks} --model {output.model} --workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {params.fdr} --gap {params.gap} {params.span_params} &> {log}'


def span_tuned_input_fun(wildcards):
    args = dict(
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output
    )

    sample = wildcards.sample
    bin = wildcards.bin

    anns_file = find_labels_for_sample(sample, config)
    assert anns_file, f"Peaks annotations file is missing for {sample}"
    args['span_markup'] = anns_file
    args['model'] = f'span/fit/{sample}_{bin}.span'
    return args

rule call_peaks_span_tuned:
    input:
        unpack(span_tuned_input_fun)
    output: 'span/{sample}_{bin}_tuned.peak'
    log: 'logs/span/{sample}_{bin}_tuned.log'

    conda: 'envs/java8.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        'java -Xmx8G -jar bin/span-0.11.0.jar analyze --model {input.model} '
        '--workdir span --threads {threads}  --labels {input.span_markup} --peaks {output} &> {log}'