from pipeline_util import *

localrules: download_omnipeak

######## Step: Peak Calling: Omnipeak ##################
rule all_omnipeak_results:
    input:
        omnipeak_peaks=expand(f'omnipeak/{{sample}}_{config["omnipeak_bin"]}_{config["omnipeak_fdr"]}.peak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


rule download_omnipeak:
    output: 'bin/omnipeak.jar'
    shell: 'wget -O {output}  https://download.jetbrains.com/biolabs/omnipeak/omnipeak-1.0.6677.jar'


def omnipeak_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        **control_args,
        omnipeak=rules.download_omnipeak.output,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_omnipeak:
    input: unpack(omnipeak_input_fun)
    output:
        peaks='omnipeak/{sample}_{bin}_{fdr}.peak'
    log: 'logs/omnipeak/{sample}_{bin}_{fdr}.log'
    conda: '../envs/java.env.yaml'
    threads: config['omnipeak_threads']
    params:
        omnipeak_params=config['omnipeak_params'],
        omnipeak_fragment=config['omnipeak_fragment'],
        omnipeak_iterations=config['omnipeak_iterations'],
        omnipeak_threshold=config['omnipeak_threshold'],
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else "",
        additional_arg=lambda wildcards, input: "",
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
         'java --add-modules=jdk.incubator.vector -Xmx{resources.mem_ram}G -jar {input.omnipeak} '
         'analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} {params.control_arg} '
         '--peaks {output.peaks} --model omnipeak/{wildcards.sample}_{wildcards.bin}.omni '
         '--workdir omnipeak --iterations {params.omnipeak_iterations} --threshold {params.omnipeak_threshold} '
         '--bin {wildcards.bin} --fragment {params.omnipeak_fragment} --fdr {wildcards.fdr} --threads {threads} '
         '{params.omnipeak_params} {params.additional_arg} &> {log}'
