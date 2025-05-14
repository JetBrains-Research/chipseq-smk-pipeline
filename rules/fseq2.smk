from pipeline_util import *


######## Step: Peak Calling: Fseq2 ##################
rule all_fseq2_results:
    input:
        fseq2_peaks=expand(f'fseq2/{{sample}}_peaks.narrowPeak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


def fseq2_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_fseq2:
    input: unpack(fseq2_input_fun)
    output:
        peaks=f'fseq2/{{sample}}_peaks.narrowPeak'
    log: f'logs/fseq2/{{sample}}.log'
    conda: '../envs/fseq2.env.yaml'
    params:
        work_dir=WORK_DIR,
        control_arg=lambda wildcards, input: \
            f" -control_file {WORK_DIR}/{input.control}" if input.get('control', None) else ""
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
        'fseq2 callpeak -v -q_thr 0.05 {params.control_arg} -chrom_size_file {input.chrom_sizes} '
        '-name {wildcards.sample} {params.work_dir}/{input.signal} &> {log}'
