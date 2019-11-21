from pipeline_util import *

localrules: all_macs2_results

######## Step: Call peaks MACS2 ##################
rule all_macs2_results:
    input:
        macs2_peaks=expand(
            'macs2/{sample}_{macs2_suffix}_peaks.{macs2_mode}Peak',
            sample=fastq_aligned_names(config),
            macs2_mode=config['macs2_mode'],
            macs2_suffix=config['macs2_suffix'],
        )

def macs2_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control'] = f'bams/{control_sample}.bam'

    return dict(
        signal=f'bams/{sample}.bam',
        **control_args
    )


rule call_peaks_macs2:
    input: unpack(macs2_input_fun)
    output: f'macs2/{{sample}}_{config["macs2_suffix"]}_peaks.{config["macs2_mode"]}Peak'
    log: f'logs/macs2_{config["macs2_suffix"]}/{{sample}}_{config["macs2_suffix"]}_{config["macs2_mode"]}.log'
    conda: '../envs/py27.env.yaml'
    params:
        macs2_params=config['macs2_params'],
        macs2_suffix=config['macs2_suffix'],
        species=macs_species(config['genome']),
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else ""
    shell:
        'macs2 callpeak -t {input.signal} {params.control_arg} --outdir {params.outdir} ' 
        '-n {wildcards.sample}_{params.macs2_suffix} -g {params.species} ' 
        '{params.macs2_params} &> {log}'