from pipeline_util import *

localrules: all_macs3_results

######## Step: Call peaks MACS3 ##################
rule all_macs3_results:
    input:
        macs3_peaks=expand(
            f'macs3/{{sample}}_{config["macs3_suffix"]}_peaks.{config["macs3_mode"]}Peak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )

def macs3_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        **control_args
    )


rule call_peaks_macs3:
    input: unpack(macs3_input_fun)
    output: f'macs3/{{sample}}_{config["macs3_suffix"]}_peaks.{config["macs3_mode"]}Peak'
    log: f'logs/macs3_{config["macs3_suffix"]}/{{sample}}_{config["macs3_suffix"]}_{config["macs3_mode"]}.log'
    conda: '../envs/macs3.env.yaml'
    params:
        macs3_params=config['macs3_params'],
        macs3_suffix=config['macs3_suffix'],
        species=macs_species(config['genome']),
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
        control_arg=lambda wildcards, input: f" -c {input.control}" if input.get('control', None) else ""
    shell:
        'macs3 callpeak -f BAM -t {input.signal} {params.control_arg} --outdir {params.outdir} ' 
        '-n {wildcards.sample}_{params.macs3_suffix} -g {params.species} ' 
        '{params.macs3_params} &> {log} || true; if [[ ! -f {output} ]]; then touch {output}; fi;'