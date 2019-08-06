from pipeline_util import *

localrules: step7_macs2_results

######## Step: Call peaks MACS2 ##################
rule step7_macs2_results:
    input:
        macs2_peaks=expand(
            'macs2/{sample}_{macs2_suffix}_peaks.{type}Peak',
            sample=fastq_aligned_names(config),
            type=config['macs2_mode'],
            macs2_suffix=config['macs2_suffix'],
        )

def macs2_input_fun(config):
    def inner(wildcards):
        sample = wildcards.sample

        args = {
            'signal': f'bams/{sample}.bam'
        }

        control_sample = sample_2_control(config)[sample]

        if control_sample:
            args['control'] = f'bams/{control_sample}.bam'

        return args

    return inner


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