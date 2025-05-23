import os
from pipeline_util import *

localrules: all_bayespeak_results

######## Step: Peak Calling: BayesPeak ##################

rule all_bayespeak_results:
    input:
        bayespeak_peaks=expand(f'bayespeak/{{sample}}.bed',
            sample=filter(lambda f: not is_control(f) and f in SAMPLE_2_CONTROL_MAP,
                aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


def bayespeak_input_fun(wildcards):
    sample = wildcards.sample
    control_sample = SAMPLE_2_CONTROL_MAP[sample]

    return dict(
        signal_pileup=f"{config['bams_dir']}/pileup/{sample}.filtered.bed",
        control_pileup=f"{config['bams_dir']}/pileup/{control_sample}.filtered.bed",
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_bayespeak:
    input: unpack(bayespeak_input_fun)
    output: 'bayespeak/{sample}.bed'
    log: 'logs/bayespeak/{sample}.log'
    params:
        bayespeak_rscript_executable=config['bayespeak_rscript_executable'],
        control_pileup=lambda wildcards, input: input['control_pileup']
    resources:
        threads = 1,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        '{params.bayespeak_rscript_executable} {workflow.basedir}/scripts/run_bayespeak.R \
        {input.signal_pileup} {input.control_pileup} {wildcards.sample}.csv &> {log} && '
        'mv {wildcards.sample}.csv {output}'
