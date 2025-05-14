import os
from pipeline_util import *

localrules: all_hotspot_results

######## Step: Peak Calling: Hotspot ##################
rule all_hotspot_results:
    input:
        hotspot_peaks=expand(f'hotspot/{{sample}}.peak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )

rule call_peaks_hotspot:
    input:
        signal_pileup=f"{config['bams_dir']}/pileup/{{sample}}.bed"
    output: 'hotspot/{sample}.peak'
    log: 'logs/hotspot/{sample}.log'
    params:
        hotspot_executable=config['hotspot_executable'],
    resources:
        threads = 1,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        '{params.hotspot_executable} -i {input.signal_pileup} -o {wildcards.sample}.hotspot &> {log} &&'
        'cat {wildcards.sample}.hotspot | tail -n +2 |\
            awk -v OFS="\t" \'{{print $1, $6, $7, ".", $8}}\' |\
            sort -k1,1 -k2,2n -k3,3n > {output} && rm {wildcards.sample}.hotspot'
