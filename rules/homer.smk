from pipeline_util import *


######## Step: Peak Calling: HOMER ##################
rule all_homer_results:
    input:
        homer_peaks=expand(f'homer/{{sample}}.peaks',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


def homer_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        **control_args,
    )


rule call_peaks_homer:
    input: unpack(homer_input_fun)
    output:
        peaks=f'homer/{{sample}}.peaks'
    log: f'logs/homer/{{sample}}.log'
    conda: '../envs/homer.env.yaml'
    params:
        work_dir=WORK_DIR,
        control_mktags = lambda wildcards, input: \
            f"makeTagDirectory control {WORK_DIR}/{input.control} && " if input.get('control', None) else "",
        control_arg=lambda wildcards, input: f" -i control" if input.get('control', None) else ""
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
        'tmp_hotspot=$(mktemp -d) && mkdir -p $tmp_hotspot && cd $tmp_hotspot &&'
        'makeTagDirectory treatment {params.work_dir}/{input.signal} && {params.control_mktags} '
        'findPeaks treatment -style histone -o auto {params.control_arg} && '
        'cat treatment/regions.txt | grep -v "#" | cut -f2- | sort -k1,1 -k2,2n -k3,3n > {params.work_dir}/{output.peaks} && '
        'cd .. && rm -rf $tmp_hotspot'
