from pipeline_util import *
from snakemake.shell import shell
# Hack to disable strict bash mode, because HOMER exit code is non-zero when no peaks are found
shell.prefix("")

######## Step: Peak Calling: HOMER ##################
rule all_homer_results:
    input:
        homer_peaks=expand(f'homer/{{sample}}_{config["homer_suffix"]}',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )

rule make_tags_dir:
    input: f"{config['bams_dir']}/{{sample}}.bam"
    output: temp(directory('homer/{sample}'))
    conda: '../envs/homer.env.yaml'
    shell: 'makeTagDirectory {output} {input}'

def homer_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = directory(f'homer/{control_sample}')

    return dict(
        signal=directory(f"homer/{sample}"),
        **control_args,
    )

rule call_peaks_homer:
    input: unpack(homer_input_fun)
    output:
        peaks=f'homer/{{sample}}_{config["homer_suffix"]}'
    log: f'logs/homer/{{sample}}.log'
    conda: '../envs/homer.env.yaml'
    params:
        work_dir=WORK_DIR,
        control_arg=lambda wildcards, input: f" -i {input.get('control', None)}" if input.get('control', None) else ""
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
        'echo'
