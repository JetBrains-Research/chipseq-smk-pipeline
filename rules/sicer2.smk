import os

from pipeline_util import *

localrules: all_sicer_results

######## Step: Peak Calling: SICER2 ##################
def sicer2_all_peaks_input():
    files = []

    window_size=config['sicer2_window']
    gap=config['sicer2_gap']

    # XXX: change significance only via config, SICER rule takes the value from
    # config, not via wildcards
    # IMPORTANT
    # We explicitly check presence of output files, otherwise missing intermediate tmp files
    # are causing even existing files re-computations.
    for sample in filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS)):
        if sample not in SAMPLE_2_CONTROL_MAP:
            # w/o control
            significance = config['sicer2_evalue']
            no_control_file = f'sicer2/{sample}-W{window_size}-G{gap}-E{significance}.scoreisland'
            if not os.path.exists(f'{WORK_DIR}/{no_control_file}'):  # Workaround for intermediate tmp
                files.append(no_control_file)
        else:
            # with control
            significance = config['sicer2_fdr']
            with_control_file = f'sicer2/{sample}-W{window_size}-G{gap}-FDR{significance}-island.bed'
            if not os.path.exists(f'{WORK_DIR}/{with_control_file}'):  # Workaround for intermediate tmp
                files.append(with_control_file)

    return files

rule all_sicer2_results:
    input: *sicer2_all_peaks_input()


def sicer2_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control_pileup'] = f"{config['bams_dir']}/pileup/{control_sample}.bed"

    return dict(
        signal_pileup = f"{config['bams_dir']}/pileup/{{sample}}.bed",
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_sicer2:
    input: unpack(sicer2_input_fun)
    output: 'sicer2/{sample}-W{width}-G{gap, \d+}-{any_suffix}'
    log: 'logs/sicer2/{sample}-W{width}-G{gap}-{any_suffix}.log'
    conda: '../envs/sicer2.env.yaml'
    shadow: "shallow"
    params:
        significance=lambda wildcards, input: \
            f"-fdr {config['sicer2_fdr']}" if input.get('control_pileup', None) else f"-e {config['sicer2_evalue']}",
        control_arg=lambda wildcards, input: \
            f"-c {WORK_DIR}/{input.control_pileup}" if input.get('control_pileup', None) else "",
        peaks_file=lambda wildcards, output: os.path.basename(output[0]),
        workdir=WORK_DIR,
        fragment=config['sicer2_fragment'],
        genome=config['genome'],
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        # sicer -t treatment.bed -c control.bed -s hg38 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000
        # sicer -t treatment.bed -s hg38
        'echo "Significance threshold: {params.significance}" > {params.workdir}/{log} &&'
        ' tmp_sicer=$(mktemp -d) && mkdir -p $tmp_sicer && cd $tmp_sicer && '
        ' sicer -t {params.workdir}/{input.signal_pileup} {params.control_arg}'
        ' -s {params.genome} -w {wildcards.width} -rt 1 -f {params.fragment} '
        ' -g {wildcards.gap} {params.significance} -cpu {threads} &>> {params.workdir}/{log} &&'
        ' ls -lah  &>> {params.workdir}/{log} &&'
        ' mv {params.peaks_file} {params.workdir}/{output} &>> {params.workdir}/{log} && '
        ' rm -rf $tmp_sicer'
