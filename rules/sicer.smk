import os
from pipeline_util import *

localrules: all_sicer_results

######## Step: Peak Calling: SICER ##################
def sicer_all_peaks_input():
    files = []

    window_size=config['sicer_window']
    gap=config['sicer_gap']

    # XXX: change significance only via config, SICER rule takes the value from
    # config, not via wildcards
    # IMPORTANT
    # We explicitly check presence of output files, otherwise missing intermediate tmp files
    # are causing even existing files re-computations.
    for sample in filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS)):
        if sample not in SAMPLE_2_CONTROL_MAP:
            # w/o control
            significance = config['sicer_evalue']
            no_control_file = f'sicer/{sample}-W{window_size}-G{gap}-E{significance}.scoreisland'
            if not os.path.exists(f'{WORK_DIR}/{no_control_file}'):  # Workaround for intermediate tmp
                files.append(no_control_file)
        else:
            # with control
            significance = config['sicer_fdr']
            with_control_file = f'sicer/{sample}-W{window_size}-G{gap}-islands-summary-FDR{significance}'
            if not os.path.exists(f'{WORK_DIR}/{with_control_file}'):  # Workaround for intermediate tmp
                files.append(with_control_file)

    return files

rule all_sicer_results:
    input: *sicer_all_peaks_input()


rule bam_to_pileup:
    input: f"{config['bams_dir']}/{{sample}}.bam"
    output: temp(f"{config['bams_dir']}/pileup/{{sample}}.bed")
    conda: '../envs/bio.env.yaml'
    shell: 'bedtools bamtobed -i {input} > {output}'


rule pileup_bed_effective_genome_fraction:
    input:
        pileup_bed=rules.bam_to_pileup.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        temp(str(rules.bam_to_pileup.output) + ".egf")
    run:
        value = effective_genome_fraction(
            config['genome'], input.chrom_sizes, input.pileup_bed
        )
        shell("echo -n '{value}' > {output}")


def sicer_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control_pileup'] = f"{config['bams_dir']}/pileup/{control_sample}.bed"

    return dict(
        signal_pileup = f"{config['bams_dir']}/pileup/{{sample}}.bed",
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
        effective_genome_fraction=rules.pileup_bed_effective_genome_fraction.output,
    )


rule call_peaks_sicer:
    input: unpack(sicer_input_fun)
    output: 'sicer/{sample}-W{width}-G{gap, \d+}-{any_suffix}'
    log: 'logs/sicer/{sample}-W{width}-G{gap}-{any_suffix}.log'
    conda: '../envs/py27.env.yaml'
    shadow: "shallow"
    params:
        significance=lambda wildcards, input: config['sicer_fdr'] if input.get('control_pileup',
            None) else config['sicer_evalue'],
        signal_pileup_bed_fname=lambda wildcards, input: os.path.basename(input.signal_pileup),
        control_arg=lambda wildcards, input: os.path.basename(input.control_pileup) if input.get('control_pileup', None) else "",
        pileups_dir=lambda wildcards, input: os.path.split(str(input.signal_pileup))[0],
        peaks_file=lambda wildcards, output: os.path.basename(output[0]),
        workdir=WORK_DIR,
        fragment=config['sicer_fragment'],
        genome=config['genome'],
        script=lambda wildcards, input: "SICER.sh" if input.get('control_pileup', None) else "SICER-rb.sh"
    resources:
        threads = 1,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        # SICER.sh ["InputDir"] ["bed file"] ["control file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] [“FDR"]
        #
        # SICER-rb.sh ["InputDir"] ["bed file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] ["E-value"]
        'echo "Significance threshold: {params.significance}" > {params.workdir}/{log} &&'
        ' tmp_sicer=$(mktemp -d) && mkdir -p $tmp_sicer && cd $tmp_sicer && '
        '  {params.script} {params.workdir}/{params.pileups_dir} {params.signal_pileup_bed_fname} {params.control_arg}'
        '    $(pwd) {params.genome} 1 {wildcards.width}'
        '    {params.fragment} $(cat "{params.workdir}/{input.effective_genome_fraction}")'
        '    {wildcards.gap} {params.significance} &>> {params.workdir}/{log} &&'
        ' ls -lah  &>> {params.workdir}/{log} &&'
        ' mv {params.peaks_file} {params.workdir}/{output} &>> {params.workdir}/{log} && '
        ' rm -rf $tmp_sicer'
