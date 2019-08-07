from pipeline_util import *

localrules: all_sicer_results

# !!!!!
# SICER doesn't support out of the box hg38 and mm10, needs to be tweaked a bit
# !!!!!

######## Step: Peak Calling: SICER ##################
def sicer_all_peaks_input():
    files = []

    window_size=config['sicer_window']
    gap=config['sicer_gap']

    # XXX: change significance only via config, SICER rule takes the value from
    # config, not via wildcards

    sample_2_config_dict = sample_2_control(config)
    for sample in fastq_aligned_names(config):
        if sample_2_config_dict[sample] is None:
            # w/o control
            significance=config['sicer_evalue']
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-E{significance}.scoreisland')
        else:
            # with control
            significance=config['sicer_fdr']
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-islands-summary-FDR{significance}')

    return files

rule all_sicer_results:
    input: *sicer_all_peaks_input()


rule bam_to_pileup:
    input: 'bams/{sample}.bam'
    output: temp('bams/pileup/{sample}.bed')

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
    control_sample = sample_2_control(config)[sample]
    if control_sample:
        control_args['control_pileup'] = f'bams/pileup/{control_sample}.bed'

    return dict(
        # pileup_bed='bams/pileup/{sample}.bed',
        signal_pileup = f'bams/pileup/{sample}.bed',
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
        fragment=config['sicer_fragment'],
        genome=config['genome'],
        script=lambda wildcards, input: "SICER.sh" if input.get('control_pileup', None) else "SICER-rb.sh"
    resources:
        threads = 4,
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
        'echo "Significance threshold: {params.significance}" > {log} &&'
        ' mkdir -p tmp_sicer &&'
        ' cd tmp_sicer && '
        '  {params.script} ../{params.pileups_dir} {params.signal_pileup_bed_fname} {params.control_arg}'
        '    $(pwd) {params.genome} 1 {wildcards.width}'
        '    {params.fragment} $(cat "../{input.effective_genome_fraction}")'
        '    {wildcards.gap} {params.significance} &>> ../{log} &&'
        ' ls -lah  &>> ../{log} &&'
        ' mv {params.peaks_file} ../{output} &>> ../{log}'
