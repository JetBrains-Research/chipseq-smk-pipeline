from pipeline_util import *


######## Step: Peak Calling: Fseq2 ##################
rule all_fseq2_results:
    input:
        fseq2_peaks=expand(f'fseq2/{{sample}}_peaks.narrowPeak',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )

rule filter_pileup_bed:
    input: f"{config['bams_dir']}/pileup/{{sample}}.bed"
    output: temp(f"{config['bams_dir']}/pileup/{{sample}}.filtered.bed")
    shell: "cat {input} | grep -v \'_\' > {output}"

def fseq2_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control_pileup'] = f"{config['bams_dir']}/pileup/{control_sample}.filtered.bed"

    return dict(
        signal_pileup = f"{config['bams_dir']}/pileup/{{sample}}.filtered.bed",
        **control_args,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_fseq2:
    input: unpack(fseq2_input_fun)
    output:
        peaks=f'fseq2/{{sample}}_peaks.narrowPeak'
    log: f'logs/fseq2/{{sample}}.log'
    conda: '../envs/fseq2.env.yaml'
    params:
        control_arg=lambda wildcards, input: \
            f" -control_file {input.control_pileup}" if input.get('control_pileup', None) else "",
        # Fseq2 is not well optimized by memory usage, artificially limit CPU for the sake of memory
        threads=config['fseq2_effective_threads']
    resources:
        mem = 48, mem_ram = 48,
        time = 60 * 120
    threads: config['fseq2_threads']
    shell:
        'mkdir -p fseq2 &&'
        'fseq2 callpeak -v -cpus {params.threads} -q_thr 0.05 {params.control_arg} -chrom_size_file {input.chrom_sizes} '
        '-o fseq2 -name {wildcards.sample} -standard_narrowpeak {input.signal_pileup} &> {log}'
