from pipeline_util import *

######## Step: Peak Calling: LanceOtron ##################
rule all_lanceotron_results:
    input:
        lanceotron_peaks=expand(f'lanceotron/{{sample}}_L-tron.bed',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


rule bam_to_rpkm_bw:
    input:
        bam=f"{config['bams_dir']}/{{sample}}.bam",
        index=f"{config['bams_dir']}/{{sample}}.bam.bai"
    output: temp(f"lanceotron/{{sample}}.bw")
    log: 'logs/lanceotron/bw/{sample}.log'
    conda: '../envs/deeptools.env.yaml'
    threads: 8
    params:
        lanceotron_bw_params=config['lanceotron_bw_params']
    resources:
        threads = 8,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output} {params.lanceotron_bw_params} &> {log}'


rule call_peaks_lanceotron:
    input: 'lanceotron/{sample}.bw'
    output: 'lanceotron/{sample}_L-tron.bed'
    log: f'logs/lanceotron/{{sample}}.log'
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
        # LanceOtron cuts extension of the bw file before first . to compute resulting name
        't=$(mktemp -u -p lanceotron -t XXXXX) && nt=$(basename $t) && cp {input} $t && '
        'lanceotron callPeaks $t -f lanceotron &> {log} && rm $t && '
        'mv lanceotron/${{nt}}_L-tron.bed lanceotron/{wildcards.sample}_L-tron.tsv && '
        'cat lanceotron/{wildcards.sample}_L-tron.tsv | '
        'tail -n +2 | awk -v OFS="\\t" \'{{print($1, $2, $3, ".", $4)}}\' > {output}'