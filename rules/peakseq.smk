from pipeline_util import *

######## Step: Peak Calling: Peakseq ##################
rule all_peakseq_results:
    input:
        peakseq_peaks=expand(f'peakseq/{{sample}}.narrowPeak',
            sample=filter(lambda f: not is_control(f) and f in SAMPLE_2_CONTROL_MAP,
                aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


rule download_peakseq_mappability:
    output: f'peakseq/{config["genome"]}.mapp'
    params:
        genome=config["genome"]
    shell:
        'wget -O {output}.raw \
        http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/Mappability_Maps/{params.genome}_100bp.mapp && '
        'cat {output}.raw | grep -v \'_\' > {output} && rm {output}.raw'

rule bam_to_sam:
    input: f"{config['bams_dir']}/{{sample}}.bam"
    output: temp(f"peakseq/{{sample}}.sam")
    conda: '../envs/bio.env.yaml'
    shell:
        'samtools view -h {input} | awk \'/^@/ || $3 ~ /^chr([0-9]+|X|Y|M)$/\' > {output}'


rule bam_preprocess:
    input: rules.bam_to_sam.output
    output: temp(directory("peakseq/{sample}_reads"))
    params:
        peakseq_executable=config['peakseq_executable']
    shell:
        'mkdir -p {output} && '
        '{params.peakseq_executable} -preprocess SAM {input} {output} && '
        'mv {output}/chr_ids.txt {output}/chr_ids.txt.raw &&'
        'cat {output}/chr_ids.txt.raw | grep -v \'_\' > {output}/chr_ids.txt && '
        'rm {output}/chr_ids.txt.raw'


def peakseq_config_dat_input_fun(wildcards):
    sample = wildcards.sample
    control_sample = SAMPLE_2_CONTROL_MAP[sample]

    return dict(
        mappability_file=rules.download_peakseq_mappability.output,
        signal=f"peakseq/{sample}_reads",
        control=f"peakseq/{control_sample}_reads",
        chrom_sizes=rules.download_chrom_sizes.output,
    )

def peakseq_input_fun(wildcards):
    sample = wildcards.sample
    pcdif = peakseq_config_dat_input_fun(wildcards)
    return dict(
        config_dat=f"peakseq/{sample}.config.dat",
        **pcdif
    )


rule prepare_peakseq_config:
    input: unpack(peakseq_config_dat_input_fun)
    output: 'peakseq/{sample}.config.dat'
    shell:
        'echo "" > {output} && '
        'echo "Experiment_id {wildcards.sample}" >> {output} && '
        'echo "Mappability_map_file {input.mappability_file}" >> {output} && '
        'echo "ChIP_Seq_reads_data_dirs {input.signal}" >> {output} && '
        'echo "Input_reads_data_dirs {input.control}" >> {output} &&'
        'echo "chromosome_list_file {input.chrom_sizes}" >> {output} && '
        'echo "narrowPeak_output_file_path peakseq/{wildcards.sample}.raw" >> {output} &&'
        'echo "Enrichment_mapped_fragment_length 200" >> {output} && '
        'echo "target_FDR 0.05" >> {output} && '
        'echo "N_Simulations 50" >> {output} && '
        'echo "Minimum_interpeak_distance 200" >> {output} && '
        'echo "Background_model Simulated" >> {output} && '
        'echo "max_Qvalue 0.05" >> {output}'


rule call_peaks_peakseq:
    input: unpack(peakseq_input_fun)
    output: 'peakseq/{sample}.narrowPeak'
    log: f'logs/peakseq/{{sample}}.log'
    params:
        peakseq_executable=config['peakseq_executable']
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
        '{params.peakseq_executable} -peak_select {input.config_dat} &> {log} && '
        'cat peakseq/{wildcards.sample}.raw | awk \'{{printf("chr%s\\n", $0)}}\' |\
         sort -k1,1 -k2,2n -k3,3n > {output} && '
        'rm peakseq/{wildcards.sample}.raw'