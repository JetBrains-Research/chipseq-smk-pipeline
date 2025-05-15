from pipeline_util import *

localrules: download_gps

######## Step: Peak Calling: GPS(GEM) ##################
rule all_gps_results:
    input:
        gps_peaks=expand(f'gps/{{sample}}.bed',
            sample=filter(lambda f: not is_control(f), aligned_names(config, FASTQ_PATHS, BAMS_PATHS))
        )


rule download_gps:
    output: 'bin/gem.jar'
    params:
        workdir=WORK_DIR
    shell:
        'mkdir -p {params.workdir}/bin && tmp_gem=$(mktemp -d) && mkdir -p $tmp_gem && cd $tmp_gem && '
        'wget https://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz && '
        'tar xvf gem.v3.4.tar.gz && mv gem/gem.jar {params.workdir}/{output} && cd .. && rm -rf $tmp_gem'

rule coordinates:
    input:
        signal=f"{config['bams_dir']}/{{sample}}.bam",
        macs2_peaks=f'macs2/{{sample}}_{config["macs2_suffix"]}_peaks.{config["macs2_mode"]}Peak',
        chromosomes_sizes=rules.download_chrom_sizes.output,
    output: temp("gps/{sample}.coords.txt")
    shell:
        'cat {input.macs2_peaks}  | awk -v OFS=\':\' \'{{print $1,int(($2 + $3) / 2)}}\' > {output}'


rule reads_distribution:
    input:
        signal=f"{config['bams_dir']}/{{sample}}.bam",
        coordinates=rules.coordinates.output,
        chromosomes_sizes=rules.download_chrom_sizes.output,
        gps=rules.download_gps.output
    output: temp("Read_Distribution_gps/{sample}.txt")
    shell:
        'mkdir -p Read_Distribution_gps/ && '
        'java -Xmx1G -cp {input.gps} edu.mit.csail.cgs.deepseq.analysis.GPS_ReadDistribution '
        '--g {input.chromosomes_sizes} --coords {input.coordinates} --chipseq {input.signal} --f SAM '
        '--name gps/{wildcards.sample} --range 250 --smooth 5 --mrc 4'


def gps_input_fun(wildcards):
    sample = wildcards.sample

    control_args = {}
    if sample in SAMPLE_2_CONTROL_MAP:
        control_sample = SAMPLE_2_CONTROL_MAP[sample]
        control_args['control'] = f"{config['bams_dir']}/{control_sample}.bam"

    return dict(
        signal=f"{config['bams_dir']}/{sample}.bam",
        signal_reads_distribution=rules.reads_distribution.output,
        **control_args,
        gps=rules.download_gps.output,
        chrom_sizes=rules.download_chrom_sizes.output,
    )


rule call_peaks_gps:
    input: unpack(gps_input_fun)
    output:
        peaks='gps/{sample}.bed'
    log: 'logs/gps/{sample}.log'
    conda: '../envs/java.env.yaml'
    params:
        control_arg=lambda wildcards, input: f" --ctrl {input.control}" if input.get('control', None) else "",
    resources:
        mem = 12, mem_ram = 8,
        time = 60 * 120
    shell:
         'java -cp {input.gps} edu.mit.csail.cgs.deepseq.discovery.GPS '
         '--d {input.signal_reads_distribution} --g {input.chrom_sizes} '
         '--expt {input.signal} {params.control_arg} --f SAM --out gps/{wildcards.sample} --q 0.05 &> {log} && '
         'cat gps/{wildcards.sample}_2_GPS_significant.txt | sed \'s#:#\t#g\' | tail -n +2 | '
         'awk \'{{printf("chr%d\t\t%d\\tt%d\\tn", $1+1, $2-1000, $2+1000)}}\' | '
         'sort -k1,1 -k2,2n -k3,3n > {output}'
