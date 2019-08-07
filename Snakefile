from pipeline_util import *

configfile: "config.yaml"
workdir: config['work_dir']


include: "rules/step1_raw_qc.smk"
include: "rules/step2_trim_fastq.smk"
include: "rules/step3_alignment.smk"
include: "rules/step4_filter_aligned_reads.smk"
include: "rules/step5_reads_coverage.smk"
include: "rules/step6_bam_quality_metrics.smk"
include: "rules/step7_macs2.smk"
include: "rules/step8_sicer.smk"
include: "rules/step9_span.smk"

wildcard_constraints:
    sample="[^/]+"

localrules: all

if not os.path.exists(config['fastq_dir']):
    raise ValueError(f"Reads directory not exists: {config['fastq_dir']}")

rule all:
    input:
        # Reads qc
        rules.step1_raw_qc_results.input,

        # Optional reads trimming, this option is controlled by setting: config[trim_reads]
        *([] if not is_trimmed(config) else rules.step2_trim_fastq_results.input),

        # Alignment
        rules.step3_alignment_results.input,

        # Filter only aligned reads: not all peak callers capable to exclude
        # unaligned reads or reads aligned with bad quality
        # Optionally deduplicated bams and save to 'deduplicated' folder
        rules.step4_filter_aligned_reads_results.input,

        # Visualization
        rules.step5_reads_coverage_results.input,

        # Optional: Quality metrics
        rules.step6_bam_quality_metrics_results.input,

        # macs2
        rules.step7_macs2_results.input,

        # sicer
        rules.step8_sicer_results.input,

        # span
        rules.step9_span.input,
        rules.step9_span_tuned.input
