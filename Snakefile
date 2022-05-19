from pipeline_util import *
from pipeline_util import _fastq_paths, _sample_2_control

# use this file as a basic config file in your working directory
# if you'd like to customise it: fix it directly or override required args
# using --config options or from --configfile file.
configfile: f"{workflow.basedir}/config.yaml"

WORK_DIR = os.getcwd()
FASTQ_PATHS = _fastq_paths(config)
SAMPLE_2_CONTROL_MAP = _sample_2_control(config, FASTQ_PATHS)

onstart:
    print(f"Working directory: {WORK_DIR}")
    print(f"Snakefile directory: {workflow.basedir}")
    print(f"Environment: TMPDIR={os.environ.get('TMPDIR', '<n/a>')}")
    print(f"Environment: PATH={os.environ.get('PATH', '<n/a>')} ")
    print(f"Config: ", *[f'{k}: {v}' for k, v in config.items()], sep = "\n  ")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')
    os.system('')

    # check shell (cond not work properly due to shell detection issues
    print("Snakemake shell check")
    shell('echo "  SNAKEMAKE VERSION: $(snakemake --version)"')

    print("FastQ Reads:", config['fastq_dir'])

    #---------------------------------------------------------------------
    # Let's create symlinks for several pipeline source dirs to simplify
    # further paths in pipeline
    for pipeline_dir in ['scripts', 'envs', 'schemas']:
        if not os.path.exists(pipeline_dir):
            src_dir = os.path.join(workflow.basedir , pipeline_dir)
            if os.path.exists(src_dir):
                print(f"Linking '{pipeline_dir}' directory:")
                shell(f"ln -sf {src_dir} {pipeline_dir}")

include: "rules/raw_qc.smk"
include: "rules/trim_fastq.smk"
include: "rules/alignment.smk"
include: "rules/deduplicated_reads.smk"
include: "rules/reads_coverage.smk"
include: "rules/bam_quality_metrics.smk"
include: "rules/macs2.smk"
include: "rules/sicer.smk"
include: "rules/span.smk"

wildcard_constraints:
    sample="[^/]+"

localrules: all

if not os.path.exists(config['fastq_dir']):
    raise ValueError(f"Reads directory not exists: {config['fastq_dir']}")

rule all:
    input:
        # Reads qc
        rules.all_raw_qc_results.input,

        # Optional reads trimming, this option is controlled by setting: config[trim_reads]
        *([] if not bool(config['trim_reads']) else rules.all_trim_fastq_results.input),

        # Alignment
        rules.all_alignment_results.input,

        # Deduplicated bams saved to 'deduplicated' folder
        rules.all_deduplicated_reads_results.input,

        # Visualization
        rules.all_reads_coverage_results.input,

        # Optional: Quality metrics
        rules.all_bam_quality_metrics_results.input,

        # macs2
        rules.all_macs2_results.input,

        # sicer
        rules.all_sicer_results.input,

        # span
        rules.all_span.input
