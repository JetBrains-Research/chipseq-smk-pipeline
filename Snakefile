from pipeline_util import *
from pipeline_util import _fastq_paths, _bams_paths, _sample_2_control

# use this file as a basic config file in your working directory
# if you'd like to customise it: fix it directly or override required args
# using --config options or from --configfile file.
configfile: f"{workflow.basedir}/config.yaml"

WORK_DIR = os.getcwd()
FASTQ_DIR = config['fastq_dir']
FASTQ_EXT = config['fastq_ext']
FASTQ_PATHS = _fastq_paths(FASTQ_DIR, FASTQ_EXT)
BAMS_DIR = config['bams_dir']
BAMS_PATHS = _bams_paths(BAMS_DIR)
SAMPLE_2_CONTROL_MAP = _sample_2_control(config, FASTQ_PATHS, BAMS_PATHS)

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

    if bool(config['start_with_bams']):
        print("Bam directory:", BAMS_DIR)
    else:
        print("Fastq directory:", FASTQ_DIR)
        print("Fastq extension:", FASTQ_EXT)

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
include: "rules/reads_bw.smk"

# Peak callers
include: "rules/macs2.smk"
include: "rules/sicer.smk"
include: "rules/sicer2.smk"
include: "rules/omnipeak.smk"
include: "rules/macs3.smk"
include: "rules/homer.smk"
include: "rules/fseq2.smk"
include: "rules/hotspot.smk"
include: "rules/peakseq.smk"
include: "rules/gps.smk"
include: "rules/bayespeak.smk"
include: "rules/lanceotron.smk"

wildcard_constraints:
    sample="[^/]+"

localrules: all

if bool(config['start_with_bams']):
    if not os.path.exists(config['bams_dir']):
        raise ValueError(f"Bam directory not exists: {config['fastq_dir']}")
else:
    if not os.path.exists(config['fastq_dir']):
        raise ValueError(f"Fastq directory not exists: {config['fastq_dir']}")


rule all:
    input:
        *([] if bool(config['start_with_bams']) else [
            # Reads qc
            *([] if not bool(config['reads_qc']) else rules.all_raw_qc_results.input),
            # Optional reads trimming
            *([] if not bool(config['trim_reads']) else rules.all_trim_fastq_results.input),
            # Alignment
            rules.all_alignment_results.input,
            # Alignment qc
            *([] if not bool(config['bams_qc']) else rules.all_alignment_qc.input),
        ]),
        # Visualization
        *([] if not bool(config['bw']) else rules.all_reads_bw_results.input),
        # Visualization tags
        *([] if not bool(config['tagsbw']) else rules.all_tags_bw_results.input),
        # macs2
        *([] if not bool(config['macs2']) else rules.all_macs2_results.input),
        # sicer
        *([] if not bool(config['sicer']) else rules.all_sicer_results.input),
        # sicer2
        *([] if not bool(config['sicer2']) else rules.all_sicer2_results.input),
        # omnipeak
        *([] if not bool(config['omnipeak']) else rules.all_omnipeak_results.input),
        # macs3
        *([] if not bool(config['macs3']) else rules.all_macs3_results.input),
        # homer
        *([] if not bool(config['homer']) else rules.all_homer_results.input),
        # fseq2
        *([] if not bool(config['fseq2']) else rules.all_fseq2_results.input),
        # hotspot
        *([] if not bool(config['hotspot']) else rules.all_hotspot_results.input),
        # peakseq
        *([] if not bool(config['peakseq']) else rules.all_peakseq_results.input),
        # gps
        *([] if not bool(config['gps']) else rules.all_gps_results.input),
        # bayespeak
        *([] if not bool(config['bayespeak']) else rules.all_bayespeak_results.input),
        # lanceotron
        *([] if not bool(config['lanceotron']) else rules.all_lanceotron_results.input)