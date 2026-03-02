import os

from pipeline_util import *

localrules: download_salmon_fagz

rule salmon_all:
    input: expand('salmon/{sample}/quant.sf', sample=get_rnaseq_samples())


REFERENCE_URLS = {
    "mm39": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz"
}

rule download_salmon_fagz:
    """
    Download transcripts FASTA for a given build (e.g. mm39),
    and decompress them into salmon-index/{build}/...
    """
    output:
        transcripts_fa = f"salmon-index/{config['genome']}/transcripts.fa",
    params:
        url = REFERENCE_URLS[config['genome']]
    log:
        f"logs/salmon-index/download_{config['genome']}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.transcripts_fa})" logs/reference

        curl -L --fail --retry 3 --retry-delay 2 "{params.url}" \
          -o "{output.transcripts_fa}.gz" &>> "{log}"
        gunzip "{output.transcripts_fa}.gz"
        """


rule salmon_index:
    """
    Build a standard (non-decoy) Salmon index
    from transcripts.fa produced by download_salmon_reference_files.
    """
    input:
        transcripts_fa = f"salmon-index/{config['genome']}/transcripts.fa"
    output:
        index = directory(f"salmon-index/{config['genome']}/salmon_index")
    threads: config.get("threads", 8)
    params:
        k = config.get("kmer_size", 31)
    log:
        f"logs/salmon/index_{config['genome']}.log"
    conda: "../envs/salmon.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"

        salmon index \
            -t {input.transcripts_fa} \
            -i {output.index} \
            -k {params.k} \
            -p {threads} \
            &> {log}
        """

rule salmon_quant:
    """
    Quantify paired-end RNA-seq reads using Salmon,
    using the index produced by rule salmon_index.
    """
    input:
        sample = paired_reads_paths(config,True),
        index = rules.salmon_index.output.index
    output:
        quant = "salmon/{sample}/quant.sf"
    threads: config.get("threads", 8)
    params:
        r1=lambda wildcards, input: input.sample[0],
        r2=lambda wildcards, input: input.sample[1],
        libtype = config.get("libtype", "A")
    log:
        "logs/salmon/{sample}.log"
    conda: "../envs/salmon.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"

        salmon quant \
            -i {input.index} \
            -l {params.libtype} \
            -1 {params.r1} \
            -2 {params.r2} \
            -p {threads} \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o "$(dirname {output.quant})" \
            &> {log}
        """