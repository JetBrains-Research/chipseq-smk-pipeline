[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)

# chipseq-smk-pipeline

[Snakemake](https://snakemake.readthedocs.io/en/stable/) based pipeline for ChIP-seq and ATAC-seq datasets processing
from raw data QC and alignment to visualization and peak calling.

Pipeline
--------
Pipeline aligned FASTQ or gzipped FASTQ reads, defined in `config.yaml`.
Reads folder is a relative path in pipeline working directory and defined by `fastq_dir` property.
FASTQ reads extension is defined by `fastq_ext` property, e.g. could be `fq`, `fq.gz`, `fastq`, `fastq.gz`.

| Path          | Description                                                                          |
|---------------|--------------------------------------------------------------------------------------|
| `config.yaml` | Default pipeline options                                                             |
| `trimmed`     | Trimmed FASTQ file, if `trim_reads` option is True.                                  |
| `bams`        | BAMs with aligned reads, `MAPQ >= 30`                                                |
| `bw`          | BAM coverage visualization using [DeepTools](https://doi.org/10.1093/nar/gku365)     |
| `macs2`       | [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137) peaks                              |
| `sicer`       | [SICER](https://doi.org/10.1093/bioinformatics/btp340) peaks                         |
| `span`        | [SPAN](https://doi.org/10.1093/bioinformatics/btab376) peaks                         |
| `qc`          | QC Reports                                                                           |
| `multiqc`     | [MultiQC](https://doi.org/10.1093/bioinformatics/btw354) reports for different steps |
| `logs`        | Shell commands logs                                                                  |

Configuration
-------------
The only tool required to launch the pipeline is `conda`.

* If `conda` is not installed,
  follow the instructions at
  [Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* Navigate to repository directory.

Create a Conda environment for `snakemake`:

```bash
$ conda env create --file environment.yaml --name snakemake
```

Activate the newly created environment:

```bash
$ source activate snakemake
```

Launch
------

Run the pipeline to start with fastq reads:

```bash
$ snakemake all [--cores <cores>] --use-conda --directory <work_dir> fastq_dir=<fastq_dir> --config genome=<genome> 
```

Use `start_with_bams=True` config option to start with existing bam files.

P.S: Use `--config` to override default options from `config.yaml` file

QSUB
----

Configure profile for qsub with Torque scheduler with name `generic_qsub`

```bash
$ mkdir -p ~/.config/snakemake
$ cd ~/.config/snakemake
$ cookiecutter https://github.com/iromeo/generic.git
```

Example of ATAC-Seq processing on qsub

```bash
$ snakemake all --use-conda --profile generic_qsub --cluster-config qsub_config.yaml --jobs 150 \
    --directory <work_dir> \
    --config fastq_dir=<fastq_dir> genome=<genome> \
    macs2=True macs2_params="-q 0.05 -f BAMPE --nomodel --nolambda -B --call-summits" \
    span=True span_params="--fragment 0" bin=100 bowtie2_params="-X 2000 --dovetail"
```

P.S: Use `--config` to override default options from `config.yaml` file

SnakeCharm
----------

ChIP-Seq processing pipeline - [snakemake](https://snakemake.readthedocs.io/en/stable/) version
of [washu](https://github.com/JetBrains-Research/washu) pipeline.\
Developed with [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin
for [PyCharm](https://www.jetbrains.com/pycharm/) IDE.

Useful links
------------

* [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin
* [PyCharm](https://www.jetbrains.com/pycharm/) IDE
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system
* JetBrains Research BioLabs [homepage](https://research.jetbrains.org/groups/biolabs)
