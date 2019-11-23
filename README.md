# chipseq-smk-pipeline
ChIP-Seq processing pipeline - [snakemake](https://snakemake.readthedocs.io/en/stable/) version of [washu](https://github.com/JetBrains-Research/washu) pipeline.\
Developed with [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin for [PyCharm](https://www.jetbrains.com/pycharm/) IDE.

Configuration 
-------------
The only tool required to launch the pipeline is `conda`.
* If `conda` is not installed,
follow the instructions at
[Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* Navigate to repository directory.

Create a Conda environment for `snakemake`:
```bash
$ conda env create --file ./environment.yaml --name snakemake
```
Activate the newly created environment:
```bash
$ source activate snakemake
```

Launch
------

Run the pipeline:
```bash
$ snakemake all [--cores <cores>] --use-conda --directory <work_dir> --config pipeline_src_path=<pipeline_src_dir> genome=<genome> fastq_dir=<fastq_dir>
```

P.S: Use `--config` to override default options from `./config.yaml` file

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
    --config pipeline_src_path=<pipeline_src_dir> fastq_dir=<fastq_dir> genome=<genome> \
    macs2_params="-q 0.05 -f BAMPE --nomodel --nolambda -B --call-summits" \
    span_params="--fragment 0" bin=100 bowtie2_params="-X 2000 --dovetail"
```

P.S: Use `--config` to override default options from `./config.yaml` file

Pipeline Description
-----
Pipeline aligned FASTQ or gzipped FASTQ reads, defined in `config.yaml`. 
Reads folder is a relative path in pipeline working directory and defined by `fastq_dir` property. FASTQ reads extension is defined by `fastq_ext` property, e.g. could be `fq`, `fq.gz`, `fastq`, `fastq.gz`.


| Folder | Description |
| --- | --- |
| `./config.yaml` | Default pipeline options |
| `./trimmed` | Trimmed FASTQ file, if `trim_reads` option is True. BAM files will contain trimmed reads in this case |
| `./bams` | BAMs with aligned reads, `MAPQ >= 30` |
| `./deduplicated` | Deduplicated version of `./bams` files |
| `./macs2` | MACS2 peaks |
| `./sicer` | SICER peaks |
| `./span` | SPAN peaks |
| `./multiqc` | MultiQC reports for different steps |
| `./qc` | QC Reports |
| `./logs` | Shell commands logs |

Useful links
------------
* [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin
* [PyCharm](https://www.jetbrains.com/pycharm/) IDE
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system
* JetBrains Research BioLabs [homepage](https://research.jetbrains.org/groups/biolabs)
