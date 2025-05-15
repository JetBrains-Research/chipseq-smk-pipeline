[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)

# chipseq-smk-pipeline

[Snakemake](https://snakemake.readthedocs.io/en/stable/) based pipeline for ChIP-seq and ATAC-seq datasets processing
from raw data QC and alignment to visualization and peak calling.

![Scheme](pipeline.png?raw=true "Pipeline")

*During peak calling steps `chipseq-smk-pipeline` automatically matches signal with control file by names proximity.*


Input
-----
**Input FASTQ files**

Pipeline aligned FASTQ or gzipped FASTQ reads, defined in `config.yaml`.<br>
Reads folder is a relative path in pipeline working directory and defined by `fastq_dir` property.<br>
FASTQ reads extension is defined by `fastq_ext` property, e.g. could be `fq`, `fq.gz`, `fastq`, `fastq.gz`.

**Input BAM files**

Use `start_with_bams=True` config option to start with existing bam files.<br>
Pipeline starts with `BAM` files in `work_dir/bams` folder.


Files
-----

| Path                 | Description                                                                          |
|----------------------|--------------------------------------------------------------------------------------|
| `config.yaml`        | Default pipeline options                                                             |
| `trimmed`            | Trimmed FASTQ file, if `trim_reads` option is True.                                  |
| `bams`               | BAMs with aligned reads, `MAPQ >= 30`                                                |
| `bw`                 | BAM coverage visualization using [DeepTools](https://doi.org/10.1093/nar/gku365)     |
| `<peak_caller_name>` | Peaks provided by peak caller tool `<peak_caller_name>`                              |
| `qc`                 | QC Reports                                                                           |
| `multiqc`            | [MultiQC](https://doi.org/10.1093/bioinformatics/btw354) reports for different steps |
| `logs`               | Shell commands logs                                                                  |

Requirements
------------
The pipeline requires `conda`.

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

On Ubuntu please ensure that `gawk` is installed:

```bash
$ sudo apt-get install gawk
```

Launch
------

Run the pipeline to start with fastq reads:

```bash
$ snakemake -p -s <chipseq-smk-pipeline>/Snakefile \
    all [--cores <cores>] --use-conda --directory <work_dir> \
    --config fastq_dir=<fastq_dir> genome=<genome> --rerun-incomplete
```

The Default pipeline doesn't perform coverage visualization and launch peak callers.<br>
Please add `bw=True`, `<peak_caller_name>=True` to create coverage bw files and call peaks with `<peak_caller_name>`.

See `config.yaml` for a complete list of parameters. Use`--config` to override default options from `config.yaml` file.

Peak callers
------------
Supported peak caller tools:

* [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137)
* [MACS3](https://macs3-project.github.io/MACS/)
* [SICER](https://doi.org/10.1093/bioinformatics/btp340)
* [SPAN](https://doi.org/10.1093/bioinformatics/btab376)
* [HOMER](https://doi.org/10.1016/j.molcel.2010.05.004)
* [FSeq2](https://doi.org/10.1093/nargab/lqab012)
* [HotSpot](https://doi.org/10.1038/ng.759)
* [PeakSeq](https://doi.org/10.1038/nbt.1518)

To launch MACS2 in `--broad` mode, use the following config:

```bash
$ snakemake -p -s <chipseq-smk-pipeline>/Snakefile \
    all [--cores <cores>] --use-conda --directory <work_dir> \
    --config fastq_dir=<fastq_dir> genome=<genome> \
    macs2=True macs2_mode=broad macs2_params="--broad --broad-cutoff 0.1" macs2_suffix=broad0.1 \
    --rerun-incomplete
```

Rules
-----
Rules DAG produced with additional command line agruments `--forceall --rulegraph | dot -Tpdf > rules.pdf`

![Rules](rules.png?raw=true "Rules DAG")

Computational cluster QSUB/LFS/QSUB
-----------------------------------

Configure profile for required cluster system with name `cluster`.

```bash
$ mkdir -p ~/.config/snakemake
$ cd ~/.config/snakemake
$ cookiecutter https://github.com/iromeo/generic.git
```

Example of ATAC-Seq processing on qsub

```bash
$ snakemake -p -s <chipseq-smk-pipeline>/Snakefile \
    all --use-conda --directory <work_dir> \
    --profile cluster --cluster-config cluster_config.yaml --jobs 150 \
    --config fastq_dir=<fastq_dir> genome=<genome> \
    bowtie2_params="-X 2000 --dovetail" \
    macs2=True macs2_params="-q 0.05 -f BAMPE --nomodel --nolambda -B --call-summits" \
    span=True span_fragment=0 --rerun-incomplete
```

P.S: Use `--config` to override default options from `config.yaml` file

Try with test data
------------------

Please download example `fastq.gz` files
from [CD14_chr15_fastq](https://artyomovlab.wustl.edu/publications/supp_materials/4Oleg/CD14_chr15_fastq/) folder.<br>
These files are filtered on human hg19 chr15 to reduce size and make computations faster.

Launch `chipseq-smk-pipeline`:

```bash
$ snakemake -p -s <chipseq-smk-pipeline>/Snakefile \
    all --use-conda --cores all --directory <work_dir> \
    --config fastq_ext=fastq.gz fastq_dir=<work_dir> bw=True genome=hg19 macs2=True sicer=True span=True \
    --rerun-incomplete
```

Useful links
------------

* Learn more about [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system
* Developed with [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin
  for [PyCharm](https://www.jetbrains.com/pycharm/) IDE by JetBrains
  Research [BioLabs](https://research.jetbrains.org/groups/biolabs)
