# chipseq-smk-pipeline
Technical ChIP-Seq snakelike pipeline - snakemake version of https://github.com/JetBrains-Research/washu

Configuration 
-------------
The only tool required to launch the pipeline is `conda`.
* If `conda` is not installed,
follow the instructions at
[Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* Navigate to repository directory.

Create a Conda environment for `snakemake`:
```bash
$ conda env create --file envs/snakemake.yaml --name snakemake
```
Activate the newly created environment:
```bash
$ source activate snakemake
```

Launch
------

Run the pipeline:
```bash
$ snakemake all [--cores <cores>] --use-conda --config work_dir=<work_dir> genome=<genome> fastq_dir=<fastq_dir>
```

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
    --config work_dir=<work_dir> fastq_dir=<fastq_dir> genome=<genome> \
    macs2_params="-q 0.05 -f BAMPE --nomodel --nolambda -B --call-summits" \
    span_params="--fragment 0" bin=100
```
