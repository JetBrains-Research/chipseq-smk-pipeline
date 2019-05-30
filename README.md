# chipseq-smk-pipeline
Technical ChIP-Seq snakelike pipeline - snakemake version of https://github.com/JetBrains-Research/washu

Launching Pipeline
------------------
The only tool required to launch the pipeline is `conda`. You'll also need `bowtie2`
indexes for the appropriate genome build.
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
Run the pipeline:
```bash
$ snakemake all [--cores <cores>] --use-conda --config work_dir=<work dir> genome=<genome build> fastq_dir=<fastq dir>
```
