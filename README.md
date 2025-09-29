[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)

# chipseq-smk-pipeline

[Snakemake](https://snakemake.readthedocs.io/en/stable/) based pipeline for ChIP-seq and ATAC-seq datasets processing
from raw data QC and alignment to visualization and peak calling.

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
Supported peak caller tools and parameters:

| Peak caller                                                 | Snakemake                              | Command line                                                                                                                                                                                                                                        |
|-------------------------------------------------------------|----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137)           | [macs2.smk](rules/macs2.smk)           | `macs2 callpeak -f BAM -t <signal.bam> -c <control.bam> -q 0.05`                                                                                                                                                                                    | 
| `MACS2 --broad`                                             | [macs2.smk](rules/macs2.smk)           | `macs2 callpeak -f BAM -t <signal.bam> -c <control.bam> --broad --broad-cutoff 0.1`                                                                                                                                                                 | 
| [MACS3](https://macs3-project.github.io/MACS/)              | [macs3.smk](rules/macs3.smk)           | `macs3 callpeak -f BAM -t <signal.bam> -c <control.bam> -q 0.05`                                                                                                                                                                                    | 
| `MACS3 --broad`                                             | [macs3.smk](rules/macs3.smk)           | `macs3 callpeak -f BAM -t <signal.bam> -c <control.bam> --broad --broad-cutoff 0.1`                                                                                                                                                                 | 
| [SICER](https://doi.org/10.1093/bioinformatics/btp340)      | [sicer.smk](rules/sicer.smk)           | `SICER.sh pileups <signal.pileup.bed> <control.pileup.bed> <outputdir> <species> 1 200 150 <effective genome fraction> 600 0.01`                                                                                                                    | 
| [SICER2](https://github.com/zanglab/SICER2)                 | [sicer2.smk](rules/sicer2.smk)         | `sicer -t <signal.pileup.bed> -c <control.pileup.bed> -s hg38 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000`                                                                                                                               | 
| [HOMER](https://doi.org/10.1016/j.molcel.2010.05.004)       | [homer.smk](rules/homer.smk)           | `findPeaks <signal tags dir> -style factor -o <output> -i <control tags dir>` \*                                                                                                                                                                    | 
| `HOMER histone`                                             | [homer.smk](rules/homer.smk)           | `findPeaks <signal tags dir> -style histone -o <output> -i <control tags dir>` \*                                                                                                                                                                   | 
| [FSeq2](https://doi.org/10.1093/nargab/lqab012)             | [fseq2.smk](rules/fseq2.smk)           | `fseq2 callpeak -v -q_thr 0.05 -control_file <control.pileup.bed> -chrom_size_file <chrom.sizes> -o fseq2 -name <name> -standard_narrowpeak <signal.pileup.bed>`                                                                                    | 
| [HotSpot](https://doi.org/10.1038/ng.759)                   | [hotspot.smk](rules/homer.smk)         | `<hotspot_executable> -i <signal.pileup.bed> -o <signal>.hotspot`                                                                                                                                                                                   | 
| [PeakSeq](https://doi.org/10.1038/nbt.1518)                 | [peakseq.smk](rules/peakseq.smk)       | `<peakseq_executable> -peak_select <config.dat>` **                                                                                                                                                                                                 | 
| [GPS](https://doi.org/10.1093/bioinformatics/btq590)        | [gps.smk](rules/gps.smk)               | `java -cp <gps.jar> edu.mit.csail.cgs.deepseq.discovery.GPS --d <input.signal_reads_distribution> --g <chrom.sizes>                                               --expt <signal.bam> --ctlr <control.bam> --f SAM --out gps/<signal> --q 0.05` *** | 
| [BayesPeak](https://doi.org/10.1186/1471-2105-10-299)       | [bayespeak.smk](rules/bayespeak.smk)   | `Rscript <run_bayespeak.R> <signal.pileup.bed> <control.pileup.bed> <signal>.csv ` ****                                                                                                                                                             | 
| [LanceOtron](https://doi.org/10.1093/bioinformatics/btac525) | [lanceotron.smk](rules/lanceotron.smk) | `lanceotron callPeaksInput <signal.bw> -i <control.bw> -f lanceotron`                                                                                                                                                                               |
| [Omnipeak](https://github.com/JetBrains-Research/omnipeak)  | [omnipeak.smk](rules/omnipeak.smk)     | `java --add-modules=jdk.incubator.vector -Xmx8G -jar <omnipeak.jar> analyze -t <signal.bam> --chrom.sizes <chrom.sizes> -c <control.bam --peaks <output.peak> --iterations 10 --threshold 1e-4 --bin 100 --fragment auto --fdr 0.05`                |
| [SPAN](https://github.com/JetBrains-Research/span)          | [span.smk](rules/span.smk)             | `java -Xmx8G -jar <span.jar> analyze -t <signal.bam> --chrom.sizes <chrom.sizes> -c <control.bam --peaks <output.peak> --iterations 10 --threshold 1e-4 --bin 100 --fragment auto --fdr 0.05`                                                       | 

\* HOMER tags preparation

```
makeTagDirectory <tags dir> <bam>
```

** PeakSeq `config.dat`

```
Experiment_id <signal>
Mappability_map_file <mappability_file>
ChIP_Seq_reads_data_dirs <signal>
Input_reads_data_dirs <control>
chromosome_list_file {chrom.sizes}
narrowPeak_output_file_path peakseq/<signal>.raw
Enrichment_mapped_fragment_length 200
target_FDR 0.05
N_Simulations 50
Minimum_interpeak_distance 200
Background_model Simulated
max_Qvalue 0.05
```

*** GPS Read Distribution

```
java -Xmx1G -cp <gps.jar> edu.mit.csail.cgs.deepseq.analysis.GPS_ReadDistribution --g <chrom.sizes> --coords <input.coords> --chipseq <signal.bam> --f SAM --name gps/<signal> --range 250 --smooth 5 --mrc 4
```

**** BayesPeak [run_bayespeak.R](scripts/run_bayespeak.R)


Peak callers installation
-------------------------
This section contains instructions for manual peak callers installation.

* BayesPeak
    1. Install R
  ```
  mamba install  -c conda-forge r-base=3.6.3
  ```
    2. In R console
  ```
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(version = "3.10")  # Explicitly set correct Bioconductor version
  BiocManager::install(c("IRanges", "GenomicRanges"))
  ```

    3. Install BayesPeak
  ```
  wget https://www.bioconductor.org/packages//2.10/bioc/src/contrib/BayesPeak_1.8.0.tar.gz
  R CMD INSTALL BayesPeak_1.8.0.tar.gz 
  ```

* Hotspot
    1. Install required dependencies
  ```
  sudo apt-get install build-essential libgsl-dev
  ```
    2. Download and make
  ```
  wget https://github.com/StamLab/hotspot/archive/refs/tags/v4.1.1.zip
  gunzip v4.1.1.zip
  cd hotspot-4.1.1/hotspot-distr/hotspot-deploy
  make
  ```
* PeakSeq<br>
  Download and make
  ```
  git clone https://github.com/gersteinlab/PeakSeq.git
  cd PeakSeq
  make
  ```

Rules
-----
Rules DAG produced with additional command line arguments `--forceall --rulegraph | dot -Tpdf > rules.pdf`

![Rules](pipeline.png?raw=true "Rules DAG")

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
    omnipeak=True omnipeak_fragment=0 --rerun-incomplete
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
    --config fastq_ext=fastq.gz fastq_dir=<work_dir> bw=True genome=hg19 macs2=True sicer=True omnipeak=True \
    --rerun-incomplete
```

Useful links
------------

* Learn more about [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system
* Developed with [SnakeCharm](https://plugins.jetbrains.com/plugin/11947-snakecharm) plugin
  for [PyCharm](https://www.jetbrains.com/pycharm/) IDE by JetBrains
  Research [BioLabs](https://research.jetbrains.org/groups/biolabs)
