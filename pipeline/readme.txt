Current usage:

$ conda env create --name pipeline-chipseq --file environment.yaml
$ conda env create --name macs2 --file environment2.yaml

$ conda activate pipeline-chipseq
$ snakemake --config work_dir=<work dir> genome=<genome build> fastq_dir=<fastq dir>