from pipeline_util import *

localrules: all_bam_quality_metrics_results

######## Step: BAMS QC Metrics ##################
rule all_bam_quality_metrics_results:
    input:
        bam_qc_phantom=expand('qc/phantom/{sample}.phantom.tsv', sample=fastq_aligned_names(FASTQ_PATHS)),
        bam_qc_pbc=expand('qc/pbc_nrf/{sample}.pbc_nrf.tsv', sample=fastq_aligned_names(FASTQ_PATHS)),


rule download_phantompeakqualtools:
    output: directory('bin/phantompeakqualtools')
    log: 'logs/phantompeakqualtools.log'

    conda: "../envs/phantom.env.yaml"
    params:
          targz='phantompeakqualtools.tar.gz'
    shell:
        'cd bin &> {log} && '
        'curl --location '
        'https://storage.googleapis.com/google-code-archive-downloads/v2/'
        'code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz '
        '--output {params.targz} &>> ../{log} && '
        'tar xvf {params.targz} &>> ../{log}'

rule install_spp:
    output: touch('flags/spp.installed')
    log: 'logs/spp/installation.log'

    conda: "../envs/phantom.env.yaml"
    script: "../scripts/spp_install.R"


# This rule requires spp R package installed
rule bam_qc_phantom:
    input:
         ppqt_dir=rules.download_phantompeakqualtools.output,
         spp_flag=rules.install_spp.output,
         bam='bams/{sample}.bam'
    output:
        tsv='qc/phantom/{sample}.phantom.tsv',
        pdf='qc/phantom/{sample}.pdf',
    log: 'logs/phantom/{sample}.phantom.tsv'

    conda: "../envs/phantom.env.yaml"
    params:
          run_spp=lambda wildcards, input: os.path.join(str(input.ppqt_dir), 'run_spp.R')

    shell:
        # > awk: line 2: function and never defined
        # on Ubuntu also gawk required: sudo apt-get install gawk
        # 'Rscript {params.run_spp} -c={input.bam} -savp -out={output.tsv} -rf'
        'Rscript --default-packages=utils,stats,grDevices,caTools,graphics {params.run_spp}'
        ' -c={input.bam} -savp -out={output.tsv} -odir=qc/phantom -rf'

rule bam_qc_pbc_nrf:
    input: 'bams/pileup/{sample}.bed'
    output: 'qc/pbc_nrf/{sample}.pbc_nrf.tsv'
    params:
          tmp_dir='tmp'
    shell: '''
mkdir -p {params.tmp_dir} &&
(T=$'\\t'
>&2 echo "TotalReadPairs${{T}}DistinctReadPairs${{T}}OneReadPair${{T}}TwoReadPairs${{T}}\
NRF=Distinct/Total${{T}}PBC1=OnePair/Distinct${{T}}PBC2=OnePair/TwoPair"

cat {input} | \
    sort -k1,1 -k3,3n -k2,2n -k6,6 -T {params.tmp_dir} | \
    awk -v OFS='\\t' '{{print $1,$2,$3,$6}}' | uniq -c | \
    awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}}
    ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}}
    END{{
        if (mt!=0){{m0_t=m0/mt}} else {{m0_t=-1.0}};
        if (m0!=0){{m1_0=m1/m0}} else {{m1_0=-1.0}};
        if (m2!=0){{m1_2=m1/m2}} else {{m1_2=-1.0}};
        printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0_t,m1_0,m1_2;
    }}') > {output}
    '''