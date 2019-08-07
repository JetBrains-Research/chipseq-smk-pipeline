import os
import re
from glob import glob

from match_control import find_control_for

def fastq_paths(config):
    fq_dir = config['fastq_dir']
    fq_ext = config['fastq_ext']
    return list(glob(os.path.join(fq_dir, '*.' + fq_ext)))


def fastq_names_wo_ext(config):
    # file name w/o ext and parent folders: supports *.fastq and *.fastq.gz
    return [_split_to_fname_and_ext(f)[0] for f in fastq_paths(config)]


def fastq_aligned_names(config):
    return _paired_fastq_samples_names(config) + _single_fastq_samples_names(config)


def sample_2_control(config):
    fq_files = fastq_paths(config)

    result = {}
    for fq_path in fq_files:
        ext = _split_to_fname_and_ext(fq_path)[1]
        control_path = find_control_for(fq_path, ext)

        if control_path:
            control_sample = _sample_by_fastq_file(control_path)
        else:
            control_sample = None

        sample = _sample_by_fastq_file(fq_path)

        result[sample] = control_sample

    return result


def _sample_by_fastq_file(fq_path):
    name = _split_to_fname_and_ext(fq_path)[0]

    if name[-2:] in ['_1', '_2']:
        # paired file
        return name[:-2]
    else:
        # single end file
        return name


def _split_to_fname_and_ext(path):
    # assumes file has valid ext
    # understands *.{ext} and *.{ext}.gz

    fname = os.path.basename(path)

    name, dot_ext = os.path.splitext(fname)
    if dot_ext == ".gz":
        name2, dot_ext2 = os.path.splitext(name)

        # remove first dot from ext:
        return name2, dot_ext2[1:] + dot_ext
    else:
        # remove first dot from ext:
        return name, dot_ext[1:]


def _single_fastq_samples_names(config):
    paired_samples = _paired_fastq_samples_names(config)
    return [name for name in fastq_names_wo_ext(config)
            if name[-2:] not in ['_1', '_2'] or name[:-2] not in paired_samples]


def _paired_fastq_samples_names(config):
    fq_names = set(fastq_names_wo_ext(config))
    paired_samples = [name[:-2] for name in fq_names if name[-2:] == '_1']
    return [sample for sample in paired_samples if sample + '_2' in fq_names]


def is_trimmed(config):
    return bool(config['trim_reads'])


def trimmed_fastq_sample_names(config):
    # here `name` could have _1 and _2 suffix in case of paired reads
    return [f"{name}_trimmed" for name in fastq_names_wo_ext(config)]


def effective_genome_fraction(genome, chrom_sizes_path, pileup_bed):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""

    # Get chr names covered with reads (e.g. if data is filtered by chromosome name
    # or some chrs excluded during alignment
    chromosomes = set()
    with open(str(pileup_bed)) as f:
        for line in f:
            chr = line.split()[0]
            chromosomes.add(chr)

    # Sized of chromosomes covered with reads
    chrom_sizes = {}
    with open(str(chrom_sizes_path)) as f:
        for line in f:
            chromosome, size = line.split()
            chrom_sizes[chromosome] = int(size)

    # Normalization if not all genome chromosomes are covered
    chromosomes_length = sum([chrom_sizes.get(c, 0) for c in chromosomes])
    genome_length = sum(chrom_sizes.values())

    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return (size / genome_length) * (1.0 * chromosomes_length / genome_length)


def macs_species(genome):
    """Convert genome to macs2 species encoding"""
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))


def tuned_peaks_input_files(config):
    span_bin = config['span_bin']

    tuned_peaks = []
    for sample in fastq_aligned_names(config):
        labels_file = find_labels_for_sample(sample, config)
        if labels_file:
            tuned_peaks.append(f'span/{sample}_{span_bin}_tuned.peak')

    return tuned_peaks


def find_labels_for_sample(sample, config):
    span_markup = config['span_markup']
    if span_markup:
        return span_markup

    labels2files_dict = labels2files(config)
    chunks = sample.split("_")
    for ch in chunks:
        labels_file = labels2files_dict.get(ch, None)
        if labels_file:
            return labels_file
    return None


def labels2files(config):
    labels2files_dict = {}
    if os.path.exists(config['span_labels_dir']):
        for f in glob(os.path.join(config['span_labels_dir'], '*_labels.bed')):
            fname = os.path.basename(f)
            label_name = fname.replace('_labels.bed', '')
            labels2files_dict[label_name] = f

    return labels2files_dict


def bowtie2_input_paths(config, paired):
    if is_trimmed(config):
        if paired:
            return [
                "trimmed/{sample}_1_trimmed.fq.gz",
                "trimmed/{sample}_2_trimmed.fq.gz"
            ]
        else:
            return [
                "trimmed/{sample}_trimmed.fq.gz",
            ]
    else:
        if paired:
            return [
                config['fastq_dir'] + f"/{{sample}}_1.{config['fastq_ext']}",
                config['fastq_dir'] + f"/{{sample}}_2.{config['fastq_ext']}",
            ]
        else:
            return [
                config['fastq_dir'] + f"/{{sample}}.{config['fastq_ext']}"
            ]
