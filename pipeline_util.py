import os
import re
from glob import glob

from match_control import find_control_for


def _fastq_paths(config):
    return list(glob(os.path.join(config['fastq_dir'], '*.' + config['fastq_ext'])))


def fastq_names_wo_ext(fastq_paths):
    # file name w/o ext and parent folders: supports *.fastq and *.fastq.gz
    return [_split_to_fname_and_ext(f)[0] for f in fastq_paths]


def fastq_aligned_names(fastq_paths):
    return _paired_fastq_samples_names(fastq_paths) + _single_fastq_samples_names(fastq_paths)


def _sample_2_control(fastq_paths):
    result = {}
    for fq_path in fastq_paths:
        ext = _split_to_fname_and_ext(fq_path)[1]
        control_path = find_control_for(fq_path, ext)

        if control_path:
            control_sample = _sample_by_fastq_file(control_path)
        else:
            control_sample = None

        sample = _sample_by_fastq_file(fq_path)

        result[sample] = control_sample

    return result


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


def _sample_by_fastq_file(fq_path):
    name = _split_to_fname_and_ext(fq_path)[0]
    if re.match('.*_R?[12]_?$', name):
        # paired file
        return re.sub('_R?[12]_?$', '', name)
    else:
        # single end file
        return name


def _single_fastq_samples_names(fastq_paths):
    paired_samples = _paired_fastq_samples_names(fastq_paths)
    return [name for name in fastq_names_wo_ext(fastq_paths)
            if not re.match('.*_R?[12]_?$', name) or re.sub('_R?[12]_?$', '', name) not in paired_samples]


def _paired_fastq_samples_names(fastq_paths):
    fq_names = set(fastq_names_wo_ext(fastq_paths))
    result = []
    for name in fq_names:
        if re.match('.*_R?[12]_?$', name):
            sample = re.sub('_R?[12]_?$', '', name)
            suffix = name.replace(sample, '').replace('1', '2')
            if f'{sample}{suffix}' in fq_names:
                result.append(sample)
    return result


def _get_paired_suffixes(config):
    # we assume that all the files share similar naming
    suffix1, suffix2 = '_1', '_2'
    for f in list(glob(os.path.join(config['fastq_dir'], '*.' + config['fastq_ext']))):
        name = _split_to_fname_and_ext(f)[0]
        if re.match('.*_R?[12]_?$', name):
            sample = re.sub('_R?[12]_?$', '', name)
            suffix1 = name.replace(sample, '')
            suffix2 = suffix1.replace('1', '2')
            break
    return suffix1, suffix2


def is_trimmed(config):
    return bool(config['trim_reads'])


def trimmed_fastq_sample_names(fastq_paths):
    # here `name` could have _1 and _2 suffix in case of paired reads
    return [f"{name}_trimmed" for name in fastq_names_wo_ext(fastq_paths)]


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


def tuned_peaks_input_files(config, fastq_paths):
    span_bin = config['span_bin']

    tuned_peaks = []
    for sample in fastq_aligned_names(fastq_paths):
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
            suffix1, suffix2 = _get_paired_suffixes(config)
            return [
                f"trimmed/{{sample}}{suffix1}_trimmed.{config['fastq_ext']}",
                f"trimmed/{{sample}}{suffix2}_trimmed.{config['fastq_ext']}"
            ]
        else:
            return [
                f"trimmed/{{sample}}_trimmed.{config['fastq_ext']}",
            ]
    else:
        if paired:
            suffix1, suffix2 = _get_paired_suffixes(config)
            return [
                config['fastq_dir'] + f"/{{sample}}{suffix1}.{config['fastq_ext']}",
                config['fastq_dir'] + f"/{{sample}}{suffix2}.{config['fastq_ext']}",
            ]
        else:
            return [
                config['fastq_dir'] + f"/{{sample}}.{config['fastq_ext']}"
            ]
