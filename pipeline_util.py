#!/usr/bin/env python
import os
import re
from glob import glob


def fastq_files(config):
    return glob(os.path.join(config['fastq_dir'], '*.f*q'))


def fastq_names(config):
    return [os.path.splitext(os.path.basename(fastq_file))[0] for fastq_file in fastq_files(config)]


def fastq_common_names_paired(config):
    basenames = {os.path.splitext(fastq_file)[0] for fastq_file in fastq_files(config)}
    paired_candidates = [basename[:-2] for basename in basenames if basename[-2:] == '_1']
    return [os.path.basename(common_name) for common_name in paired_candidates if common_name + '_2' in basenames]


def fastq_names_single(config):
    common_names = fastq_common_names_paired(config)
    return [fastq_name for fastq_name in fastq_names(config)
            if fastq_name[-2:] not in ['_1', '_2'] or fastq_name[:-2] not in common_names]


def fastq_aligned_names(config):
    return fastq_common_names_paired(config) + fastq_names_single(config)


def macs_species(genome):
    """Convert genome to macs2 species encoding"""
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))


def effective_genome_fraction(genome, chrom_sizes_path):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""
    with open(chrom_sizes_path, 'r') as chrom_sizes:
        chrom_length = sum([int(line.split('\t')[1]) for line in chrom_sizes if 'chr_' not in line])
    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return size / chrom_length


def find_control_for(file, ext="bam"):
    bam_name = os.path.basename(file).lower()
    if 'input' in bam_name:
        return ''

    # Find all the files within folder
    parent_dir = os.path.dirname(file)
    names = [os.path.basename(n) for n in glob(f'{parent_dir}/*.{ext}')]
    input_name = _find_control_for_target(bam_name, names)
    if input_name is None:
        return ''
    else:
        return input_name


def _find_control_for_target(target, candidates):
    if _is_control(target):
        return None

    def sort_function(x):
        return _lcs(str(target), x.lower())

    controls = [str(name) for name in candidates if _is_control(name)]
    if len(controls) > 0:
        return max(controls, key=sort_function)
    else:
        return None


def _is_control(c):
    return re.match('.*input.*', re.sub('.*/', '', str(c)), flags=re.IGNORECASE) is not None


def _lcs(x, y):
    """
    Finds longest common subsequence
    Code adopted from https://en.wikibooks.org/wiki/Algorithm_Implementation/
    Strings/Longest_common_subsequence#Python
    """
    m = len(x)
    n = len(y)
    # An (m+1) times (n+1) matrix
    c = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if x[i - 1] == y[j - 1]:
                c[i][j] = c[i - 1][j - 1] + 1
            else:
                c[i][j] = max(c[i][j - 1], c[i - 1][j])

    def back_track(i, j):
        if i == 0 or j == 0:
            return ""
        elif x[i - 1] == y[j - 1]:
            return back_track(i - 1, j - 1) + x[i - 1]
        else:
            if c[i][j - 1] > c[i - 1][j]:
                return back_track(i, j - 1)
            else:
                return back_track(i - 1, j)

    return len(back_track(m, n))


def macs2_input_fun(config):
    def inner(wildcards):
        sample = wildcards.sample

        args = {
            'signal': f'bams/{sample}.bam'
        }

        control_sample = sample_2_control(config)[sample]
        if control_sample:
            args['control'] = f'bams/{control_sample}.bam'

        return args

    return inner


def sample_2_control(config):
    fq_files = fastq_files(config)

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
        return name[-2:]
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


def sicer_all_peaks_input(config):
    files = []

    window_size=config['sicer_window']
    gap=config['sicer_gap']

    # XXX: change significance only via config, SICER rule takes the value from
    # config, not via wildcards

    sample_2_config_dict = sample_2_control(config)
    for sample in fastq_aligned_names(config):
        if sample_2_config_dict[sample] is None:
            # w/o control
            significance=config['sicer_evalue']
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-E{significance}.scoreisland')
        else:
            # with control
            significance=config['sicer_fdr']
            files.append(f'sicer/{sample}-W{window_size}-G{gap}-islands-summary-FDR{significance}')

    return files
