import os
import re
from glob import glob


####################
# Fastq processing #
####################

def is_control(c):
    return re.match('.*(input|control|wce).*', re.sub('.*/', '', str(c)), flags=re.IGNORECASE) is not None


def control_not_required(c):
    return re.match('.*(atac|dnase|dhs).*', re.sub('.*/', '', str(c)), flags=re.IGNORECASE) is not None


def _fastq_paths(config):
    return list(glob(os.path.join(config['fastq_dir'], '*.' + config['fastq_ext'])))


def _bams_paths(bams_dir):
    return list(glob(os.path.join(bams_dir, '*.bam')))


def fastq_names_wo_ext(fastq_paths):
    # file name w/o ext and parent folders: supports *.fastq and *.fastq.gz
    return [_split_to_fname_and_ext(f)[0] for f in fastq_paths]


def aligned_names(config, fastq_paths, bams_paths):
    if not bool(config['start_with_bams']):
        return fastq_aligned_names(config, fastq_paths)
    # file name w/o ext and parent folders: supports *.fastq and *.fastq.gz
    return [_split_to_fname_and_ext(f)[0] for f in bams_paths]


def fastq_aligned_names(config, fastq_paths):
    return _paired_fastq_samples_names(config, fastq_paths) + _single_fastq_samples_names(config, fastq_paths)


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


def _sample_by_file(config, fq_path):
    name = _split_to_fname_and_ext(fq_path)[0]
    if not bool(config['fastq_single_end_only']) and re.match('.*_R?[12]_?$', name):
        # paired file
        return re.sub('_R?[12]_?$', '', name)
    else:
        # single end file
        return name


def _single_fastq_samples_names(config, fastq_paths):
    if bool(config['fastq_single_end_only']):
        return fastq_names_wo_ext(fastq_paths)
    paired_samples = _paired_fastq_samples_names(config, fastq_paths)
    return [name for name in fastq_names_wo_ext(fastq_paths)
            if not re.match('.*_R?[12]_?$', name) or re.sub('_R?[12]_?$', '', name) not in paired_samples]


def _paired_fastq_samples_names(config, fastq_paths):
    if bool(config['fastq_single_end_only']):
        return []
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


def bowtie2_input_paths(config, paired):
    if bool(config['trim_reads']):
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


def trimmed_fastq_sample_names(fastq_paths):
    # here `name` could have _1 and _2 suffix in case of paired reads
    return [f"{name}_trimmed" for name in fastq_names_wo_ext(fastq_paths)]


#################
# Control reads #
#################

def _sample_2_control(config, fastq_paths, bams_paths):
    result = {}
    paths = fastq_paths if not bool(config['start_with_bams']) else bams_paths
    for path in paths:
        ext = _split_to_fname_and_ext(path)[1]
        control_path = find_control_for(path, ext)
        control_sample = _sample_by_file(config, control_path) if control_path else None
        result[_sample_by_file(config, path)] = control_sample
    return result


def find_control_for(file, ext="bam"):
    if is_control(file) or control_not_required(file):
        return None
    # Find all the files within folder
    controls = [os.path.basename(n) for n in glob(f'{os.path.dirname(file)}/*.{ext}') if is_control(n)]
    return _lcs_or_parts(os.path.basename(file), controls)


def _lcs_or_parts(bam_name, controls):
    if len(controls) == 0:
        return ''
    # Compute lcs for _ and - separated names, otherwise compute for all symbols
    lcs_parts = [_lcs(re.split('[\-_\s\.]+', str(bam_name.lower())),
                      re.split('[\-_\s\.]+', x.lower())) for x in controls]
    cps = [i for i in range(len(lcs_parts)) if lcs_parts[i] == max(lcs_parts)]
    if len(cps) == 1:
        return controls[cps[0]]
    # Otherwise return most similar for string
    return max(controls, key=lambda x: _lcs(str(bam_name.lower()), x.lower())) if len(controls) > 0 else ''


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


####################
# MACS2 parameters #
####################


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


# Small test
if __name__ == '__main__':
    print(_lcs_or_parts('GSM646340_H1_H3K36me3_rep2.bam', [
        'GSM646390_HMEC_Input_rep2.bam',
        'GSM646351_H1_Input_rep1.bam',
        'GSM646352_H1_Input_rep2.bam',
    ]))
