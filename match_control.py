import re
import os
import glob


def find_control_for(file, ext="bam"):
    bam_name = os.path.basename(file).lower()
    if 'input' in bam_name:
        return ''

    # Find all the files within folder
    parent_dir = os.path.dirname(file)
    names = [os.path.basename(n) for n in glob.glob(f'{parent_dir}/*.{ext}')]
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
