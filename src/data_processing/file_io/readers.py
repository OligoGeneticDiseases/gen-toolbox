import os

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk

import shutil
import re

duplicates = 0

def find_file(main_dir, filename):
    """
    Will find the first file with the exact name match in a dir tree. Otherwise returns None.
    :param main_dir: The upmost directory to start the search in, will walk through subdirectories
    :param filename: The filename to be searched.
    :return: tuple(filename, full path)
    """
    assert os.path.exists(main_dir), "Path {} does not exist.".format(main_dir)

    for (dirpath, dirnames, files) in walk(main_dir):
        for name in files:
            if name == filename:
                return name, os.path.join(dirpath, name)
    return None, None


def find_filetype(dir, filetype, findunique=False, verbose=True):
    """
    Will find all files of a certain type (e.g. .vcf or .bam files) in a directory. Method will enter every
    subdirectory. Can look for only a single filetype at a time.
    :param verbose: If verbose, report repeating filenames in terminal
    :param dir: String of directory to walk.
    :param filetype: String of filetype to search for (e.g. .vcf or .bam)
    :param findunique: By default find all files of a given type,
    setting it to True will skip files with the same name but different location
    :return: list of tuple of file name and file directory
    """
    assert os.path.exists(dir), "Path {} does not exist.".format(dir)
    duplicates = 0

    unique_files = list(())
    for (dirpath, dirnames, files) in walk(dir):
        for name in files:
            if name.endswith(filetype):
                if name not in map(lambda x: x[0], unique_files):
                    unique_files.append((name, os.path.join(dirpath, name)))
                else:
                    duplicates += 1
                    if verbose:
                        print("Duplicate filename {0}".format(name))
                    if not findunique:
                        # Append anyway if findunique is set to False
                        unique_files.append((name, os.path.join(dirpath, name)))
    print("Duplicate filenames: {0}".format(duplicates))
    return unique_files

def find_vcfs(dir):
    return find_filetype(dir, '.vcf')


def find_bams(dir):
    return find_filetype(dir, '.bam')

def find_type(dir, extension):
    return find_filetype(dir, extension)


def find_prefixes(dir, extension):
    """
    Finds a list of prefixes with a given file extension e.g. E00001.vcf.gz --> E00001 in a directory

    :param dir: Directory to search in.
    :param extension: File extension string, must not include the seperator '.' ("vcf" or "vcf.gz" not ".vcf")
    """
    samples = find_filetype(dir, extension)
    clean_prefixes = list()
    for prefix in samples:
        prefix_clean = prefix[0].rsplit(".")
        # We only want exact matches, so E00001.genome.vcf.gz would be skipped
        # This tiny math ensures that long extensions don't confuse us.
        if len(prefix_clean) < 2 + len(extension.split(".")):
            prefix_clean = prefix_clean[0]
            clean_prefixes.append(prefix_clean)
    return clean_prefixes





def file_len(fname):
    if os.path.exists(fname):
        with open(fname) as f:
            i, l = 0, 0
            for i, l in enumerate(f):
                pass
        return i + 1
    return None


def count_unique_names(infile, col, seperator="\t"):
    names = []
    if os.path.exists(infile):
        with open(infile) as f:
            for i, l in enumerate(f):
                cols = l.split(seperator)
                if not col > len(cols):
                    names.append(cols[col])
                else:
                    raise IndexError("Is your seperator correct? "
                                     "There weren't enough columns "
                                     "after splitting the line #{0}\n{1}!".format(i, l))
    unique = set(names)
    return len(unique)
