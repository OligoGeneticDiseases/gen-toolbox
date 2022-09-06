import argparse
import itertools
import os
from distutils.version import LooseVersion as Version

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk

import shutil
import sys

duplicates = 0

reflines = None


class Interval():
    def __init__(self, chrom, start, stop, symbol):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.symbol = symbol

    def between(self, range):
        if self.chrom == range[0]:
            if self.start <= range[1] and self.stop >= range[2]:
                return True
        return False

    def __str__(self):
        return "\t".join((self.chrom, str(self.start), str(self.stop), str(self.symbol)))

    def compare(self, other):
        return self.symbol == other.symbol and self.chrom == other.chrom


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
    :return: list of tuples
    of file name and file directory
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


def write_filelist(dir, file_name, file_paths, include_duplicates=False, verbose=False):
    try:
        with open(os.path.join(dir, file_name), "w") as f:
            all_paths = list()
            unique_prefixes = list()
            duplicate_prefixes = list()

            file_paths.sort(key=lambda pair: pair[0])
            for file_path_pair in file_paths:
                full_path = file_path_pair[1]
                all_paths.append(full_path)
                prefix = trim_prefix(file_path_pair[0])  # E0000000 /trimmed by "." and "_"

                # Is unique
                if prefix not in unique_prefixes:
                    unique_prefixes.append(prefix)  # it is unique and added to the list
                    f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included

                # Not unique
                else:
                    # But write anyway
                    if include_duplicates:
                        f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                    else:
                        if verbose:
                            print("Duplicate prefix {0}\t{1}".format(prefix, full_path))
                        # otherwise append to list, which will be printed if duplicates are to be excluded
                        duplicate_prefixes.append([file_path_pair[0], prefix, file_path_pair[1]])
            print("Writing filelist to {}".format(os.path.join(dir, file_name)))
            print("Total files: {}".format(len(file_paths)))
            print("Duplicate prefixes: {}".format(len(duplicate_prefixes)))
            if not include_duplicates:  # write the duplicates to a separate file
                duplicate_filename = os.path.join(dir, "duplicates_" + file_name)
                print("Creating duplicate sample list in {0}".format(duplicate_filename))
                f_duplicates = open(duplicate_filename, "w")
                for duplicate_name in duplicate_prefixes:
                    f_duplicates.write("\t".join(duplicate_name) + "\n")

        return all_paths

    except (NameError, KeyboardInterrupt, SyntaxError, AssertionError):
        if NameError or AssertionError:
            print("Unable to find {}".format(dir))
            raise
        if KeyboardInterrupt:
            print("Quitting.")


def trim_prefix(filename, seperator="."):
    prefix_clean = filename.rsplit(".")[0].rsplit("_")[0]
    return prefix_clean


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


def write_prefixes_list(dir, output):
    samples = find_vcfs(dir)
    samples.sort(key=lambda pair: pair[0])
    clean_prefixes = list()

    with open(output, "wb+") as prefixes:
        for prefix in samples:
            prefix_clean = prefix[0].rsplit('.')
            if len(prefix_clean) == 2:
                prefix_clean = prefix_clean[0]
                prefixes.write(prefix_clean + "\n")
                clean_prefixes.append(prefix_clean)

            else:
                pass
                # This file might be an intermediary vcf file

    return clean_prefixes


def write_vcfs_list(dir, output):
    return write_filelist(dir, output, find_vcfs(dir))


def write_bams_list(dir, output):
    return write_filelist(dir, output, find_bams(dir))


def rename_file_idx(name_path, idx):
    basename = os.path.basename(name_path).split(".")
    basename[0] = basename[0] + "_re" + str(idx)
    new_name = ".".join(basename)
    return new_name


def copy_vcf(files, dest, overwrite=False):
    try:
        os.makedirs(dest)
    except OSError:
        if not os.path.isdir(dest):
            raise
    nr = len(files)
    unique_copied = []
    renamed = []
    not_found = []

    print("Starting copying of {0} files to destination folder {1}".format(nr, dest))

    for path in files:
        path = path.rstrip()
        if not os.path.isfile(os.path.join(dest, os.path.basename(path))) or overwrite:
            if os.path.isfile(path):
                shutil.copy(path, dest)
                fname = os.path.join(dest, os.path.basename(path))
                unique_copied.append(fname)
            else:
                print("File does not exist: " + str(path))
                not_found.append(str(path))
        else:
            re_idx = 1

            while os.path.isfile(os.path.join(dest, rename_file_idx(path, re_idx))):
                re_idx += 1
            fname = rename_file_idx(path, re_idx)
            if os.path.isfile(path):
                shutil.copy(path, os.path.join(dest, fname))

            renamed.append(os.path.join(dest, fname))
    print("Finished copying of {0} files. Renamed {1} files. ".format(len(unique_copied), len(renamed)))
    return unique_copied, renamed, not_found


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

