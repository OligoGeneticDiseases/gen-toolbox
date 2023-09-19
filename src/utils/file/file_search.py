import os

try:
    from scandir import walk
except ImportError:
    from os import walk


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
    for dirpath, dirnames, files in os.walk(dir):
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
