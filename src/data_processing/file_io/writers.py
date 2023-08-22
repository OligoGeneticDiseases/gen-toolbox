import os
import shutil
import re

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk

duplicates = 0

def write_filelist(dir, file_name, file_paths, include_duplicates=False, verbose=1, regex=None):
    '''
    This function writes lists of files that have been inputted as arrays.
    :param dir: Directory to be written to.
    :param file_name: Name of the file to be written to.
    :param file_paths: The filepaths as an array of tuples. Filename and filepath.
    :param include_duplicates: Whether to include duplicate files into the list. False by default.
    If False, write duplicate filename-paths into a separate file duplicate_.*
    :param verbose level: Verbosity of stdout. Prints filepaths that are looked at. Default 1. Show
    Total values but no lines in the terminal, 0 turns off all terminal response. TODO: Convert to logging
    :param regex: regex expression to be evaluated. For example, try to match only certain files containing a string.
    Regex is done after duplicate prefix check.
    :return: Returns tuple of (<all paths written to the main file>, <duplicate files written out or empty list>)
    '''
    try:
        regex_rules = None
        if regex is not None:
            # It is more efficient to compile once and reuse the regex object.
            regex_rules = re.compile(regex, flags=re.I)
            #print(regex_rules)
        with open(os.path.join(dir, file_name), "w") as f:
            all_paths = list()
            unique_prefixes = list()
            duplicate_prefixes = list()
            regexed_out = 0
            file_paths.sort(key=lambda pair: pair[0])
            for file_path_pair in file_paths:
                full_path = file_path_pair[1]
                all_paths.append(full_path)
                prefix = trim_prefix(file_path_pair[0])  # E0000000 /trimmed by "." and "_"
                # Is unique
                if prefix not in unique_prefixes:
                    if regex is not None:
                        if eval_regex(full_path,
                                      regex_rules) is not None:  # there is a match with the corresponding regex
                            unique_prefixes.append(prefix)  # it is unique and added to the list
                            f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                            if verbose > 1:
                                print("Unique regexed: {0}".format(full_path))
                        else:
                            regexed_out += 1
                            if verbose > 1:
                                print("No regex match for: {0}".format(full_path))
                    else:
                        unique_prefixes.append(prefix)  # it is unique and added to the list
                        f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                        if verbose > 1:
                            print("Unique unregexed: {0}".format(full_path))
                    all_paths.append(full_path)
                else:
                    duplicate_prefixes.append([prefix, file_path_pair[0], file_path_pair[1]])
                    if include_duplicates:
                        if regex is not None:
                            if eval_regex(full_path,
                                          regex_rules) is not None:  # there is a match with the corresponding regex
                                f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                                if verbose > 1:
                                    print("Not unique regexed: {0}".format(full_path))

                        else:
                            f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                            if verbose > 1:
                                print("Not unique unregexed: {0}".format(full_path))
                        all_paths.append(full_path)


            if verbose > 0:
                print("Writing list to {}".format(os.path.join(dir, file_name)))
                print("Total lines: {}".format(len(file_paths)))
                print("Duplicate prefixes: {}".format(len(duplicate_prefixes)))
                print("Regex matches {0}, regex filtered out {1} lines.".format(len(unique_prefixes), regexed_out))
            if not include_duplicates:  # write the duplicates to a separate file
                duplicate_filename = os.path.join(dir, "duplicates_" + file_name)
                if len(duplicate_prefixes) > 0:
                    if verbose > 0:
                        print("Creating duplicate sample list in {0}".format(duplicate_filename))
                    f_duplicates = open(duplicate_filename, "w")
                    for duplicate_name in duplicate_prefixes:
                        f_duplicates.write("\t".join(duplicate_name) + "\n")

        return all_paths, duplicate_prefixes

    except (NameError, KeyboardInterrupt, SyntaxError, AssertionError):
        if NameError or AssertionError:
            print("Unable to find {}".format(dir))
            raise
        if KeyboardInterrupt:
            print("Quitting.")


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

def rename_file_idx(name_path, idx):
    basename = os.path.basename(name_path).split(".")
    basename[0] = basename[0] + "_re" + str(idx)
    new_name = ".".join(basename)
    return new_name

def write_vcfs_list(dir, output):
    return write_filelist(dir, output, find_vcfs(dir))

def write_bams_list(dir, output):
    return write_filelist(dir, output, find_bams(dir))

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

def write_filelist(input_dir, file_name, file_paths, include_duplicates=False, verbose=1, regex=None):
    '''
    This function writes lists of files that have been inputted as arrays.
    :param input_dir: Directory to be written to.
    :param file_name: Name of the file to be written to.
    :param file_paths: The filepaths as an array of tuples. Filename and filepath.
    :param include_duplicates: Whether to include duplicate files into the list. False by default.
    If False, write duplicate filename-paths into a separate file duplicate_.*
    :param verbose level: Verbosity of stdout. Prints filepaths that are looked at. Default 1. Show
    Total values but no lines in the terminal, 0 turns off all terminal response. TODO: Convert to logging
    :param regex: regex expression to be evaluated. For example, try to match only certain files containing a string.
    Regex is done after duplicate prefix check.
    :return: Returns tuple of (<all paths written to the main file>, <duplicate files written out or empty list>)
    '''
    try:
        regex_rules = None
        if regex is not None:
            regex_rules = re.compile(regex, flags=re.I)
        with open(os.path.join(input_dir, file_name), "w") as f:
            all_paths = list()
            unique_prefixes = list()
            duplicate_prefixes = list()
            regexed_out = 0
            file_paths.sort(key=lambda pair: pair[0])
            for file_path_pair in file_paths:
                full_path = file_path_pair[1]
                all_paths.append(full_path)
                prefix = trim_prefix(file_path_pair[0])  # E0000000 /trimmed by "." and "_"
                # Is unique
                if prefix not in unique_prefixes:
                    if regex is not None:
                        if eval_regex(full_path,
                                      regex_rules) is not None:  # there is a match with the corresponding regex
                            unique_prefixes.append(prefix)  # it is unique and added to the list
                            f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                            if verbose > 1:
                                print("Unique regexed: {0}".format(full_path))
                        else:
                            regexed_out += 1
                            if verbose > 1:
                                print("No regex match for: {0}".format(full_path))
                    else:
                        unique_prefixes.append(prefix)  # it is unique and added to the list
                        f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                        if verbose > 1:
                            print("Unique unregexed: {0}".format(full_path))
                    all_paths.append(full_path)
                else:
                    duplicate_prefixes.append([prefix, file_path_pair[0], file_path_pair[1]])
                    if include_duplicates:
                        if regex is not None:
                            if eval_regex(full_path,
                                          regex_rules) is not None:  # there is a match with the corresponding regex
                                f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                                if verbose > 1:
                                    print("Not unique regexed: {0}".format(full_path))

                        else:
                            f.write(prefix + "\t" + full_path + "\n")  # write lines only if duplicates included
                            if verbose > 1:
                                print("Not unique unregexed: {0}".format(full_path))
                        all_paths.append(full_path)

            if verbose > 0:
                print("Writing list to {}".format(os.path.join(input_dir, file_name)))
                print("Total lines: {}".format(len(file_paths)))
                print("Duplicate prefixes: {}".format(len(duplicate_prefixes)))
                print("Regex matches {0}, regex filtered out {1} lines.".format(len(unique_prefixes), regexed_out))
            if not include_duplicates:  # write the duplicates to a separate file
                duplicate_filename = os.path.join(input_dir, "duplicates_" + file_name)
                if len(duplicate_prefixes) > 0:
                    if verbose > 0:
                        print("Creating duplicate sample list in {0}".format(duplicate_filename))
                    f_duplicates = open(duplicate_filename, "w")
                    for duplicate_name in duplicate_prefixes:
                        f_duplicates.write("\t".join(duplicate_name) + "\n")

        return all_paths, duplicate_prefixes

    except (NameError, KeyboardInterrupt, SyntaxError, AssertionError):
        if NameError or AssertionError:
            print("Unable to find {}".format(input_dir))
            raise
        if KeyboardInterrupt:
            print("Quitting.")


def trim_prefix(filename, seperator="."):
    prefix_clean = filename.rsplit(".")[0].rsplit("_")[0]
    return prefix_clean


def eval_regex(text, regex):
    """
    Will return the parsed string with the defined regex rule
    :param text: Text to be parsed
    :param regex: Regex ruleset
    :return: Parsed string or None if not able to parse.
    """
    result = regex.search(text)
    return result