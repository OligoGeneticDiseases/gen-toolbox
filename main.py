import argparse
import os.path

import hail as hl

import file_utility


def vcfs_to_matrixtable(f, destination, write=True):
    if type(f) is list:
        files = list()
        with open(f) as vcflist:
            for vcfpath in vcflist:
                assert os.path.exists(vcfpath)
                files.append(vcfpath)
        f = files  # Send the list of .VCF paths to the importer

    else:
        assert os.path.exists(f), "Path {} does not exist.".format(f)

    table = hl.import_vcf(f, force=True, reference_genome='GRCh38')
    if write:
        table.write(destination)
    return table


if __name__ == '__main__':
    try:

        parser = argparse.ArgumentParser(prog="Annotation pipeline command-line file tool.")
        subparsers = parser.add_subparsers(title="commands", dest="command")
        findtype = subparsers.add_parser("Findtype", help="Find all specific files of a given filetype.")
        findtype.add_argument("-s", "--source", help="Directory to be searched.", action="store", type=str)
        findtype.add_argument("-d", "--directory", help="Directory to be saved to.", nargs='?', const=".",
                              action="store", type=str)         # TODO: Default value returns None type
        findtype.add_argument("-t", "--type", help="Filetype to be used.", action="store", type=str)
        findtype.add_argument("-r", "--regex", help="Filter strings by including parsed regex.", action="store",
                              type=str)
        readvcfs = subparsers.add_parser("Readvcfs", help="Turn VCF file(s) into a Hail MatrixTable.")
        readvcfs.add_argument("-f", "--file",
                              help="The VCF file to be parsed or file containing VCF pathnames on each line, LF (Unix)",
                              action="store", type=str)  # TODO: Default value returns None type
        readvcfs.add_argument("-d", "--dest", help="Destination folder to write the Hail MatrixTable files.",
                              nargs='?', const=os.path.abspath("."))

        args = parser.parse_args()
        if args.command is not None:

            if str.lower(args.command) == "findtype":
                # print(findtypes(args.directory, args.type))
                files = file_utility.find_filetype(args.source, args.type)

                file_utility.write_filelist(args.directory, "{0}.{1}.txt".format(os.path.basename(
                    os.path.normpath(args.source)), args.type), files, regex=args.regex)
                # Convert the directory into a name for the file, passing found files with
                # Regex in files matching only with a matching regex (e.g. *.vep.vcf wildcard).
                # Unique files only, duplicates written to duplicates_*.txt
            elif str.lower(args.command) == "readvcfs":
                print("Turning {0} into MatrixTable in directory {1}".format(args.file, args.dest))
                mt = vcfs_to_matrixtable(args.file, args.dest)
                mt.show()

        else:
            parser.print_usage()
            print("Invalid input, quitting.")

    except KeyboardInterrupt:
        print("Quitting.")
        raise

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
