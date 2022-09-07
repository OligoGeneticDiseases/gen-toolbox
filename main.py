import argparse
import os.path
import sys
import time

from pipeline_utility import file_utility

sys.path.append('pipeline_utility')


def findtypes(dir, type):
    start = time.time()
    excels = file_utility.find_filetype(args.directory, args.type)
    f = open("excels.txt", "w")
    for pair in excels:
        f.write("{0}\t{1}\n".format(pair[0], pair[1]))
    end = time.time()
    difference = end - start
    print("Created {0} with {1} lines in {2} seconds".format(f.name, len(excels), difference))
    return excels


if __name__ == '__main__':
    try:

        parser = argparse.ArgumentParser(prog="Annotation pipeline command-line file tool.")
        subparsers = parser.add_subparsers(title="commands", dest="command")
        findtype = subparsers.add_parser("Findtype", help="Find all specific files of a given filetype.")
        findtype.add_argument("-s", "--source", help="Directory to be searched.", action="store", type=str)
        findtype.add_argument("-d", "--directory", help="Directory to be saved to.", nargs='?', const=".",action="store", type=str)
        findtype.add_argument("-t", "--type", help="Filetype to be used.", action="store", type=str)
        findtype.add_argument("-r", "--regex", help="Filter strings by including parsed regex.", action="store",
                              type=str)

        args = parser.parse_args()
        if args.command is None:
            parser.print_usage()
            exit()

        if str.lower(args.command) == "findtype":
            #print(findtypes(args.directory, args.type))
            files = file_utility.find_filetype(args.source, args.type)

            file_utility.write_filelist(args.directory, "{0}.{1}.txt".format(os.path.basename(
                os.path.normpath(args.source)), args.type), files, regex=args.regex)
            # Convert the directory into a name for the file, passing found files with
            # Regex in files matching only with a matching regex (e.g. *.vep.vcf wildcard).
            # Unique files only, duplicates written to duplicates_*.txt
        else:
            print("Invalid command, quitting.")

    except KeyboardInterrupt:
        print("Quitting.")
        raise


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
