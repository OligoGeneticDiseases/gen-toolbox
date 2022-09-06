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
        findtype.add_argument("-d", "--directory", help="Directory to be searched.", action="store", type=str)
        findtype.add_argument("-t", "--type", help="Filetype to be used.", action="store", type=str)

        args = parser.parse_args()

        if str.lower(args.command) == "findtype":
            #print(findtypes(args.directory, args.type))
            excels = file_utility.find_filetype(args.directory, args.type)
            file_utility.write_filelist(".", "{0}.excels.txt".format(os.path.basename(
                os.path.normpath(args.directory))), excels)  # Convert the directory into a name-modifier
        else:
            print("Invalid command, quitting.")

    except KeyboardInterrupt:
        print("Quitting.")
        raise


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
