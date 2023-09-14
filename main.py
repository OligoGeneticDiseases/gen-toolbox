from pathlib import Path
import argparse
import hail as hl
from pyspark import SparkConf, SparkContext
import sys

from src.cli.command_factory import CommandFactory
from src.cli.command_methods import CommandHandler

hail_home = Path(hl.__file__).parent.__str__()
if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description=f"Gnomad frequency table burden analysis pipeline command-line "
                                                     f"tool using Hail\n{hl.cite_hail()}")
        # Command Factory Init
        cf = CommandFactory(parser=parser)
        cf.create_find_type_command()
        cf.create_read_vcfs_command()
        cf.create_load_db_command()
        cf.create_check_relatedness_command()
        cf.create_pca_command()
        args = cf.parser.parse_args()

        if args.command is None:
            parser.print_usage()
        else:
            # Command Handler init
            ch = CommandHandler(args=args)

            if str.lower(ch.args.command) == "findtype":  # Handle findtype command
                ch.handle_find_type_command()
            else:
                conf = SparkConf()
                conf.set('spark.sql.files.maxPartitionBytes', '60000000000')
                conf.set('spark.sql.files.openCostInBytes', '60000000000')
                conf.set('spark.submit.deployMode', u'client')
                conf.set('spark.app.name', u'HailTools-TSHC')
                conf.set('spark.executor.memory', "4g")
                conf.set('spark.driver.memory', "56g")
                conf.set("spark.jars", "{0}/backend/hail-all-spark.jar".format(hail_home))
                conf.set("spark.executor.extraClassPath", "./hail-all-spark.jar")
                conf.set("spark.driver.extraClassPath", "{0}/backend/hail-all-spark.jar".format(hail_home))
                conf.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
                conf.set("spark.kryo.registrator", "is.hail.kryo.HailKryoRegistrator")
                conf.set("spark.driver.bindAddress", "127.0.0.1")

                # TODO: import spark config from json file
                # TODO: Add a temp dir option: https://discuss.hail.is/t/hail-doesnt-respect-tmp-dir/2319

                if str.lower(args.command) == "readvcfs":  # Handle readvcfs command
                    sc = SparkContext(conf=conf)
                    conf.set("spark.local.dir", "{0}".format(args.dest))
                    hl.init(backend="spark", sc=sc, min_block_size=4096, tmp_dir=args.temp, local_tmpdir=args.temp)
                    ch.handle_read_vcfs_command()
                elif str.lower(args.command) == "loaddb":  # Handle readvcfs command
                    conf.set("spark.local.dir", "{0}".format(args.dest))
                    sc = SparkContext(conf=conf)
                    hl.init(backend="spark", sc=sc, min_block_size=64)
                    ch.handle_load_db_command()
                elif str.lower(args.command) == "relatedness2":
                    sc = SparkContext(conf=conf)
                    conf.set("spark.local.dir", "{0}".format(args.dest))
                    hl.init(backend="spark", sc=sc, min_block_size=4096, tmp_dir=args.temp, local_tmpdir=args.temp)
                    ch.handle_check_relatedness()
                    hl.utils.info("Finished with relatedness2 command!")
                elif str.lower(args.command) == "pca":
                    ch.handle_pca()

    except Exception as e:
        print("ERROR: " + str(e))
        sys.exit(1)
