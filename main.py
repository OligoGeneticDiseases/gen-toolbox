from pathlib import Path
import argparse
import hail as hl
from pyspark import SparkConf, SparkContext

from CommandFactory import CommandFactory
from CommandHandler import CommandHandler

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

                if str.lower(args.command) == "readvcfs":  # Handle readvcfs command
                    sc = SparkContext(conf=conf)
                    hl.init(backend="spark", sc=sc, min_block_size=128)
                    ch.handle_read_vcfs_command()
                elif str.lower(args.command) == "loaddb":  # Handle readvcfs command
                    conf.set("spark.local.dir", "{0}".format(args.out))
                    sc = SparkContext(conf=conf)
                    hl.init(backend="spark", sc=sc, min_block_size=128)
                    ch.handle_load_db_command()

    except Exception as exp:
        print("Quitting.")
        raise
