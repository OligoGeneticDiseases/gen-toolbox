import argparse
import json
from pathlib import Path

import hail as hl
from pyspark import SparkContext, SparkConf
from src.cli.command_factory import CommandFactory
from src.cli.command_methods import CommandHandler

hail_home = Path(hl.__file__).parent.__str__()


def setup_parser():
    parser = argparse.ArgumentParser(
        description=f"Gnomad frequency table burden analysis pipeline command-line tool using Hail\n{__import__('hail').cite_hail()}"
    )
    cf = CommandFactory(parser=parser)
    cf.create_find_type_command()
    cf.create_read_vcfs_command()
    cf.create_load_db_command()
    cf.create_check_relatedness_command()
    cf.create_pca_command()
    return parser


def setup_spark_config(args):
    path = Path(__file__).parent / "../config/spark_conf.json"
    with path.open() as f:
        conf_data = json.load(f)

    conf_data["spark.jars"] = conf_data["spark.jars"].format(hail_home=hail_home)
    conf_data["spark.driver.extraClassPath"] = conf_data[
        "spark.driver.extraClassPath"
    ].format(hail_home=hail_home)

    conf = SparkConf().setAll(conf_data.items())
    return conf


def command_handlers(args, conf):
    handlers = {
        "findtype": CommandHandler(args).handle_find_type_command,
        "readvcfs": lambda: init_spark_and_run(
            args, conf, CommandHandler(args).handle_read_vcfs_command
        ),
        "loaddb": lambda: init_spark_and_run(
            args, conf, CommandHandler(args).handle_load_db_command
        ),
        "relatedness2": lambda: init_spark_and_run(
            args, conf, CommandHandler(args).handle_check_relatedness
        ),
        "pca": CommandHandler(args).handle_pca,
    }
    return handlers


def init_spark_and_run(args, conf, func):
    sc = SparkContext(conf=conf)

    conf.set("spark.local.dir", f"{args.dest}")
    conf.set("spark.jars", f"{hail_home}/backend/hail-all-spark.jar")
    conf.set("spark.driver.extraClassPath", f"{hail_home}/backend/hail-all-spark.jar")

    if args.command == "readvcfs" or args.command == "relatedness2":
        hl.init(
            backend="spark",
            sc=sc,
            min_block_size=4096,
            tmp_dir=args.temp,
            local_tmpdir=args.temp,
        )
    else:
        hl.init(backend="spark", sc=sc, min_block_size=64)

    func()
