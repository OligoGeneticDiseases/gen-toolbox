import sys
from src.cli.command_setup import setup_parser, command_handlers, setup_spark_config


def main():
    try:
        parser = setup_parser()
        args = parser.parse_args()
        if args.command is None:
            parser.print_usage()
            sys.exit(1)

        conf = setup_spark_config(args)
        handlers = command_handlers(args, conf)
        command = str.lower(args.command)
        if command in handlers:
            handlers[command]()
        else:
            print(f"Unknown command: {command}")
            sys.exit(1)

    except Exception as e:
        print("ERROR: " + str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
