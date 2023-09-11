import argparse
import sys
from src.cli.command_factory import CommandFactory
from src.cli.command_methods import CommandHandler

def main():
    parser = argparse.ArgumentParser(description='Genetic toolbox for variant analysis')
    parser.add_argument('command', type=str, help='Command to execute')
    parser.add_argument('args', nargs='*', help='Arguments for the command')
    args = parser.parse_args()

    command = CommandFactory.create_command(args.command, args.args)
    if command is None:
        print(f'Unknown command: {args.command}')
        sys.exit(1)

    handler = CommandHandler()
    handler.handle_command(command)


if __name__ == '__main__':
    main()