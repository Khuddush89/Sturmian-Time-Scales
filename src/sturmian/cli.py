"""Command-line entry points for verification and figure generation."""

import argparse
from pathlib import Path

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="sturmian", description="Reproduce finite-grid and hybrid Sturmian computations."
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("verify", "run all exact and numerical assertions"),
        ("figures", "generate seven PNG figures after verification"),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument(
            "--output",
            type=Path,
            default=Path("outputs/reference"),
            help="output directory, relative to the current working directory",
        )
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            from .verification import run
        else:
            from .figures import run
        run(args.output)
    except (AssertionError, FileNotFoundError, ValueError) as error:
        parser.exit(1, f"sturmian: {error}\n")
    return 0
