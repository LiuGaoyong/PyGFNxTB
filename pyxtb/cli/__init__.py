# ruff: noqa: D103
"""Run XTB as Command Line Interface(CLI)."""

from sys import argv

from pyxtb.exe import run_xtb


def __run_xtb() -> str:
    """Run XTB as CLI."""
    _, out, err, is_success, _ = run_xtb(*argv[1:])
    if not is_success:
        raise RuntimeError(err)
    else:
        return out


def main() -> None:
    raise SystemExit(__run_xtb())
