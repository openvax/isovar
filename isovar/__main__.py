# PYTHON_ARGCOMPLETE_OK
"""Allow the same CLI through ``python -m isovar``."""

from .cli.commands import run

if __name__ == "__main__":
    run()
