"""Subcommands reuse legacy parsers/handlers without changing their defaults."""

from importlib import import_module
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

from isovar import __version__
from isovar.cli import commands


@pytest.mark.parametrize("name", commands.COMMANDS)
def test_dispatch_passes_options_unchanged_to_the_existing_handler(name, monkeypatch):
    calls = []
    module_name = commands.COMMANDS[name][0]
    handler = import_module("isovar.cli." + module_name)
    monkeypatch.setattr(handler, "run", lambda args, **kw: calls.append((args, kw)))
    options = ["--bam", "plot", "--output", "run", "--min-mapping-quality", "7"]
    commands.run([name] + options)
    assert calls == [(options, {"prog": "isovar " + name})]


@pytest.mark.parametrize("name", commands.COMMANDS)
def test_subcommand_help_uses_correct_program_and_preserves_parser(name, capsys):
    handler = import_module("isovar.cli." + commands.COMMANDS[name][0])
    original = handler.parser.prog if hasattr(handler, "parser") else None
    with pytest.raises(SystemExit) as error:
        commands.run([name, "--help"])
    assert error.value.code == 0
    help_text = capsys.readouterr().out
    assert help_text.startswith("usage: isovar " + name + " ")
    assert ("--input" if name in {"fusion", "sv-rna"} else "--vcf") in help_text
    if original is not None:
        assert handler.parser.prog == original


@pytest.mark.parametrize("name", commands.COMMANDS)
def test_subcommand_errors_keep_exit_status_and_program(name, capsys):
    with pytest.raises(SystemExit) as error:
        commands.run([name, "--not-an-option"])
    assert error.value.code == 2
    assert ("isovar " + name + ": error:") in capsys.readouterr().err


def test_option_first_interface_is_unchanged(monkeypatch):
    from isovar.cli import isovar_main

    args = ["--bam", "plot", "--output", "run"]
    calls = []
    monkeypatch.setattr(isovar_main, "run", lambda args, **kw: calls.append((args, kw)))
    commands.run(args)
    assert calls == [(args, {"prog": "isovar"})]


def test_root_help_version_unknown_command_and_no_args(capsys):
    commands.run([])
    text = capsys.readouterr().out
    assert all(name in text for name in commands.COMMANDS)
    for args in [["--help"], ["--version"]]:
        with pytest.raises(SystemExit) as error:
            commands.run(args)
        assert error.value.code == 0
    assert __version__ in capsys.readouterr().out
    with pytest.raises(SystemExit) as error:
        commands.run(["plt"])
    assert error.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_only_selected_handler_is_imported(monkeypatch):
    imported = []

    def load(name, package):
        imported.append((name, package))
        return SimpleNamespace(run=lambda args, **kw: None)

    monkeypatch.setattr(commands, "import_module", load)
    commands.run(["allele-counts", "--bam", "rna.bam"])
    assert imported == [(".isovar_allele_counts", "isovar.cli")]


def test_entry_points_retain_all_legacy_aliases():
    config = Path("pyproject.toml").read_text()
    entries = dict(re.findall(r'^([a-z-]+) = "(isovar\.cli\.[^" ]+)"$', config, re.M))
    assert entries.pop("isovar") == "isovar.cli.commands:run"
    assert entries == {"isovar-" + name: "isovar.cli." + module + ":run"
                       for name, (module, _) in commands.COMMANDS.items() if name not in {"run", "fusion", "sv-rna"}}


def test_module_help_does_not_load_plotting_or_command_handlers():
    code = (
        "import sys; from isovar.cli.commands import run; run([]); "
        "assert 'matplotlib' not in sys.modules; "
        "assert not any(x.startswith('isovar.cli.isovar_') for x in sys.modules)"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    result = subprocess.run([sys.executable, "-m", "isovar", "--version"], capture_output=True, text=True)
    assert result.returncode == 0 and result.stdout.strip() == "isovar " + __version__


@pytest.mark.parametrize("name", ["run", "allele-counts", "reference-contexts"])
def test_real_subcommand_csv_matches_legacy_handler(name, tmp_path):
    from .test_cli import vcf_args, args_with_bam

    args = vcf_args if name == "reference-contexts" else args_with_bam
    module_name = commands.COMMANDS[name][0]
    legacy, subcommand = tmp_path / "legacy.csv", tmp_path / "subcommand.csv"
    import_module("isovar.cli." + module_name).run(args + ["--output", str(legacy)])
    commands.run([name] + args + ["--output", str(subcommand)])
    assert legacy.read_bytes() == subcommand.read_bytes()
