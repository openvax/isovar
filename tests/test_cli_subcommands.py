"""Subcommands reuse legacy parsers/handlers without changing their defaults."""

import inspect
from importlib import import_module
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

from isovar import __version__
from isovar.cli import commands


# Structured-input commands: no --vcf and no hyphenated legacy script.
JSON_COMMANDS = {"variant-rna-paths", "fusion", "sv-rna", "allele-interpretations", "sv-audit"}

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
    expected_input = ("--manifest" if name in {"calibrate-phasing", "sv-audit", "cohort", "normal-evidence"} else "--genome-fasta" if name == "viral" else "reconstruct" if name == "isoforms" else "--gene" if name == "uorf-evidence" else "import-sid" if name == "audit-bundle" else "--sample-id" if name == "antigen-report"
                      else "--input" if name in JSON_COMMANDS | {"proteome"} else "--vcf")
    assert expected_input in help_text
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
                       for name, (module, _) in commands.COMMANDS.items() if name not in {"run", *JSON_COMMANDS}}


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


def test_every_command_is_in_exactly_one_help_group():
    grouped = [name for _, names in commands.COMMAND_GROUPS for name in names]
    assert sorted(grouped) == sorted(commands.COMMANDS)


def test_top_level_help_lists_commands_by_group_with_run_first(capsys, monkeypatch):
    monkeypatch.setenv("COLUMNS", "200")
    commands.run([])
    help_text = capsys.readouterr().out
    positions = [help_text.index(title + ":") for title, _ in commands.COMMAND_GROUPS]
    assert positions == sorted(positions)
    assert commands.COMMAND_GROUPS[0][1][0] == "run"
    # Long descriptions wrap, so compare with whitespace collapsed.
    collapsed = " ".join(help_text.split())
    for name, (_, description) in commands.COMMANDS.items():
        assert "%s %s" % (name, description) in collapsed
    assert "SUPPRESS" not in help_text


def test_top_level_help_fits_an_80_column_terminal(capsys, monkeypatch):
    monkeypatch.setenv("COLUMNS", "80")
    commands.run([])
    lines = capsys.readouterr().out.splitlines()
    assert max(len(line) for line in lines) <= 80


def test_table_exports_group_is_exactly_the_commands_that_write_tables():
    def writes_table(name):
        source = inspect.getsource(import_module("isovar.cli." + commands.COMMANDS[name][0]))
        return "run_dataframe_command" in source or "add_output_args" in source
    groups = dict(commands.COMMAND_GROUPS)
    assert all(writes_table(name) for name in groups["Table exports"])
    others = [n for title, names in commands.COMMAND_GROUPS if title != "Table exports" for n in names]
    assert not [n for n in others if writes_table(n) and n != "run"]


def test_unknown_command_suggests_the_closest_name_without_listing_all(capsys):
    with pytest.raises(SystemExit) as error:
        commands.run(["variant-sequnces"])
    assert error.value.code == 2
    err = capsys.readouterr().err
    assert "invalid choice: 'variant-sequnces' (did you mean 'variant-sequences'?)" in err
    assert "isovar --help" in err
    assert "sv-audit" not in err
    with pytest.raises(SystemExit):
        commands.run(["zzzzzz"])
    assert "did you mean" not in capsys.readouterr().err


def test_command_line_doc_lists_every_command_in_help_group_order():
    doc = (Path(__file__).parent.parent / "docs" / "command-line.md").read_text()
    for title, names in commands.COMMAND_GROUPS:
        section = doc[doc.index("### " + title + "\n"):]
        positions = [section.index("`isovar %s`" % name) for name in names]
        assert positions == sorted(positions)
