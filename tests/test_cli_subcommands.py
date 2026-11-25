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
JSON_COMMANDS = {"rna-origin", "rna-features", "variant-rna-paths", "fusion-rna", "sv-rna", "allele-interpretations", "sv-audit"}

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
    expected_input = ("--hypotheses" if name == "rna-origin" else "--features" if name == "rna-features" else "--manifest" if name in {"calibrate-phasing", "sv-audit", "cohort-report", "normal-evidence"} else "--genome-fasta" if name == "viral-rna" else "reconstruct" if name == "isoforms" else "--sample-id" if name == "antigen-report"
                      else "--input" if name in JSON_COMMANDS | {"rna-proteome"} else "--vcf")
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


def test_entry_point_is_single_unified_isovar_cli():
    config = Path("pyproject.toml").read_text()
    entries = dict(re.findall(r'^([a-z-]+) = "(isovar\.cli\.[^" ]+)"$', config, re.M))
    assert entries == {"isovar": "isovar.cli.commands:run"}


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


def _leaf_parsers():
    """Every command's actual parsers, including nested subcommands such as `isoforms reconstruct`."""
    import argparse
    leaves = {}

    class Captured(Exception):
        pass

    def walk(parser, path):
        children = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
        if not children:
            leaves[" ".join(path)] = parser
        for action in children:
            for name, child in action.choices.items():
                walk(child, path + [name])

    def capture(self, args=None, namespace=None):
        walk(self, [self.prog])
        raise Captured

    original = argparse.ArgumentParser.parse_known_args
    argparse.ArgumentParser.parse_known_args = capture
    try:
        for name, (module, _) in commands.COMMANDS.items():
            try:
                import_module("isovar.cli." + module).run(["--help"], prog="isovar " + name)
            except (Captured, SystemExit):
                pass
    finally:
        argparse.ArgumentParser.parse_known_args = original
    return leaves


def test_every_option_is_documented_hyphenated_and_states_its_default():
    """#660: options keep help text, one flag spelling style and visible defaults."""
    import argparse
    leaves = _leaf_parsers()
    assert {path.split()[1] for path in leaves} == set(commands.COMMANDS)
    problems = []
    for path, parser in leaves.items():
        for action in parser._actions:
            if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
                continue
            if action.help == argparse.SUPPRESS:
                continue
            name = "/".join(action.option_strings) or action.dest
            if not action.help:
                problems.append("%s %s: no help" % (path, name))
            for option in action.option_strings:
                if option.startswith("--") and (option != option.lower() or "_" in option):
                    problems.append("%s %s: not lowercase-hyphenated" % (path, option))
            trivial = action.default in (None, False, True, argparse.SUPPRESS, [], ())
            if action.option_strings and not trivial and "default" not in (action.help or "").lower():
                problems.append("%s %s: default %r not stated" % (path, name, action.default))
    assert not problems, "\n".join(problems)


def test_every_command_controls_its_log_level():
    leaves = _leaf_parsers()
    missing_level = [
        path for path, parser in leaves.items()
        if not any("--log-level" in action.option_strings for action in parser._actions)
    ]
    missing_quiet = [
        path for path, parser in leaves.items()
        if not any(
            "-q" in action.option_strings and "--quiet" in action.option_strings
            for action in parser._actions
        )
    ]
    assert not missing_level
    assert not missing_quiet


def test_quiet_flag_raises_log_level_to_warning():
    parser = commands.load_subcommand_parser("allele-counts")
    default_args = parser.parse_args(["--bam", "rna.bam"])
    assert default_args.log_level == "INFO" and default_args.quiet is False
    quiet_args = parser.parse_args(["--bam", "rna.bam", "-q"])
    assert quiet_args.log_level == "WARNING" and quiet_args.quiet is True
    err_quiet_args = parser.parse_args(["--bam", "rna.bam", "--log-level", "ERROR", "--quiet"])
    assert err_quiet_args.log_level == "ERROR" and err_quiet_args.quiet is True


def test_output_and_short_flag_conventions_are_consistent():
    """#655, #658, #663: -o/--output is always a file, --output-dir a dir, short flags minimal."""
    import argparse

    allowed_short_flags = {"-h", "-o", "-f", "-q"}
    deprecated_visible = {
        "--plot-dir",
        "--summary-json",
        "--evidence-json",
        "--per-call",
        "--reference-fasta",
        "--library-cache-dir",
        "--sv-scratch-dir",
    }
    leaves = _leaf_parsers()
    problems = []
    for path, parser in leaves.items():
        for action in parser._actions:
            if action.help == argparse.SUPPRESS:
                continue
            for opt in action.option_strings:
                if opt.startswith("-") and not opt.startswith("--") and opt not in allowed_short_flags:
                    problems.append("%s: unexpected short flag %s" % (path, opt))
                if opt in deprecated_visible:
                    problems.append("%s: deprecated flag %s should be a hidden alias" % (path, opt))
            if "--output" in action.option_strings:
                if "-o" not in action.option_strings:
                    problems.append("%s: --output missing -o short flag" % path)
                if action.metavar != "FILE":
                    problems.append("%s: --output metavar is %r, expected 'FILE'" % (path, action.metavar))
            if "--output-dir" in action.option_strings and action.metavar != "DIR":
                problems.append("%s: --output-dir metavar is %r, expected 'DIR'" % (path, action.metavar))
    assert not problems, "\n".join(problems)


def test_version_option_has_explicit_help_text():
    import argparse

    parser = commands.build_parser()
    version_actions = [a for a in parser._actions if isinstance(a, argparse._VersionAction)]
    assert len(version_actions) == 1
    assert "version" in (version_actions[0].help or "").lower()
    assert "isovar" in (version_actions[0].help or "").lower()


@pytest.mark.parametrize("alias,canonical", sorted(commands.COMMAND_ALIASES.items()))
def test_deprecated_command_aliases_dispatch_with_stderr_warning(alias, canonical, monkeypatch, capsys):
    calls = []
    module_name = commands.COMMANDS[canonical][0]
    handler = import_module("isovar.cli." + module_name)
    monkeypatch.setattr(handler, "run", lambda args, **kw: calls.append((args, kw)))
    commands.run([alias, "--help"])
    captured = capsys.readouterr()
    assert captured.out == ""
    assert ("'isovar %s' is deprecated; use 'isovar %s'" % (alias, canonical)) in captured.err
    assert calls == [(["--help"], {"prog": "isovar " + canonical})]


@pytest.mark.parametrize(
    "command,legacy_args,canonical_args,attr",
    [
        ("fusion-rna", ["--input", "in.json", "-o", "out.json", "--plot-dir", "figs"], ["--input", "in.json", "--output", "out.json", "--output-dir", "figs"], "plot_dir"),
        ("antigen-report", ["--sample-id", "s1", "-o", "out.pdf", "--summary-json", "s.json"], ["--sample-id", "s1", "--output", "out.pdf", "--output-json", "s.json"], "summary_json"),
        ("evidence-report", ["--sample-id", "s1", "--bam", "a.bam", "-o", "out.pdf", "--per-call", "c.tsv", "--evidence-json", "e.json"], ["--sample-id", "s1", "--bam", "a.bam", "--output", "out.pdf", "--output-tsv", "c.tsv", "--output-json", "e.json"], "per_call"),
        ("rna-proteome", ["--input", "in.json", "--reference-context", "ref.json", "-o", "out.json", "--fasta", "out.fa"], ["--input", "in.json", "--reference-context", "ref.json", "--output", "out.json", "--output-fasta", "out.fa"], "fasta"),
        ("sv-rna", ["--input", "in.json", "--bam", "a.bam", "-o", "out.json", "--sv-scratch-dir", "/tmp"], ["--input", "in.json", "--bam", "a.bam", "--output", "out.json", "--scratch-dir", "/tmp"], "sv_scratch_dir"),
        ("sv-audit", ["--manifest", "m.json", "--output-dir", "out", "--library-cache-dir", "cache"], ["--manifest", "m.json", "--output-dir", "out", "--library-cache", "cache"], "library_cache_dir"),
        ("variant-rna-paths", ["--input", "in.json", "--bam", "a.bam", "--source", "L1", "-o", "out.json", "--reference-fasta", "ref.fa"], ["--input", "in.json", "--bam", "a.bam", "--source", "L1", "--output", "out.json", "--genome-fasta", "ref.fa"], "reference_fasta"),
        ("rna-origin", ["--bam", "a.bam", "--hypotheses", "h.json", "--sample-id", "s1", "--source-id", "L1", "-o", "out.json"], ["--bam", "a.bam", "--hypotheses", "h.json", "--sample-id", "s1", "--source", "L1", "--output", "out.json"], "source_id"),
        ("rna-features", ["--bam", "a.bam", "--features", "f.json", "--sample-id", "s1", "--source-id", "L1", "-o", "out.json"], ["--bam", "a.bam", "--features", "f.json", "--sample-id", "s1", "--source", "L1", "--output", "out.json"], "source_id"),
    ],
)
def test_hidden_flag_aliases_match_canonical_flags(command, legacy_args, canonical_args, attr):
    parser = commands.load_subcommand_parser(command)
    legacy = parser.parse_args(legacy_args)
    canonical = parser.parse_args(canonical_args)
    assert getattr(legacy, attr) == getattr(canonical, attr)
    assert legacy.output == canonical.output if hasattr(legacy, "output") else True


@pytest.mark.parametrize(
    "command,argv",
    [
        ("rna-origin", ["--bam", "a.bam", "--hypotheses", "h.json", "--sample-id", "s", "--source", "L", "-o", "out.json", "-f"]),
        ("rna-features", ["--bam", "a.bam", "--features", "f.json", "--sample-id", "s", "--source", "L", "-o", "out.json", "--force"]),
        ("variant-rna-paths", ["--bam", "a.bam", "--input", "i.json", "--source", "L", "-o", "out.json", "-f"]),
        ("viral-rna", ["--bam", "a.bam", "--genome-fasta", "g.fa", "--sample-id", "s", "--output-dir", "out", "-f"]),
        ("normal-evidence", ["--manifest", "m.json", "-o", "out.json", "-f"]),
        ("isoforms", ["reconstruct", "--bam", "a.bam", "--gene", "MAGEA1", "--sample-id", "s", "--genome-fasta", "g.fa", "--gtf", "g.gtf", "--annotation-version", "47", "--output-dir", "out", "-f"]),
    ],
)
def test_force_flag_is_accepted_on_overwrite_guarded_commands(command, argv):
    parser = commands.load_subcommand_parser(command)
    parsed = parser.parse_args(argv)
    assert parsed.force is True


@pytest.mark.parametrize("typo, suggestion", [("--verison", "--version"), ("--vfc", "--vcf"), ("--bogus", None)])
def test_mistyped_top_level_options_are_named_not_sent_to_the_pipeline(typo, suggestion, capsys):
    with pytest.raises(SystemExit) as error:
        commands.run([typo])
    assert error.value.code == 2
    err = capsys.readouterr().err
    assert err.startswith("isovar: error: unrecognized option '%s'" % typo) and "--bam" not in err
    assert ("did you mean '%s'" % suggestion in err) if suggestion else "did you mean" not in err


def test_help_command_shows_the_listing_or_a_command_help(capsys):
    with pytest.raises(SystemExit) as error:
        commands.run(["help"])
    assert error.value.code == 0 and "Pipeline:" in capsys.readouterr().out
    with pytest.raises(SystemExit) as error:
        commands.run(["help", "fusion-rna"])
    assert error.value.code == 0 and capsys.readouterr().out.startswith("usage: isovar fusion-rna ")
