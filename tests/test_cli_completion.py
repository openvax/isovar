# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Shell tab completion (argcomplete), lazy CLI startup and subtle help styling (#684)."""

import io
import os
from pathlib import Path
import subprocess
import sys

import pytest

from isovar.cli import commands
from isovar.cli.completion import completing_subcommand, is_completing
from isovar.cli.help import HelpArgumentParser, color_enabled, style_help, style_table_header

REPO_ROOT = Path(__file__).resolve().parents[1]

# Heavy runtime libraries that must never load during CLI parser construction.
HEAVY_RUNTIME_MODULES = (
    "pandas",
    "numpy",
    "pysam",
    "pyarrow",
    "pyensembl",
    "matplotlib",
    "reportlab",
    "orfdata",
    "oncoref",
    "gtfparse",
    "biotite",
    "torch",
)

# Subcommands whose parser construction does not call varcode.cli.variant_args.
NON_VARCODE_SUBCOMMANDS = (
    "allele-interpretations",
    "antigen-report",
    "calibrate-phasing",
    "cohort-report",
    "fusion-rna",
    "isoforms",
    "normal-evidence",
    "rna-features",
    "rna-origin",
    "rna-proteome",
    "sv-audit",
    "sv-rna",
    "sv-source",
    "variant-rna-paths",
    "viral-rna",
)


def _subprocess_env(**extra):
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    parts = [str(REPO_ROOT)]
    if existing:
        parts.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(parts)
    env.update(extra)
    return env


def test_entry_modules_carry_argcomplete_marker_within_first_1024_bytes():
    for rel in ("isovar/cli/commands.py", "isovar/__main__.py"):
        head = (REPO_ROOT / rel).read_bytes()[:1024]
        assert b"PYTHON_ARGCOMPLETE_OK" in head, f"Missing marker in {rel}"


def test_top_level_import_and_parser_avoid_heavy_modules():
    forbidden = HEAVY_RUNTIME_MODULES + ("varcode",)
    code = (
        "import sys\n"
        "import argcomplete\n"
        "import isovar\n"
        "import isovar.cli\n"
        "import isovar.cli.commands as cmds\n"
        "cmds.build_parser()\n"
        "cmds.build_completion_parser({})\n"
        f"loaded = [m for m in {forbidden!r} if m in sys.modules]\n"
        "assert not loaded, f'Heavy modules loaded at top-level CLI import: {loaded}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("subcommand", sorted(commands.COMMANDS))
def test_subcommand_parser_loading_avoids_heavy_runtime_modules(subcommand):
    forbidden = list(HEAVY_RUNTIME_MODULES)
    if subcommand in NON_VARCODE_SUBCOMMANDS:
        forbidden.append("varcode")
    code = (
        "import sys\n"
        "from isovar.cli import commands\n"
        f"parser = commands.load_subcommand_parser({subcommand!r})\n"
        "assert parser is not None\n"
        f"loaded = [m for m in {tuple(forbidden)!r} if m in sys.modules]\n"
        f"assert not loaded, f'Heavy modules loaded for {subcommand}: {{loaded}}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "comp_line,expected_items",
    [
        ("isovar ", set(commands.COMMANDS)),
        ("isovar sv-", {"sv-source", "sv-rna", "sv-audit"}),
        ("isovar --", {"--help", "--version"}),
        ("isovar run --", {"--vcf", "--bam", "--output", "--min-alt-rna-reads", "--log-level", "--quiet"}),
        ("isovar fusion-rna --", {"--input", "--output", "--peptide-lengths", "--min-fragments", "--output-dir", "--dpi"}),
        ("isovar sv-rna --", {"--input", "--output", "--bam", "--min-anchor-bases", "--dense-support", "--scratch-dir"}),
        ("isovar sv-audit --", {"--manifest", "--output-dir", "--resume", "--event-id", "--library-cache"}),
        ("isovar allele-counts --", {"--vcf", "--bam", "--cell-umi-labels", "--min-base-quality"}),
        ("isovar cohort-report --", {"--manifest", "--output-dir", "--resume", "--render-only"}),
        ("isovar isoforms ", {"reconstruct", "quantify"}),
        ("isovar isoforms reconstruct --", {"--bam", "--gene", "--gene-set", "--gtf", "--output-dir"}),
    ],
)
def test_argcomplete_returns_expected_completions(tmp_path, comp_line, expected_items):
    out_file = tmp_path / "argcomplete_out.txt"
    env = _subprocess_env(
        _ARGCOMPLETE="1",
        _ARGCOMPLETE_IFS="\n",
        _ARGCOMPLETE_STDOUT_FILENAME=str(out_file),
        COMP_LINE=comp_line,
        COMP_POINT=str(len(comp_line)),
    )
    result = subprocess.run(
        [sys.executable, "-c", "from isovar.cli.commands import run; run()"],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    completions = {
        line.strip()
        for line in out_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    assert expected_items <= completions, (
        f"Expected {expected_items} in completions {completions}"
    )


def test_completion_helpers():
    assert is_completing({"_ARGCOMPLETE": "1"})
    assert not is_completing({})
    assert completing_subcommand(commands.COMMANDS, {}) is None
    assert (
        completing_subcommand(
            commands.COMMANDS,
            {"_ARGCOMPLETE": "1", "COMP_LINE": "isovar fusion-rna --in", "COMP_POINT": "22"},
        )
        == "fusion-rna"
    )
    assert (
        completing_subcommand(
            commands.COMMANDS,
            {"_ARGCOMPLETE": "1", "COMP_LINE": "isovar fus", "COMP_POINT": "10"},
        )
        is None
    )


def test_all_subcommand_parsers_use_help_argument_parser():
    for name in sorted(commands.COMMANDS):
        parser = commands.load_subcommand_parser(name)
        assert isinstance(parser, HelpArgumentParser), (
            f"Subcommand {name!r} parser is {type(parser).__name__}, expected HelpArgumentParser"
        )


def test_help_styling_plain_in_format_help_and_cyan_only_on_tty(monkeypatch):
    class _TTYBuffer(io.StringIO):
        def isatty(self):
            return True

    parser = HelpArgumentParser(prog="isovar-test", description="Test CLI")
    parser.add_argument("--vcf", help="Input ``VCF`` path")
    plain = parser.format_help()
    assert "\033[" not in plain
    assert "``" not in plain

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    non_tty = io.StringIO()
    assert not color_enabled(non_tty)
    parser.print_help(file=non_tty)
    assert "\033[" not in non_tty.getvalue()

    tty = _TTYBuffer()
    assert color_enabled(tty)
    parser.print_help(file=tty)
    styled = tty.getvalue()
    assert "\033[1musage:\033[0m" in styled
    assert "\033[36m--vcf\033[0m" in styled
    assert "\033[1;36m" in styled
    assert style_table_header("contig  start\n7       100\n", tty).startswith("\033[1;36m")

    monkeypatch.setenv("NO_COLOR", "1")
    assert not color_enabled(tty)
    assert "\033[" not in style_help(plain, tty)

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert not color_enabled(tty)
