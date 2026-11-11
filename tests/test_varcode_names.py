"""Isovar uses only varcode names that survive varcode 12 (#691)."""

import ast
from pathlib import Path

# varcode 12 removes these (openvax/varcode#693); isovar pins varcode<12, so
# without this check a reintroduced use would surface only at the upgrade.
REMOVED_KEYWORDS = {"ensembl": "genome="}
REMOVED_ATTRIBUTES = {"ensembl": ".genome"}
REMOVED_FUNCTIONS = {"load_vcf_fast": "load_vcf"}


def test_sources_use_varcode_genome_names():
    root = Path(__file__).parents[1]
    paths = [p for d in ("isovar", "tests", "examples", "benchmarks") for p in (root / d).rglob("*.py")]
    uses = []
    for path in paths:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.keyword) and node.arg in REMOVED_KEYWORDS:
                uses.append((path, node.value.lineno, node.arg + "=", REMOVED_KEYWORDS[node.arg]))
            elif isinstance(node, ast.Attribute) and node.attr in REMOVED_ATTRIBUTES:
                uses.append((path, node.lineno, "." + node.attr, REMOVED_ATTRIBUTES[node.attr]))
            name = getattr(node, "id", None) or getattr(node, "attr", None) or getattr(node, "name", None)
            if name in REMOVED_FUNCTIONS:
                uses.append((path, getattr(node, "lineno", 0), name, REMOVED_FUNCTIONS[name]))
    assert not uses, "\n".join(f"{p.relative_to(root)}:{line}: {old} -> use {new}" for p, line, old, new in uses)
