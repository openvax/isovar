"""Scoped, hashed references to the RNA reads behind an exported measurement.

Exports never carry read names. A sequenced segment, identified by its
(read group, QNAME, mate bits), is referred to by the SHA-256 of that identity
together with the evidence scope ``[sample_id, source]``. Two exports with the
same scope give the same ID to the same read, so their support can be combined
without counting a read twice. IDs from different scopes never match, which
says nothing about whether their reads overlap.
"""

from hashlib import sha256
import json


def canonical_json(value):
    """Compact, key-sorted, ASCII JSON used for every exported identifier."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def content_identifier(kind, value):
    """``kind`` plus ``_`` and the SHA-256 of the canonical JSON of ``[kind, value]``."""
    return kind + "_" + sha256(canonical_json([kind, value]).encode()).hexdigest()


def segment_support(scope, identities):
    """Counts and hashed IDs for a set of sequenced segments in one scope.

    Parameters
    ----------
    scope : list
        ``[sample_id, source]``, the evidence scope of the export.
    identities : iterable of tuple
        Segment identities ``(read group, QNAME, mate bits)``; repeats count once.

    Returns
    -------
    dict
        ``evidence_scope``; ``segments`` and ``fragments`` counts; their hashed
        ``segment_ids`` and ``fragment_ids``; and ``evidence_set_id``, one ID
        for the whole set of segments.
    """
    identities = sorted({tuple(identity) for identity in identities})
    fragments = sorted({identity[:2] for identity in identities})
    segment_ids = [content_identifier("segment", [scope, identity]) for identity in identities]
    return dict(
        evidence_scope=list(scope), segments=len(identities), fragments=len(fragments),
        segment_ids=segment_ids,
        fragment_ids=[content_identifier("fragment", [scope, fragment]) for fragment in fragments],
        evidence_set_id=content_identifier("evidence", sorted(segment_ids)))


def union_rna_support(supports):
    """Support of several measurements together, counting each read once.

    Use this to combine evidence from different hypotheses, or from different
    exports of the same reads, instead of adding their counts.

    Parameters
    ----------
    supports : iterable of dict
        Read sets, each with ``evidence_scope``, ``segment_ids`` and
        ``fragment_ids``: entries of a protein hypothesis export's
        ``evidence_sets``, or an SV ORF export's ``rna_support``.

    Returns
    -------
    dict
        The same fields as `segment_support` for the union.

    Raises
    ------
    ValueError
        If there are no supports, if they come from different evidence scopes
        (whose overlap cannot be known), or if any lacks read identities.
    """
    supports = list(supports)
    if not supports:
        raise ValueError("No RNA supports to combine")
    if any(support.get("segment_ids") is None for support in supports):
        raise ValueError("RNA supports without read identities cannot be combined; pass evidence "
                         "sets, such as an export's evidence_sets entries")
    scopes = {canonical_json(support["evidence_scope"]) for support in supports}
    if len(scopes) > 1:
        raise ValueError("RNA supports from different evidence scopes cannot be combined: "
                         "their overlap is unknown")
    segment_ids = sorted({i for support in supports for i in support["segment_ids"]})
    fragment_ids = sorted({i for support in supports for i in support["fragment_ids"]})
    return dict(
        evidence_scope=list(supports[0]["evidence_scope"]), segments=len(segment_ids),
        fragments=len(fragment_ids), segment_ids=segment_ids, fragment_ids=fragment_ids,
        evidence_set_id=content_identifier("evidence", segment_ids))
