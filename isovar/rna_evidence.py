"""The RNA support record, and scoped, hashed references to the reads behind it.

Every count of supporting RNA in Isovar's outputs is one record (`rna_support`):
``reads`` (sequenced segments: a mate, a single-end read or a long read),
``fragments`` (templates: both mates count once), and, when cell/UMI labels
were assessed, ``umis`` and ``cells`` with their completeness.

Exports never carry read names. A sequenced segment, identified by its
(read group, QNAME, mate bits), is referred to by the SHA-256 of that identity
together with the evidence scope ``[sample_id, source]``. Two exports with the
same scope give the same ID to the same read, so their support can be combined
without counting a read twice. IDs from different scopes never match, which
says nothing about whether their reads overlap.
"""

from hashlib import sha256
import json

from .read_identity import count_reads, fragment_ids, source_read_ids


def canonical_json(value):
    """Compact, key-sorted, ASCII JSON used for every exported identifier."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def content_identifier(kind, value):
    """``kind`` plus ``_`` and the SHA-256 of the canonical JSON of ``[kind, value]``."""
    return kind + "_" + sha256(canonical_json([kind, value]).encode()).hexdigest()


CELL_UMI_FIELDS = ("umis", "cells", "umis_complete", "cells_complete", "unlabeled_reads",
                   "unknown_library_reads", "label_statuses")

# The support record's columns in every table, in order.
SUPPORT_COLUMNS = ("reads", "fragments", "umis", "cells", "umis_complete", "cells_complete",
                   "unlabeled_reads", "unknown_library_reads")


def support_columns(support, prefix=""):
    """A support record's `SUPPORT_COLUMNS` as table cells, prefixed; null is blank."""
    return {prefix + key: "" if support[key] is None else support[key] for key in SUPPORT_COLUMNS}


def rna_support(identities, label_row=None):
    """
    The RNA support record for a set of reads.

    Parameters
    ----------
    identities : iterable of tuple
        Read identities ``(read group, QNAME, mate bits)``; repeats count once.
    label_row : callable, optional
        A read identity's cell/UMI label row (see `isovar.cell_umi`). Without
        it the cell/UMI fields are null: not assessed.

    Returns
    -------
    dict
        ``reads``, ``fragments``, ``umis``, ``cells``, ``umis_complete``,
        ``cells_complete``, ``unlabeled_reads``, ``unknown_library_reads`` and
        ``label_statuses``.
    """
    from .cell_umi import cell_umi_counts

    identities = sorted({tuple(identity) for identity in identities})
    support = dict(reads=len(identities), fragments=len({identity[:2] for identity in identities}))
    support.update(dict.fromkeys(CELL_UMI_FIELDS) if label_row is None
                   else cell_umi_counts(label_row(identity) for identity in identities))
    return support


def evidence_set(scope, identities):
    """
    Hashed IDs for a set of reads in one evidence scope.

    Parameters
    ----------
    scope : list
        ``[sample_id, source]``, the evidence scope of the export.
    identities : iterable of tuple
        Read identities ``(read group, QNAME, mate bits)``; repeats count once.

    Returns
    -------
    dict
        ``evidence_scope``; ``reads`` and ``fragments`` counts; their hashed
        ``read_ids`` and ``fragment_ids``; and ``evidence_set_id``, one ID
        for the whole set.
    """
    identities = sorted({tuple(identity) for identity in identities})
    fragments = sorted({identity[:2] for identity in identities})
    read_ids = [content_identifier("segment", [scope, identity]) for identity in identities]
    return dict(
        evidence_scope=list(scope), reads=len(identities), fragments=len(fragments), read_ids=read_ids,
        fragment_ids=[content_identifier("fragment", [scope, fragment]) for fragment in fragments],
        evidence_set_id=content_identifier("evidence", sorted(read_ids)))


class EvidenceSets:
    """
    RNA support records for read sets in one evidence scope, storing each
    set's hashed read IDs once, keyed by its evidence set ID.

    Parameters
    ----------
    scope : list
        ``[sample_id, source]``.
    """

    def __init__(self, scope):
        self.scope, self.sets = list(scope), {}

    @staticmethod
    def record(reads, labels=None):
        """
        The RNA support record of ``reads`` (Isovar read objects), and their
        identities, or None when any read lacks one; such reads are counted
        without cell/UMI fields. Cell/UMI fields come from ``labels`` (an
        `isovar.cell_umi.CellUmiEvidence`), or are null.
        """
        reads = list(reads)
        identities = [source_read_ids(read) for read in reads]
        if not all(identities):
            return dict(reads=count_reads(reads), fragments=len(fragment_ids(reads)),
                        **dict.fromkeys(CELL_UMI_FIELDS)), None
        keys = {key for keys in identities for key in keys}
        return (rna_support(keys) if labels is None else labels.support(keys)), keys

    def support(self, reads, labels=None):
        """
        The RNA support record of ``reads`` plus ``evidence_set_id``, whose
        hashed read IDs are stored in ``sets``. The ID is null only when a
        read lacks an identity; no reads is the empty set.
        """
        record, keys = self.record(reads, labels)
        if keys is None:
            return dict(record, evidence_set_id=None)
        evidence = evidence_set(self.scope, keys)
        self.sets[evidence["evidence_set_id"]] = evidence
        return dict(record, evidence_set_id=evidence["evidence_set_id"])


def union_rna_support(supports):
    """Support of several measurements together, counting each read once.

    Use this to combine evidence from different hypotheses, or from different
    exports of the same reads, instead of adding their counts.

    Parameters
    ----------
    supports : iterable of dict
        Read sets, each with ``evidence_scope``, ``read_ids`` and
        ``fragment_ids``: entries of a protein hypothesis export's
        ``evidence_sets``, or an SV ORF export's ``rna_support``.

    Returns
    -------
    dict
        The same fields as `evidence_set` for the union.

    Raises
    ------
    ValueError
        If there are no supports, if they come from different evidence scopes
        (whose overlap cannot be known), or if any lacks read identities.
    """
    supports = list(supports)
    if not supports:
        raise ValueError("No RNA supports to combine")
    if any(support.get("read_ids") is None for support in supports):
        raise ValueError("RNA supports without read identities cannot be combined; pass evidence "
                         "sets, such as an export's evidence_sets entries")
    scopes = {canonical_json(support["evidence_scope"]) for support in supports}
    if len(scopes) > 1:
        raise ValueError("RNA supports from different evidence scopes cannot be combined: "
                         "their overlap is unknown")
    read_ids = sorted({i for support in supports for i in support["read_ids"]})
    fragment_ids = sorted({i for support in supports for i in support["fragment_ids"]})
    return dict(
        evidence_scope=list(supports[0]["evidence_scope"]), reads=len(read_ids),
        fragments=len(fragment_ids), read_ids=read_ids, fragment_ids=fragment_ids,
        evidence_set_id=content_identifier("evidence", read_ids))
