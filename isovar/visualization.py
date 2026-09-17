"""Auditable, static mutation-evidence figures; Matplotlib is an optional extra.

The data builder uses the existing protein pipeline without changing its
defaults. Genomic, cDNA and protein coordinates are deliberately separate.
"""

from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import inspect
import json
from pathlib import Path
import re

from .default_parameters import (
    PLOT_COMPARE_ASSEMBLY, PLOT_DPI, PLOT_MAX_ROWS, PLOT_OVERVIEW_DPI, PLOT_VIEW, PLOT_VIEWS, PLOT_WIDTH,
    PROTEIN_SEQUENCE_LENGTH,
)
from .dna import reverse_complement_dna
from .effect_prediction import predicted_effects_for_variant
from .protein_sequence_creator import ProteinSequenceCreator
from .protein_sequence_helpers import mutant_peptide_window_count
from .read_identity import fragment_ids, observation_groups
from .variant_helpers import base0_interval_for_variant, interbase_range_affected_by_variant_on_transcript


def _creator_settings(creator):
    return {name: getattr(creator, name)
            for name in inspect.signature(ProteinSequenceCreator).parameters}


def _transcript_data(transcript, rna_supported=False):
    return dict(id=transcript.id, name=getattr(transcript, "name", None), gene=transcript.gene_name,
                strand=transcript.strand, rna_supported=rna_supported,
                exons=sorted((int(a) - 1, int(b)) for a, b in transcript.exon_intervals))


def _reference_predictions(variant, whitelist, length, transcripts):
    """Reference-plus-one-edit predictions, never used to fill missing RNA.

    Anchor the display to the genomic edit's codon, as Isovar does, rather
    than Varcode's repeat-normalized protein event. Keep every transcript's
    annotation and unavailability reason, even when local sequences coincide.
    """
    predictions = []
    for effect in predicted_effects_for_variant(variant, transcript_id_whitelist=whitelist):
        transcript = effect.transcript
        if transcript is None:
            continue
        row = dict(transcript_id=transcript.id, effect=type(effect).__name__,
                   description=effect.short_description, protein=None)
        predictions.append(row)
        transcripts.setdefault(transcript.id, _transcript_data(transcript))
        sequence = effect.mutant_protein_sequence
        if sequence is None or not transcript.complete:
            row["unavailable_reason"] = "No complete protein prediction"
            continue
        try:
            start, end = interbase_range_affected_by_variant_on_transcript(variant, transcript)
            cds = min(transcript.start_codon_spliced_offsets)
            if start < cds or end - start != len(variant.ref):
                raise ValueError("Edit is not contiguous within the coding transcript")
        except ValueError as error:
            row["unavailable_reason"] = str(error)
            continue
        anchor = (start - cds) // 3
        left = max(0, anchor - length // 2)
        right = min(len(sequence), left + length)
        frameshift = (len(variant.ref) - len(variant.alt)) % 3 != 0
        altered_end = right if frameshift else (start - cds + len(variant.alt) + 2) // 3
        row["protein"] = dict(amino_acids=sequence[left:right], mutation_start=anchor - left,
                              mutation_end=min(right, altered_end) - left,
                              protein_start_1based=left + 1, anchor_1based=anchor + 1,
                              frameshift=frameshift, ends_with_stop_codon=False)
    return sorted(predictions, key=lambda p: (p["transcript_id"], p["description"]))


def _witness_data(translation):
    """Keep one real assembly, in transcript orientation, with all its support."""
    sequence = translation.untrimmed_variant_sequence
    reverse = translation.reference_context.strand == "-"
    prefix = sequence.suffix if reverse else sequence.prefix
    suffix = sequence.prefix if reverse else sequence.suffix
    left, right = -len(prefix), len(sequence.alt) + len(suffix)
    spans = Counter()
    genomic_blocks = set()
    junctions = Counter()
    groups = observation_groups(sequence.reads)
    for group in groups:
        before = max(len(read.suffix if reverse else read.prefix) for read in group)
        after = max(len(read.prefix if reverse else read.suffix) for read in group)
        spans[(max(left, -before), min(right, len(sequence.alt) + after))] += 1
        observed_junctions = set()
        for read in group:
            # Clip mappings to the sequence entering this witness. Blocks
            # remain in forward genomic coordinates, even on minus strands.
            query_start = max(0, len(read.prefix) - len(sequence.prefix))
            query_end = min(len(read.sequence), len(read.prefix) + len(read.allele) + len(sequence.suffix))
            blocks = []
            for q0, q1, r0, _ in read.reference_blocks:
                start, end = max(q0, query_start), min(q1, query_end)
                if start < end:
                    blocks.append((r0 + start - q0, r0 + end - q0))
            genomic_blocks.update(blocks)
            gaps = {(a[1], b[0]) for a, b in zip(blocks, blocks[1:])}
            observed_junctions.update(set(read.splice_junctions).intersection(gaps))
        junctions.update(observed_junctions)
    coverage = sequence.coverage()
    if reverse:
        coverage = coverage[::-1]
    return dict(
        strand=translation.reference_context.strand,
        cdna=(reverse_complement_dna(sequence.sequence) if reverse else sequence.sequence),
        start=left, end=right, alt_length=len(sequence.alt),
        coverage=[int(x) for x in coverage],
        spans=[dict(start=a, end=b, observations=count) for (a, b), count in sorted(spans.items())],
        observations=len(groups), templates=len(fragment_ids(sequence.reads)),
        spanning_observations=spans.get((left, right), 0),
        genomic_blocks=sorted(genomic_blocks),
        junctions=[dict(start=a, end=b, observations=n) for (a, b), n in sorted(junctions.items())],
        transcript_ids=sorted(t.id for t in translation.reference_context.transcripts),
        orf=dict(cdna=translation.variant_orf.cdna_sequence,
                 codon_offset=translation.variant_orf.offset_to_first_complete_codon,
                 variant_start=translation.variant_orf.variant_cdna_interval_start,
                 variant_end=translation.variant_orf.variant_cdna_interval_end),
    )


def collect_visualization_data(
        variant, read_evidence, creator_kwargs=None,
        compare_assembly=PLOT_COMPARE_ASSEMBLY, transcript_id_whitelist=None,
        include_reference_predictions=True):
    """Collect JSON-serializable evidence without sampling analysis inputs.

    Parameters
    ----------
    variant : varcode.Variant
    read_evidence : isovar.ReadEvidence
        Already collected observations, shared unchanged between both modes.
    creator_kwargs : dict, optional
        Arguments for ProteinSequenceCreator; omitted values use API defaults.
    compare_assembly : bool
        Compare assembly on/off, overriding only that setting in each mode.
    transcript_id_whitelist : set of str, optional
    include_reference_predictions : bool
        Add Varcode predictions on the same annotation, independently of RNA.

    Returns
    -------
    dict
        Returned protein alternatives, one supporting cDNA witness per protein,
        settings, local transcript models, and counts. ``protein`` retains the
        top result for existing consumers; ``proteins`` retains all results
        allowed by the creator's cap. Set ``max_protein_sequences_per_variant``
        to None in creator_kwargs to inspect every recovered alternative.
        No read names are exported.
    """
    from . import __version__
    from varcode import __version__ as varcode_version

    creator_kwargs = dict(creator_kwargs or {})
    requested = ProteinSequenceCreator(**creator_kwargs)
    modes = [True, False] if compare_assembly else [requested.variant_sequence_assembly]
    results, transcripts = [], {}
    for assembly in modes:
        creator = ProteinSequenceCreator(**dict(creator_kwargs, variant_sequence_assembly=assembly))
        proteins = creator.sorted_protein_sequences_for_variant(
            variant, read_evidence, transcript_id_whitelist=transcript_id_whitelist)
        entry = dict(assembly=assembly, settings=_creator_settings(creator), protein=None, proteins=[])
        selected_ids = {t.id for p in proteins[:1] for translation in p.translations
                        for t in translation.reference_context.transcripts}
        for protein in proteins:
            # A protein may aggregate synonymous RNA sequences or several
            # frames. Display one witness, never their synthetic union.
            witness = min(protein.translations, key=lambda t: (
                -len(fragment_ids(t.untrimmed_variant_sequence.reads)),
                t.untrimmed_variant_sequence.sequence,
                len(t.untrimmed_variant_sequence.prefix),
                tuple(sorted(x.id for x in t.reference_context.transcripts))))
            ids = set()
            for translation in protein.translations:
                for transcript in translation.reference_context.transcripts:
                    ids.add(transcript.id)
                    supported = transcript.id in selected_ids or transcripts.get(transcript.id, {}).get("rna_supported", False)
                    transcripts[transcript.id] = _transcript_data(transcript, rna_supported=supported)
            frames = []
            for translation in protein.translations:
                orf = translation.variant_orf
                frames.append(dict(
                    strand=translation.reference_context.strand,
                    transcript_ids=sorted(t.id for t in translation.reference_context.transcripts),
                    cdna=orf.cdna_sequence,
                    codon_offset=orf.offset_to_first_complete_codon,
                    variant_start=orf.variant_cdna_interval_start,
                    variant_end=orf.variant_cdna_interval_end,
                    variant_codon_phase=(orf.variant_cdna_interval_start - orf.offset_to_first_complete_codon) % 3))
            entry["proteins"].append(dict(
                amino_acids=protein.amino_acids,
                mutation_start=protein.mutation_start_idx, mutation_end=protein.mutation_end_idx,
                frameshift=protein.frameshift, ends_with_stop_codon=protein.ends_with_stop_codon,
                contains_mutation=protein.contains_mutation,
                templates=protein.num_supporting_fragments,
                peptide_windows=mutant_peptide_window_count(protein, creator.protein_context_peptide_length),
                transcript_ids=sorted(ids), translation_count=len(protein.translations),
                frames=frames, witness=_witness_data(witness)))
        entry["protein"] = next(iter(entry["proteins"]), None)
        results.append(entry)
    predictions = (_reference_predictions(
        variant, transcript_id_whitelist, requested.protein_sequence_length or PROTEIN_SEQUENCE_LENGTH, transcripts)
        if include_reference_predictions else [])
    start, end = base0_interval_for_variant(variant)
    counts = {name: dict(observations=len(observation_groups(getattr(read_evidence, name + "_reads"))),
                        templates=len(fragment_ids(getattr(read_evidence, name + "_reads"))))
              for name in ("ref", "alt", "other")}
    return dict(
        schema_version=1, isovar_version=__version__,
        result_filters_applied=False,
        variant=dict(contig=variant.contig, start=variant.start, ref=variant.ref, alt=variant.alt,
                     reference=variant.reference_name, interval=[start, end]),
        genes=sorted({t["gene"] for t in transcripts.values()}),
        counts=counts, modes=results, reference_predictions=predictions,
        varcode_version=varcode_version,
        transcript_whitelist=(None if transcript_id_whitelist is None else sorted(transcript_id_whitelist)),
        transcripts=[transcripts[k] for k in sorted(transcripts)],
        limitations=[
            "Protein candidates are shown before run_isovar result-level filters.",
            "Protein selection is not independent biological validation or a clinical recommendation.",
            "Detail panels show the top protein; proteins retains every result allowed by the recorded creator cap.",
            "Each protein has one selected cDNA witness and all contributing frame contexts; these are not all biological ORFs.",
            "Varcode tracks assume the reference transcript plus only the nominated edit; they are not RNA evidence.",
            "Spans and coverage count deduplicated post-merge observations, not independent molecules.",
            "Distinct read names do not establish molecular independence; no coordinate/UMI deduplication is inferred.",
        ])


def timestamped_run_directory(output_dir):
    """Create a fresh UTC run directory; never reuse/overwrite an earlier run."""
    directory = Path(output_dir).expanduser() / datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S-%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def variant_directory_name(data):
    variant = data["variant"]
    label = "-".join(data["genes"]) or "variant"
    slug = re.sub(r"[^A-Za-z0-9_.-]", "_", "%s-%s-%s-%s-%s" % (
        label, variant["contig"], variant["start"], variant["ref"] or "empty", variant["alt"] or "del"))
    return slug if len(slug) <= 120 else slug[:100] + "-" + sha256(slug.encode()).hexdigest()[:12]


def _plot_imports():
    try:
        from matplotlib.figure import Figure
        from matplotlib import rc_context
        from matplotlib.patches import Rectangle
    except ImportError as error:
        raise ImportError("Figures require the optional plotting extra: pip install 'isovar[plot]'") from error
    return Figure, rc_context, Rectangle


# High-contrast, colorblind-friendly marks on an opaque white background.
BLUE, ORANGE, GRAY, INK = "#0072B2", "#D55E00", "#707070", "#202020"


def _style_axis(ax, title):
    ax.set_facecolor("white")
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=16)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#B0B0B0")
    ax.tick_params(axis="both", labelsize=11, length=3, colors=INK)
    ax.set_yticks([])


def _side_note(ax, text, y=1, transform=None):
    """Reserve the right margin for metadata, outside the plotting area."""
    return ax.text(1.035, y, text, transform=transform or ax.transAxes,
                   ha="left", va="top", fontsize=10, color=GRAY, linespacing=1.6)


def _prediction_groups(data):
    groups = {}
    for prediction in data.get("reference_predictions", []):
        p = prediction["protein"]
        if p is not None:
            key = (p["amino_acids"], p["mutation_start"], p["mutation_end"])
            groups.setdefault(key, []).append(prediction)
    return list(groups.values())


def _protein_disagreements(proteins):
    """Offsets with different residues/stops in overlapping displayed tracks.

    Offsets are relative to each nominated edit, not a protein homology
    alignment. Missing sequence is lack of coverage, never a disagreement.
    """
    residues = {}
    for protein in proteins:
        if protein is None:
            continue
        sequence = protein["amino_acids"] + ("*" if protein["ends_with_stop_codon"] else "")
        for offset, aa in enumerate(sequence, -protein["mutation_start"]):
            residues.setdefault(offset, set()).add(aa)
    return {offset for offset, observed in residues.items() if len(observed) > 1}


def _draw_protein_rows(ax, rows, rectangle, differences=None):
    """Shared residue geometry for single-product and paginated comparisons."""
    left, right = -1, 1
    if differences is None:
        differences = _protein_disagreements(p for _, p, _, _ in rows)
    for i, (label, p, color, note) in enumerate(rows):
        y = len(rows) - i
        ax.text(-0.025, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11)
        _side_note(ax, note, y=y + .15, transform=ax.get_yaxis_transform())
        if p is None:
            ax.text(0, y, "No translated protein", va="center", fontsize=10, color=GRAY)
            continue
        a, b = p["mutation_start"], p["mutation_end"]
        left, right = min(left, -a - 0.7), max(right, len(p["amino_acids"]) - a + 0.7)
        letters = len(p["amino_acids"]) <= 80
        if not letters:
            ax.plot([-a, len(p["amino_acids"]) - a], [y, y], color=color, lw=4)
            ax.plot([0, b - a], [y, y], color=ORANGE, lw=5)
        end = len(p["amino_acids"]) - a + int(p["ends_with_stop_codon"])
        for offset in sorted(differences):
            if -a <= offset < end:
                ax.plot([offset - .4, offset + .4], [y - .27, y - .27],
                        color="#AA3377", lw=2.5, label="_protein_difference")
        for j, aa in enumerate(p["amino_acids"] if letters else ""):
            changed = a <= j < b
            if changed:
                ax.add_patch(rectangle((j - a - .46, y - .20), .92, .40, color=ORANGE, alpha=.12, lw=0))
            ax.text(j - a, y, aa, ha="center", va="center", fontfamily="DejaVu Sans Mono", fontsize=12,
                    color=ORANGE if changed else color, fontweight="bold" if changed else "normal")
        if a == b:
            ax.plot([-.5, -.5], [y - .23, y + .23], color=ORANGE, lw=2)
        if p["ends_with_stop_codon"]:
            ax.text(len(p["amino_acids"]) - a, y, "*", ha="center", va="center", fontsize=11)
    ax.axvline(-.5, color=ORANGE, alpha=.3, lw=.8, zorder=0)
    ax.set(xlim=(left, right), ylim=(.3, len(rows) + .45),
           xlabel="Amino-acid offset from mutation (N → C)")


def _protein_panel(ax, data, rectangle, max_rows):
    _style_axis(ax, "Protein context")
    rows = []
    for mode in sorted(data["modes"], key=lambda m: m["assembly"]):
        p = mode["protein"]
        note = ("%d aa · %d × %d-mers\n%d templates" % (
            len(p["amino_acids"]), p["peptide_windows"], mode["settings"]["protein_context_peptide_length"],
            p["templates"]) if p else "%d alternate read objects" % data["counts"]["alt"]["observations"])
        if p and len(p["amino_acids"]) > 80:
            note += "\nSequence in evidence.json"
        rows.append(("Isovar\nAssembly " + ("on" if mode["assembly"] else "off"), p,
                     BLUE if mode["assembly"] else GRAY, note))
    groups = _prediction_groups(data)
    for i, group in enumerate(groups[:max_rows]):
        descriptions = sorted({p["description"] for p in group})
        note = "\n".join(descriptions[:2])
        if len(descriptions) > 2:
            note += "\n+%d effect labels in notes" % (len(descriptions) - 2)
        rows.append(("Varcode %d\n(%d transcript%s)" % (i + 1, len(group), "s" if len(group) != 1 else ""),
                     group[0]["protein"], INK, note))
    _draw_protein_rows(ax, rows, rectangle)
    if groups:
        label = "Varcode: reference + nominated edit only"
        if len(groups) > max_rows:
            label += " (%d/%d sequence groups shown)" % (max_rows, len(groups))
        ax.text(0, 1.005, label, transform=ax.transAxes, fontsize=10, color=GRAY)
    ax.text(1.035, 1.005, "Orange: mutation\nMagenta: track difference\nBlank: no sequence", transform=ax.transAxes,
            fontsize=10, color=GRAY, va="bottom")


def _selected_witness(data):
    return next((m["protein"]["witness"] for m in data["modes"] if m["protein"]), None)


def _coverage_panel(ax, data):
    from matplotlib.ticker import MaxNLocator

    _style_axis(ax, "RNA coverage")
    if not any(mode["protein"] for mode in data["modes"]):
        ax.text(.5, .5, "No reconstructed cDNA meeting the current thresholds", transform=ax.transAxes, ha="center")
        ax.set_axis_off()
        return
    for mode in data["modes"]:
        if not mode["protein"]:
            continue
        w = mode["protein"]["witness"]
        ax.step(range(w["start"], w["end"]), w["coverage"], where="mid",
                color=BLUE if mode["assembly"] else GRAY, lw=1.5,
                linestyle="-" if mode["assembly"] else "--",
                label="Assembly on" if mode["assembly"] else "Assembly off")
    floor = data["modes"][0]["settings"]["min_variant_sequence_coverage"]
    ax.axhline(floor, color=GRAY, lw=1, ls=":", label="Coverage floor %d" % floor)
    peak = max([floor, 1] + [max(m["protein"]["witness"]["coverage"], default=0)
                            for m in data["modes"] if m["protein"]])
    ax.set_ylim(0, peak * 1.15)
    ax.set_ylabel("Read objects", fontsize=11)
    ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
    ax.legend(frameon=False, fontsize=10, loc="upper left", bbox_to_anchor=(1.025, 1), borderaxespad=0)
    ax.set_xlabel("cDNA offset from alternate allele (nt; 5′ → 3′)")


def _assembly_panel(ax, data, max_rows):
    _style_axis(ax, "Read overlaps")
    w = _selected_witness(data)
    if w is None:
        ax.text(.5, .5, "No translated assembly to display", transform=ax.transAxes, ha="center")
        ax.set_axis_off()
        return
    spans = w["spans"]
    color = BLUE if next(m for m in data["modes"] if m["protein"])["assembly"] else GRAY
    # Display endpoints evenly across the ordered list when a cap is needed;
    # analysis and coverage still use every observation.
    indices = (range(len(spans)) if len(spans) <= max_rows else
               sorted({round(i * (len(spans) - 1) / (max_rows - 1)) for i in range(max_rows)}))
    selected = [spans[i] for i in indices]
    for i, span in enumerate(selected):
        y = len(selected) - i
        ax.plot([span["start"], span["end"]], [y, y], color=color, lw=3, solid_capstyle="butt")
        if span["observations"] > 1:
            ax.text(span["end"] + 1.5, y, "×%d" % span["observations"], fontsize=10, va="center", color=GRAY)
    ax.plot([w["start"], w["end"]], [0, 0], color=INK, lw=5, solid_capstyle="butt")
    ax.text(-.025, 0, "Reconstructed cDNA", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11)
    ax.axvspan(0, max(w["alt_length"], .3), color=ORANGE, alpha=.16)
    ax.axvline(0, color=ORANGE, lw=1)
    note = "%d templates\n%d spanning read objects" % (w["templates"], w["spanning_observations"])
    if w["observations"] != w["templates"]:
        note += "\n%d total read objects" % w["observations"]
    if len(selected) < len(spans):
        note += "\n%d of %d span groups shown" % (len(selected), len(spans))
    _side_note(ax, note)
    ax.set(xlim=(w["start"] - 2, w["end"] + 9), ylim=(-.7, len(selected) + .6),
           xlabel="cDNA offset from alternate allele (nt; 5′ → 3′)")


def _genomic_projection(data):
    """Piecewise-linear genomic axis, retaining exon lengths and capping gaps."""
    w = _selected_witness(data)
    variant_start, variant_end = data["variant"]["interval"]
    observed = w["genomic_blocks"] if w else []
    lo = min([variant_start] + [a for a, _ in observed]) - 5
    hi = max([variant_end, variant_start + 1] + [b for _, b in observed]) + 5
    exons = sorted((max(a, lo), min(b, hi)) for t in data["transcripts"] for a, b in t["exons"]
                   if a < hi and b > lo)
    merged = []
    for a, b in exons:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
        else:
            merged.append((a, b))
    segments, position, x = [], lo, 0
    for a, b in merged + [(hi, hi)]:
        if a > position:
            width = min(a - position, 25)
            segments.append((position, a, x, x + width))
            x += width
        if b > a:
            segments.append((a, b, x, x + b - a))
            x += b - a
        position = b

    def project(value):
        for a, b, x0, x1 in segments:
            if a <= value <= b:
                return x0 + (value - a) * (x1 - x0) / (b - a)
        raise ValueError("Coordinate outside displayed genomic window")

    return lo, hi, segments, project


def _transcript_panel(ax, data, max_rows, rectangle):
    _style_axis(ax, "Transcript models")
    if not data["transcripts"]:
        ax.text(.5, .5, "No contributing transcript model", transform=ax.transAxes, ha="center")
        ax.set_axis_off()
        return
    lo, hi, segments, project = _genomic_projection(data)
    shown = data["transcripts"][:max_rows]
    for i, transcript in enumerate(shown):
        y = len(shown) - i
        exons = [(max(a, lo), min(b, hi)) for a, b in transcript["exons"] if a < hi and b > lo]
        # GTF exons were converted from closed 1-based to half-open 0-based
        # bounds in collect_visualization_data. Introns connect exon end to
        # next exon start on the forward genomic axis, for either strand.
        for left, right in zip(exons, exons[1:]):
            a, b = left[1], right[0]
            if a < b:
                x0, x1 = project(a), project(b)
                ax.plot([x0, x1], [y, y], color=GRAY, lw=1,
                        label="_annotated_intron")
        for a, b in exons:
            ax.add_patch(rectangle((project(a), y - .17), project(b) - project(a), .34,
                                   color=INK if transcript.get("rna_supported", True) else "#B0B0B0", lw=0, zorder=3))
        arrow = "→" if transcript["strand"] == "+" else "←"
        label = transcript["id"] + " " + arrow
        name = transcript.get("name")
        if name and name != transcript["id"]:
            label += "\n(" + name + ")"
        ax.text(-.025, y, label, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=10, linespacing=1.3)
    w = _selected_witness(data)
    for junction in (w["junctions"] if w else []):
        a, b = junction["start"], junction["end"]
        if lo <= a < b <= hi:
            x0, x1 = project(a), project(b)
            ax.plot([x0, x1], [0, 0], color=INK, lw=1, marker="|", markersize=8,
                    label="_observed_junction")
            ax.text((x0 + x1) / 2, -.22, str(junction["observations"]), ha="center", fontsize=8)
    ax.text(-.025, 0, "RNA junctions", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=10)
    if not w or not w["junctions"]:
        ax.text(.5, 0, "No retained junction evidence", transform=ax.get_yaxis_transform(),
                ha="center", va="center", fontsize=8, color=GRAY)
    a, b = data["variant"]["interval"]
    ax.axvline(project(a), color=ORANGE, lw=1.2, zorder=5)
    if b > a:
        ax.axvspan(project(a), project(b), color=ORANGE, alpha=.15, zorder=4)
    ticks = sorted({lo, hi, a})
    ax.set_xticks([project(p) for p in ticks], [f"{p:,}" for p in ticks])
    for start, end, x0, x1 in segments:
        if end - start > x1 - x0:
            ax.text((x0 + x1) / 2, len(shown) + .55, "//", ha="center", fontsize=10, color=GRAY)
    ax.set(xlim=(project(lo) - 2, project(hi) + 2), ylim=(-.65, len(shown) + .95),
           xlabel="Genomic position (0-based, forward strand)")
    count = "%d model%s" % (len(shown), "s" if len(shown) != 1 else "")
    if len(shown) < len(data["transcripts"]):
        count = "%d of %d models shown" % (len(shown), len(data["transcripts"]))
    note = count + "\nIsoform not established"
    if any(not t.get("rna_supported", True) for t in shown):
        note += "\nBlack: contributes RNA protein\nGray: prediction only"
    if any(b - a > x1 - x0 for a, b, x0, x1 in segments):
        note += "\n// compressed genomic gap"
    _side_note(ax, note)


def plot_variant_evidence(data, view=PLOT_VIEW, max_rows=PLOT_MAX_ROWS):
    """Render an opaque-white Figure from collected evidence (no GUI required).

    Parameters
    ----------
    data : dict
        Output of collect_visualization_data.
    view : {'all', 'protein', 'coverage', 'reads', 'assembly', 'transcripts'}
        'assembly' retains the combined coverage/read view.
    max_rows : int
        Maximum displayed span/model rows, never an analysis-input limit.

    Returns
    -------
    matplotlib.figure.Figure
        Caller may save or further customize the figure.
    """
    if view not in PLOT_VIEWS:
        raise ValueError("Unknown visualization view: %s" % view)
    if isinstance(max_rows, bool) or not isinstance(max_rows, int) or max_rows < 2:
        raise ValueError("max_rows must be an integer >= 2")
    Figure, rc_context, rectangle = _plot_imports()
    witness = _selected_witness(data)
    heights, panels, names = [], [], []
    if view in {"all", "protein"}:
        heights.append(1 + .85 * (len(data["modes"]) + min(max_rows, len(_prediction_groups(data)))))
        names.append("protein")
        panels.append(lambda ax: _protein_panel(ax, data, rectangle, max_rows))
    if view in {"all", "assembly", "coverage"}:
        heights.append(2)
        names.append("coverage")
        panels.append(lambda ax: _coverage_panel(ax, data))
    if view in {"all", "assembly", "reads"}:
        heights.append(2 + .22 * min(max_rows, len(witness["spans"]) if witness else 1))
        names.append("reads")
        panels.append(lambda ax: _assembly_panel(ax, data, max_rows))
    if view in {"all", "transcripts"}:
        heights.append(2 + .48 * min(max_rows, len(data["transcripts"])))
        names.append("transcripts")
        panels.append(lambda ax: _transcript_panel(ax, data, max_rows, rectangle))
    with rc_context({"font.family": "DejaVu Sans", "font.size": 12, "text.color": INK,
                     "axes.labelcolor": INK, "svg.fonttype": "none", "axes.formatter.useoffset": False}):
        height = sum(heights) + 1.7 + .7 * (len(panels) - 1)
        figure = Figure(figsize=(PLOT_WIDTH, height), facecolor="white")
        axes = figure.subplots(len(panels), 1, squeeze=False, gridspec_kw={"height_ratios": heights})[:, 0]
        for name, draw, ax in zip(names, panels, axes):
            ax.set_label(name)
            draw(ax)
        if witness:
            witnesses = [m["protein"]["witness"] for m in data["modes"] if m["protein"]]
            bounds = (min(w["start"] for w in witnesses) - 2, max(w["end"] for w in witnesses) + 9)
            for name, ax in zip(names, axes):
                if name in {"coverage", "reads"}:
                    ax.set_xlim(bounds)
        variant = data["variant"]
        def allele_label(bases):
            return (bases or "-") if len(bases) <= 16 else bases[:10] + "… (%d nt)" % len(bases)

        title = "%s   %s:%s %s>%s" % (
            ", ".join(data["genes"][:3]) or "Mutation evidence", variant["contig"], variant["start"],
            allele_label(variant["ref"]), allele_label(variant["alt"]))
        sample = data.get("provenance", {}).get("sample_label")
        if sample:
            title += " · " + sample
        figure.suptitle(title, x=.035, y=1 - .18 / height, ha="left", fontsize=18, fontweight="bold")
        if view not in {"coverage", "protein", "all"}:
            figure.text(.815, 1 - .7 / height, "Orange: mutation", fontsize=10, color=ORANGE)
        figure.subplots_adjust(left=.17, right=.79, top=1 - 1.05 / height, bottom=.65 / height, hspace=.5)
    return figure


def _figure_caption(data):
    caption = (
        "# Figure notes\n\nReference: " + data["variant"]["reference"] + ".\n\n"
        "Orange marks the mutation or deletion boundary; blue denotes assembly on, gray assembly off. "
        "Protein tracks are ordered no assembly, assembly, then Varcode. Magenta underlines mark "
        "different residues (or a stop versus a residue) at shared mutation-relative offsets, including "
        "positions outside the mutation. This is an edit-anchored comparison, not a homology alignment. "
        "Missing coverage is not a sequence difference. Protein peptide counts are "
        "mutation-overlapping windows of the displayed length.\n\n"
        "Reconstructed cDNA is one actual RNA-derived sequence producing the displayed protein. "
        "Other reconstructions can produce the same protein; their reads are not combined into this track. "
        "The read-overlap panel shows the first mode with a protein, normally assembly on.\n\n"
        "Read objects are post-mate-merge observations, not independent molecules. Spanning read objects "
        "cover the entire displayed reconstructed cDNA; spans are clipped to that sequence. "
        "Identical spans are grouped, with multiplicity marked only when greater than one. "
        "Coverage uses every supporting object, including groups hidden by the display limit. "
        "Protein template counts include all translations contributing that protein; cDNA counts refer "
        "only to the displayed reconstruction. Junction counts are retained CIGAR N observations.\n\n"
        "The title uses normalized 1-based variant coordinates. Genomic tracks use forward-strand "
        "0-based, half-open coordinates, with compressed genomic gaps marked //; cDNA and protein offsets "
        "are transcript-oriented. Straight gray connectors join adjacent annotated exon boundaries; "
        "black connectors on the RNA junction row are observed splice junctions. Transcript names "
        "come from the same annotation as the ENST IDs. Models do not establish a unique isoform.\n\n"
        "These are candidates before result-level filters, not independent validation or a clinical "
        "recommendation. Complete sequences, settings, support and limitations are in evidence.json.\n"
    )
    if data.get("reference_predictions"):
        caption += (
            "\n## Varcode baseline\n\nThese tracks predict the reference transcript plus the nominated variant only, "
            "not observed expression or a patient haplotype. The display is anchored to the genomic edit's codon; "
            "Varcode effect labels may use a different repeat-normalized protein boundary. Equal local predicted "
            "sequences are grouped, not ranked. No predicted residues fill missing Isovar context. Reference "
            "tracks do not mark a terminal stop; their ends may simply be the display boundary. Black transcript "
            "boxes contribute to a displayed RNA protein; gray boxes are prediction-only candidates. Compatibility "
            "can use read alignment outside the displayed retained cDNA; this local track is not a full explanation "
            "of every model exclusion.\n\n"
            "| Track | Transcript (name) | Varcode effect |\n| --- | --- | --- |\n")
        models = {t["id"]: t for t in data["transcripts"]}
        for i, group in enumerate(_prediction_groups(data), 1):
            for p in group:
                t = models[p["transcript_id"]]
                caption += "| Varcode %d | %s (%s) | %s |\n" % (i, t["id"], t.get("name") or "unnamed", p["description"])
        for p in data["reference_predictions"]:
            if p["protein"] is None:
                caption += "| Unavailable | %s | %s: %s |\n" % (
                    p["transcript_id"], p["description"], p["unavailable_reason"])
    return caption


def save_variant_figures(data, output_dir, view=PLOT_VIEW, max_rows=PLOT_MAX_ROWS, dpi=PLOT_DPI):
    """Write PNG/SVG panels, a vector PDF, notes and evidence to a fresh directory.

    'all' exports the overview and four separate panels; 'assembly' exports
    the combined RNA view and its two separate panels. Other views export one.
    all-figures.pdf contains one standalone panel per page, without an overview.
    The overview PNG is a preview capped at PLOT_OVERVIEW_DPI; its SVG is vector.
    """
    if isinstance(dpi, bool) or not isinstance(dpi, int) or dpi < 72:
        raise ValueError("dpi must be an integer >= 72")
    figure = plot_variant_evidence(data, view=view, max_rows=max_rows)
    directory = Path(output_dir) / variant_directory_name(data)
    directory.mkdir(parents=True, exist_ok=False)
    views = (["all", "protein", "coverage", "reads", "transcripts"] if view == "all" else
             ["assembly", "coverage", "reads"] if view == "assembly" else [view])
    _, rc_context, _ = _plot_imports()
    from matplotlib.backends.backend_pdf import PdfPages

    with rc_context({"svg.fonttype": "none", "svg.hashsalt": "isovar", "pdf.fonttype": 42}), \
            PdfPages(directory / "all-figures.pdf", metadata={"CreationDate": None, "ModDate": None}) as pdf:
        for i, panel in enumerate(views):
            if i:
                figure = plot_variant_evidence(data, view=panel, max_rows=max_rows)
            name = "overview" if panel == "all" else panel
            figure.savefig(directory / (name + ".svg"), facecolor="white", transparent=False, metadata={"Date": None})
            render_dpi = min(dpi, PLOT_OVERVIEW_DPI) if panel == "all" else dpi
            figure.savefig(directory / (name + ".png"), facecolor="white", transparent=False, dpi=render_dpi)
            if panel not in {"all", "assembly"}:
                pdf.savefig(figure, facecolor="white", transparent=False)
            figure.clear()
    payload = dict(data, rendering=dict(view=view, panels=views, max_rows=max_rows, dpi=dpi,
                                       overview_dpi=min(dpi, PLOT_OVERVIEW_DPI),
                                       width_inches=PLOT_WIDTH, background="white"))
    (directory / "evidence.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    (directory / "caption.md").write_text(_figure_caption(data))
    return directory
