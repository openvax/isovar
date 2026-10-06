"""Publish DNA/RNA anatomy and every sequence/support row from a frozen Sid screen.

This annotates existing results; it does not acquire reads or infer translation.
Genomic intervals in the ledger are zero-based, half-open. Markdown loci are
one-based, inclusive; SV cuts and insertion locations are explicitly interbase.
"""

import argparse
from collections import Counter
import csv
from pathlib import Path
import textwrap

from examples.sid_sv_audit.inventory import digest, identity, read_json, write_json


PRODUCTS = {
    "40a3ceff825db98f490f303c67b1fad64e5594e16670e9ed43d7dc1c14be58fa": "Bulk SARC0277 T2",
    "58f81602ec2f87b86a032cf869f1abefc9f4235dbb1ce082e0b6bbc085b155de": "ONT tagged T2",
    "6e9147690525529f9591602afd1191f9862f4df16117a865c6726654456c9368": "10x Kamil T2",
    "8275c5f211640cb978ef5f54e7cf22e71dd96efc8a8f84fb14039c097182989b": "PacBio T1",
}


def literal_edit(allele):
    """Remove VCF padding without changing the original allele or HGVS placement."""
    contig, position, ref, alt = allele
    if not ref or not alt or (set(ref + alt) - set("ACGT")):
        return None
    start = position - 1
    while ref and alt and ref[0] == alt[0]:
        ref, alt, start = ref[1:], alt[1:], start + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return dict(contig=contig, start=start, end=start + len(ref), ref=ref, alt=alt,
                length_delta=len(alt) - len(ref), placement="VCF_padding_removed_not_HGVS_normalized")


def rna_edit_description(event):
    """Describe the local strand-oriented edit without inventing altered splicing."""
    edit = event["literal_edit"]
    if edit is None:
        return "The nominated duplication has no resolved literal allele; its altered RNA sequence and protein are unknown."
    strands = sorted({a["strand"] for a in event["anatomy"] if a["location"] == "exon"
                      and not any(e["at_exon_boundary"] for e in a["exons"])})
    changes = []
    for strand in strands:
        rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        ref, alt = (edit["ref"], edit["alt"]) if strand == "+" else (rc(edit["ref"]), rc(edit["alt"]))
        changes.append("%s-strand local allele, if the affected exon is retained: `%s` → `%s`" % (
            "Plus" if strand == "+" else "Minus", ref or "∅", alt or "∅"))
    observed = any(o.get("proteins", 0) for o in event["outcomes"])
    description = "; ".join(changes) + ". " if changes else (
        "The event has no wholly exonic, non-boundary placement in the supplied models; the mature RNA edit is unknown. ")
    description += ("Coding windows below establish local RNA sequence in the stated products; they do not establish the full mature transcript. "
                    if observed else "No resolved coding window establishes an altered protein in these inputs. ")
    description += "An intronic or exon-boundary event can alter splicing; its mature-RNA outcome remains unresolved without an observed altered splice path. "
    if edit["length_delta"]:
        description += "The genomic length change applies to a coding frame only where that literal edit is retained in CDS; intronic bases are not automatically translated."
    return description.rstrip()


def anatomy(model, start, end):
    """Locate an interval or interbase insertion in transcript order, on either strand."""
    exons = sorted(model["exons"], reverse=model["strand"] == "-")
    offset, hits = 0, []
    cds_start, cds_end = model.get("cds_start"), model.get("cds_end")
    for number, (a, b) in enumerate(exons, 1):
        # Insertions at exon edges are boundaries, not asserted mature-RNA insertions.
        lo, hi = max(start, a), min(end, b)
        inserted = start == end and a < start < b
        boundary = start == end and start in (a, b)
        if lo < hi or inserted or boundary:
            q0 = offset + (lo - a if model["strand"] == "+" else b - hi)
            q1 = q0 + max(0, hi - lo)
            if start == end:
                q0 = q1 = offset + (start - a if model["strand"] == "+" else b - start)
            if cds_start is None:
                features = ["exonic; CDS unresolved/noncoding"]
            else:
                features = []
                if q0 < cds_start:
                    features.append("5′ UTR")
                if q1 > cds_start and q0 < cds_end or q0 == q1 and cds_start <= q0 < cds_end:
                    features.append("CDS")
                if q1 > cds_end or q0 == q1 and q0 >= cds_end:
                    features.append("3′ UTR")
            hits.append(dict(exon=number, exon_interval=[a, b], cdna_interval=[q0, q1],
                             features=features, at_exon_boundary=boundary))
        offset += b - a
    introns = [[b, c] for (_, b), (c, _) in zip(sorted(exons), sorted(exons)[1:])
               if max(start, b) < min(end, c) or start == end and b < start < c]
    within_span = min(a for a, _ in exons) <= start < max(b for _, b in exons)
    return dict(transcript_id=model["transcript_id"], gene=model.get("gene"),
                gene_id=model.get("gene_id"), contig=model["contig"], strand=model["strand"],
                exons=hits, overlapping_introns=introns, total_exons=len(exons),
                location="exon_and_intron" if hits and introns else "exon" if hits else
                "intron" if introns or within_span else "outside_transcript",
                cds_start=cds_start, cds_end=cds_end)


def clip_blocks(blocks, interval):
    """Map an ORF interval to genomic sequence, respecting reverse-strand blocks."""
    start, end = interval
    clipped = []
    for block in blocks:
        lo, hi = max(start, block["query_start"]), min(end, block["query_end"])
        if lo >= hi:
            continue
        delta = lo - block["query_start"]
        if block["strand"] == "+":
            a = block["reference_start"] + delta
            b = a + hi - lo
        else:
            b = block["reference_end"] - delta
            a = b - (hi - lo)
        clipped.append(dict(contig=block["contig"], strand=block["strand"],
                            reference_start=a, reference_end=b, query_start=lo, query_end=hi))
    return clipped


def dna_calls(nominations):
    """Retain only explicit DNA VCF evidence; RNA callers do not establish DNA events."""
    calls = {}
    for nomination in nominations:
        original = nomination["original"]
        for call in original.get("original_calls", []):
            vcf = call.get("vcf")
            if not vcf or call.get("source") != "ESVEE_PURPLE_LINX_PASS":
                continue
            key = identity([vcf.get(k) for k in ("source_url", "sample", "chrom", "pos", "vcf_id")])
            calls[key] = call
    return [calls[k] for k in sorted(calls)]


def group_key(row):
    return [row.get("geometry_id") or row["target"][0], row["kind"], row["amino_acids"],
            row.get("nucleotide_sequence"), row.get("mutation_interval"), row["ends_with_stop_codon"]]


def crossed_joins(occurrence):
    """Occurrence junctions retain their original path indices, not list offsets."""
    joins = [j for j in occurrence["junctions"] if j["junction_index"] in occurrence["crossed_junctions"]]
    if {j["junction_index"] for j in joins} != set(occurrence["crossed_junctions"]):
        raise ValueError("Crossed junction missing from occurrence")
    return joins


def transcript_model(genome, tid):
    t = genome.transcript_by_id(tid)
    return dict(transcript_id=tid, gene=t.gene_name or t.gene_id, gene_id=t.gene_id,
                contig="chr" + t.contig, strand=t.strand,
                exons=sorted((a - 1, b) for a, b in t.exon_intervals),
                cds_start=min(t.start_codon_spliced_offsets) if t.complete else None,
                cds_end=max(t.stop_codon_spliced_offsets) + 1 if t.complete else None)


def build(screen_path, audits, checks_path):
    from pyensembl import EnsemblRelease

    screen = read_json(screen_path)
    pins = dict(screen["input_pins"])
    pins[str(screen_path)] = digest(screen_path)
    pins[str(checks_path)] = digest(checks_path)
    # Frozen ledgers already include the annotation and consumed indexes. Do not
    # silently annotate against a newer local database or changed reconstruction.
    for path, expected in pins.items():
        if digest(path) != expected:
            raise ValueError("Changed frozen input: " + path)
    genome = EnsemblRelease(115)
    genome.gene_ids_at_locus("1", 1, 1)  # Initialize the local annotation paths.
    for path in [genome.gtf_path, *genome.transcript_fasta_paths]:
        if path not in pins:
            raise ValueError("Active annotation source was not pinned by the screen: " + path)
    events, source_metadata, details, result_cache = {}, {}, {}, {}
    model_cache = {}

    def model(tid):
        if tid not in model_cache:
            model_cache[tid] = transcript_model(genome, tid)
        return model_cache[tid]

    for directory in audits:
        inventory_path, refs_path = directory / "inventory.json.gz", directory / "references.json.gz"
        inventory, refs = read_json(inventory_path), read_json(refs_path)
        for path in (inventory_path, refs_path):
            pins[str(path)] = digest(path)
        if refs["inventory_sha256"] != pins[str(inventory_path)]:
            raise ValueError("References belong to another inventory")
        for path in sorted((directory / "results").glob("*/*.json.gz")):
            result = read_json(path)
            if pins.get(str(path.resolve())) != digest(path):
                raise ValueError("Unpinned SV result: " + str(path))
            request = result["request"]
            gid, sid = request["geometry_id"], request["source_id"]
            source_metadata[sid] = inventory["sources"][sid]
            geometry = inventory["geometries"][gid]
            if gid not in events:
                nominations = [dict(id=n, **inventory["nominations"][n]) for n in geometry["nominations"]]
                calls = dna_calls(nominations)
                locations = []
                for side, endpoint in enumerate(geometry["breakends"], 1):
                    pos = endpoint["position"] - (endpoint["retained_side"] == "left")
                    for tid in refs["assignments"][gid]["transcripts"]:
                        m = refs["models"][tid]
                        if m["contig"] == endpoint["contig"]:
                            a = anatomy(dict(m, gene=model(tid)["gene"], gene_id=model(tid)["gene_id"]), pos, pos + 1)
                            if a["location"] != "outside_transcript":
                                locations.append(dict(side=side, retained_base_one_based=pos + 1, **a))
                events[gid] = dict(id=gid, kind="SV", names=sorted({n["original"].get("name", n["id"]) for n in nominations}),
                    breakends=geometry["breakends"], nominations=nominations, original_dna_calls=calls,
                    dna_origin="DNA_VCF_linked" if calls else "unknown_RNA_only_nomination",
                    anatomy=locations, excluded_transcripts=refs["assignments"][gid]["excluded_transcripts"],
                    unsupported_cds={t: refs["unsupported_cds"][t] for t in refs["assignments"][gid]["transcripts"]
                                     if t in refs["unsupported_cds"]}, outcomes=[])
            events[gid]["outcomes"].append(dict(source_id=sid, orientation=request["orientation"],
                status=result["status"], candidates=len(result["candidates"]) if "candidates" in result else None,
                records=result.get("records"), paths=result.get("paths"), event_paths=result.get("event_paths"),
                acquisition_status=result["acquisition_status"], input_limitations=result["input_limitations"],
                discovery=result.get("discovery"), support_acquisition=result.get("support_acquisition"),
                reason=result.get("reason"), error=result.get("error")))
            result_cache[str(path.resolve())] = {c["candidate_id"]: c for c in result.get("candidates", [])}
            if result.get("candidates"):
                detail_path = directory / result["reconstruction"]["path"]
                expected = result["reconstruction"]["sha256"]
                if digest(detail_path) != expected:
                    raise ValueError("Changed reconstruction")
                pins[str(detail_path)] = expected
                details[str(path.resolve())] = {p["path_id"]: p for p in read_json(detail_path)["paths"]}
        selection_path = directory / "priority-selection.json"
        if not selection_path.exists():
            continue
        pins[str(selection_path)] = digest(selection_path)
        selection = read_json(selection_path)
        for tid, entry in selection["small_variants"].items():
            variant = entry["variant"]
            annotations = variant["annotations"]
            allele = variant["alleles"][0]
            if (len(allele[2]) == len(allele[3]) and "splice_" not in annotations["consequence"]
                    and variant["gene"] != "MUC3A"):
                continue
            edit = literal_edit(allele)
            # The nonliteral duplication has only a nominated locus, not a span.
            a, b = (edit["start"], edit["end"]) if edit else (allele[1] - 1, allele[1])
            gene_ids = genome.gene_ids_at_locus(allele[0].removeprefix("chr"), a + 1, max(a + 1, b))
            tids = sorted({t for g in gene_ids for t in genome.transcript_ids_of_gene_id(g)})
            record = annotations.get("source_record") or {}
            dna = [v for v in record.get("vafs", []) if v["assay"] in ("WGS", "WES")]
            events[tid] = dict(id=tid, kind="small_variant",
                names=[variant["gene"] or "Unassigned %s:%s" % (allele[0], allele[1])],
                original_allele=allele, literal_edit=edit, variant_metadata=variant,
                original_dna_support=dna,
                dna_origin="catalogued_allele; original_DNA_VAF_metadata" if dna else
                           "catalogued_allele; original_DNA_support_unavailable",
                anatomy=[anatomy(model(t), a, b) for t in tids], outcomes=[])
            for sid in selection["source_ids"]:
                source_metadata[sid] = inventory["sources"][sid]
            if edit is None:
                events[tid]["outcomes"] = [dict(source_id=sid, status="nonliteral_unassessable", proteins=None,
                                                allele_support=None) for sid in selection["source_ids"]]
    for outcome in screen["small_variant_outcomes"]:
        events[outcome["target_id"]]["outcomes"].append(outcome)
    rows, groups = [], {}
    for frozen in screen["candidates"]:
        row = dict(frozen)
        eid = row.get("geometry_id") or row["target"][0]
        row["event_id"] = eid
        key = "seq-" + identity(group_key(row))[:20]
        row["report_sequence_id"] = key
        if row["kind"] == "SV_RNA_ORF":
            candidate = result_cache[row["result_path"]][row["candidate_id"]]
            paths = details[row["result_path"]]
            occurrences = []
            for occurrence in candidate["occurrences"]:
                path = paths[occurrence["path_id"]]
                q0, q1 = occurrence["query_interval"]
                if path["sequence"][q0:q1] != row["nucleotide_sequence"]:
                    raise ValueError("ORF nucleotide sequence differs from its RNA path")
                occurrences.append(dict(occurrence,
                    genomic_blocks=clip_blocks(path["blocks"], [q0, q1]),
                    rna_path_blocks=path["blocks"], rna_path_length=len(path["sequence"]),
                    rna_path_end_reasons=path["end_reasons"],
                    rna_path_unplaced_intervals=path["unplaced_intervals"]))
            row["rna_occurrences"] = occurrences
            row["support_details"] = {k: v for k, v in candidate["rna_support"].items()
                                      if k not in ("fragment_ids", "read_ids", "cell_barcodes")}
        rows.append(row)
        groups.setdefault(key, dict(id=key, event_id=eid, kind=row["kind"], amino_acids=row["amino_acids"],
            nucleotide_sequence=row.get("nucleotide_sequence"), mutation_interval=row.get("mutation_interval"),
            ends_with_stop_codon=row["ends_with_stop_codon"], row_indices=[]))["row_indices"].append(len(rows) - 1)
    return dict(scope=screen["scope"], annotation="Ensembl 115 / GRCh38",
        coordinates="genomic_zero_based_half_open; cDNA_zero_based_half_open; VCF_original_one_based",
        events=events, candidates=rows, sequences=groups, sources=source_metadata,
        counts=dict(events=len(events), source_hypotheses=len(rows), sequence_groups=len(groups),
                    kinds=dict(Counter(r["kind"] for r in rows))), input_pins=pins,
        native_window_checks=read_json(checks_path), proteome_comparison=screen["comparison"],
        code_sha256=digest(__file__), translation_observed=False, tumor_specificity_assessed=False)


def table(headers, rows):
    def cell(value):
        return str(value if value is not None else "unknown").replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(map(cell, headers)) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |"] +
                     ["| " + " | ".join(map(cell, row)) + " |" for row in rows]) + "\n"


def anatomy_text(a):
    bits = ["exon %s/%s %s; cDNA [%s,%s)%s" % (
        e["exon"], a["total_exons"], "/".join(e["features"]), *e["cdna_interval"],
        "; exon boundary" if e["at_exon_boundary"] else "") for e in a["exons"]]
    bits += ["intron %s:%s–%s" % (a["contig"], x + 1, y) for x, y in a["overlapping_introns"]]
    return "; ".join(bits) or a["location"]


def locus(block):
    return "%s:%s–%s(%s)" % (block["contig"], block["reference_start"] + 1,
                               block["reference_end"], block["strand"])


def junction_text(j):
    left, right = j["left"], j["right"]
    return "%s:%s(%s) → %s:%s(%s); unplaced `%s`; %s" % (
        left[0], left[1] + 1, left[2], right[0], right[1] + 1, right[2],
        j["unplaced_bases"] or "∅", j["relation"])


def start_text(evidence):
    regions = sorted({a["region"] for a in evidence["assessments"] if a.get("region")})
    return "%s; tier %s; %s; initiation unobserved" % (
        evidence["status"], evidence["tier"] or "unresolved", "/".join(regions) or "no reference anchor")


def event_markdown(event, ledger):
    rows, sequences = ledger["candidates"], ledger["sequences"]
    eid = event["id"]
    seqs = [s for s in sequences.values() if s["event_id"] == eid]
    seqs.sort(key=lambda s: (-max(rows[i]["fragments"] for i in s["row_indices"]), s["id"]))
    out = ["# " + " / ".join(event["names"]), "", "Event ID: `" + eid + "`.", "",
           "[Report and counting definitions](../../../sid-neoorf-event-report.md). "
           "[Complete JSON ledger](../events.json.gz) retains all placements and provenance.", "", "## Original DNA event", "",
           "Origin status: `" + event["dna_origin"] + "`.", ""]
    if event["kind"] == "SV":
        out += [table(["Side", "Nominated interbase cut", "Retained side", "Adjacent retained base (one-based)"],
                      [(i, "%s:%s" % (b["contig"], b["position"]), b["retained_side"],
                        b["position"] + (b["retained_side"] == "right"))
                       for i, b in enumerate(event["breakends"], 1)])]
        calls = event["original_dna_calls"]
        if calls:
            out += [table(["DNA sample / caller", "Original VCF allele", "VF / SF / DF", "PURPLE AF / JCN", "Call / cluster"], [
                (c["vcf"]["sample"] + " / " + c["vcf"]["caller"],
                 "%s:%s `%s` → `%s`" % tuple(c["vcf"][k] for k in ("chrom", "pos", "ref", "alt")),
                 " / ".join(str(c["vcf"]["info"].get(k, "unavailable")) for k in ("VF", "SF", "DF")),
                 "%s / %s" % (c["vcf"]["info"].get("PURPLE_AF"), c["vcf"]["info"].get("PURPLE_JCN")),
                 "%s / %s / %s" % (c["vcf"]["vcf_id"], c["row"].get("svtype"), c["row"].get("cluster_resolved_type")))
                for c in calls])]
            out += ["Original VCF sources:", ""] + ["- [%s](%s)" % (c["vcf"]["sample"], c["vcf"]["source_url"]) for c in calls] + [""]
        else:
            out += ["No original DNA VCF call is linked to this selected nomination. The RNA join alone does not identify a somatic DNA rearrangement.", ""]
        out += ["Nomination aliases: " + ", ".join("`%s`" % n["id"] for n in event["nominations"]) + ".", ""]
    else:
        chrom, pos, ref, alt = event["original_allele"]
        annotations = event["variant_metadata"]["annotations"]
        record = annotations.get("source_record") or {}
        out += ["Original GRCh38 VCF-style allele: **%s:%s** `%s` → `%s`." % (chrom, pos, ref, alt), "",
                "Source consequence: `%s`; source protein label: `%s`; source transcript: `%s`; cDNA label: `%s`." % (
                    annotations.get("consequence") or "unavailable", annotations.get("protein_change") or "unavailable",
                    record.get("refseq_id") or "unavailable", record.get("genomic_change_on_cdna") or "unavailable"), ""]
        edit = event["literal_edit"]
        if edit:
            out += ["After removing VCF padding: genomic **[%s,%s)**; `%s` → `%s`; net length change **%+d nt**. "
                    "This is genomic placement, not HGVS right normalization or a demonstrated mature-RNA edit." % (
                        edit["start"], edit["end"], edit["ref"] or "∅", edit["alt"] or "∅", edit["length_delta"]), ""]
        else:
            out += ["The ALT is nonliteral. Its span/sequence is unresolved; anatomy below describes only the nominated locus. RNA support is unassessable, not zero.", ""]
        if record.get("consistency_warnings"):
            out += ["Source consistency warning: `" + str(record["consistency_warnings"]) + "`.", ""]
        dna = event["original_dna_support"]
        if dna:
            primary = [v for v in dna if "redux" in v["library"] or "tumor_bg" in v["library"]]
            out += ["Original catalogue DNA depth/VAF (selected primary library labels; all processing copies are retained in the ledger). These are source metadata, not a new DNA recount.", "",
                    table(["DNA library", "Assay / timepoint / tissue", "Depth", "VAF"],
                          [(v["library"], "%s / %s / %s" % (v["assay"], v["timepoint"], v["tissue"]),
                            v["depth"], v["vaf"]) for v in primary or dna])]
        else:
            out += ["Original DNA depth/VAF is unavailable in the frozen nomination.", ""]
    out += ["## Transcript anatomy and event location", "",
            "Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. "
            "Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.", "",
            table(["Side", "Gene / transcript", "Strand", "Placement"],
                  [(a.get("side", "allele"), "%s / %s" % (a["gene"], a["transcript_id"]), a["strand"], anatomy_text(a))
                   for a in event["anatomy"]]) if event["anatomy"] else
            "No overlapping transcript model at this nominated event location; gene assignment and coding consequence remain unresolved.", "",
            "## RNA screen outcomes", ""]
    if event["kind"] == "small_variant":
        out += [rna_edit_description(event), ""]
    else:
        missing = {1, 2} - {a["side"] for a in event["anatomy"]}
        for side in sorted(missing):
            out += ["Side %s has no overlapping supplied Ensembl 115 transcript at its adjacent retained base; gene/CDS assignment there is unresolved." % side, ""]
    if event["kind"] == "SV":
        out += [table(["Product", "Orientation", "Status", "Records / paths / event paths", "ORF hypotheses"],
                      [(PRODUCTS[o["source_id"]], o["orientation"], o["status"],
                        " / ".join(str(o[k]) if o[k] is not None else "unknown" for k in ("records", "paths", "event_paths")),
                        o["candidates"]) for o in event["outcomes"]])]
        for outcome in event["outcomes"]:
            if outcome.get("reason") or outcome.get("error"):
                out += ["%s / %s: `%s`; %s. Counts unavailable from this outcome remain unknown." % (
                    PRODUCTS[outcome["source_id"]], outcome["orientation"], outcome["status"],
                    outcome.get("error") or outcome["reason"]), ""]
    else:
        out += [table(["Product", "Status", "RNA ref / alt / other fragments", "Coding windows"],
                      [(PRODUCTS[o["source_id"]], o["status"],
                        " / ".join(str(o["allele_support"][k]["fragments"]) for k in ("ref", "alt", "other"))
                        if o.get("allele_support") else "unassessable", o["proteins"]) for o in event["outcomes"]])]
    if not seqs:
        out += ["No protein sequence is resolved for this event in these recorded outcomes and search bounds. "
                "For splice nominations, alternate exon use, intron retention and exonization have not been comprehensively reconstructed. "
                "Protein sequence and altered transcript structure remain unknown; reference anatomy above is the available description.", ""]
    else:
        out += ["## Every protein sequence / coding window", "",
                "Each sequence is an RNA-supported hypothesis or local coding window. A terminal `*` means an in-window stop, "
                "not a demonstrated full-length protein. SV initiation/frame/polarity may be unresolved. "
                "In-frame windows are controls. Counts below are never summed across products, starts, windows or geometry aliases.", ""]
    for seq in seqs:
        rr = [rows[i] for i in seq["row_indices"]]
        out += ["### " + seq["id"], "", "%s; **%s aa**; stop: **%s**." % (
                    seq["kind"], len(seq["amino_acids"]), seq["ends_with_stop_codon"]), "",
                "```text", "\n".join(textwrap.wrap(seq["amino_acids"] + ("*" if seq["ends_with_stop_codon"] else ""), 80)), "```", ""]
        if seq["mutation_interval"]:
            lo, hi = seq["mutation_interval"]
            out += ["Changed/shifted interval: aa **[%s,%s)**; sequence `%s`. "
                    "This is an offset in this window, not the full reference protein." % (lo, hi, seq["amino_acids"][lo:hi]), ""]
        out += [table(["Product / orientation", "Complete fragments", "Q20 complete", "Contributing allele fragments", "Barcode labels", "Missing qualities"], [
            (PRODUCTS[r["source_id"]] + " / " + r.get("orientation", "coding"), r["fragments"],
             r.get("q20_full_window_fragments", "not assessed"), r.get("contributing_fragments", "see junction/path ledger"),
             r.get("cells", "not counted for these windows"), r.get("missing_quality_reads", r.get("fragments_with_unavailable_qualities"))) for r in rr])]
        if seq["kind"] == "SV_RNA_ORF":
            placements = {}
            for r in rr:
                for o in r["rna_occurrences"]:
                    crossed = crossed_joins(o)
                    blocks = [{k: b[k] for k in ("contig", "strand", "reference_start", "reference_end")}
                              for b in o["genomic_blocks"]]
                    key = identity([blocks, [junction_text(j) for j in crossed], o["frame_status"], start_text(o["start_evidence"])])
                    placements[key] = (o, crossed)
            out += ["Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). "
                    "Blocks describe this ORF interval, not a proven full mature transcript:", "",
                    table(["ORF genomic blocks in RNA order", "Crossed RNA join(s)", "Frame / initiation"],
                          [(" → ".join(map(locus, o["genomic_blocks"])), "<br>".join(map(junction_text, js)),
                            "%s / %s" % (o["frame_status"], start_text(o["start_evidence"]))) for o, js in placements.values()]),
                    "Flags: `" + str(sorted({f for r in rr for f in r["uncertainty_flags"]})) + "`.", ""]
        else:
            out += ["Compatible reference transcripts: " + ", ".join(sorted({t for r in rr for t in r["transcript_ids"]})) + ".", ""]
    return "\n".join(out)


def publish(ledger, output):
    output.mkdir(parents=True, exist_ok=True)
    (output / "events").mkdir(exist_ok=True)
    write_json(output / "events.json.gz", ledger)
    rows = ledger["candidates"]
    index = []
    with (output / "proteins.fasta").open("w") as fa, (output / "sv-orfs.fna").open("w") as nt:
        for seq in sorted(ledger["sequences"].values(), key=lambda s: s["id"]):
            fa.write(">%s event=%s kind=%s terminal_stop=%s\n%s\n" % (
                seq["id"], seq["event_id"], seq["kind"], seq["ends_with_stop_codon"],
                "\n".join(textwrap.wrap(seq["amino_acids"] + ("*" if seq["ends_with_stop_codon"] else ""), 80))))
            if seq["nucleotide_sequence"]:
                nt.write(">%s event=%s\n%s\n" % (seq["id"], seq["event_id"],
                         "\n".join(textwrap.wrap(seq["nucleotide_sequence"], 80))))
    with (output / "support.tsv").open("w") as handle:
        fields = ["event_id", "report_sequence_id", "kind", "source_id", "product", "orientation", "amino_acids",
                  "ends_with_stop_codon", "mutation_interval", "fragments", "q20_full_window_fragments",
                  "contributing_fragments", "cells", "missing_quality_reads", "transcript_ids", "result_path"]
        fields.insert(fields.index("transcript_ids"), "fragments_with_unavailable_qualities")
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: PRODUCTS[row["source_id"]] if k == "product" else row.get(k, "") for k in fields})
    for eid, event in sorted(ledger["events"].items()):
        (output / "events" / (eid + ".md")).write_text(event_markdown(event, ledger))
        event_rows = [r for r in rows if r["event_id"] == eid]
        index.append(("[%s](events/%s.md)" % (" / ".join(event["names"]), eid), event["kind"], event["dna_origin"],
                      len(event_rows), len({r["report_sequence_id"] for r in event_rows})))
    (output / "README.md").write_text(
        "# Sid DNA → RNA → ORF event sheets\n\n"
        "[Interpretation and counting definitions](../../sid-neoorf-event-report.md). "
        "All %d nominations/geometry entries, including unresolved and zero-candidate outcomes, are below. " % len(ledger["events"]) +
        "Sequence groups distinguish event, frame/window, stop and SV nucleotide sequence; they are not independent neoORFs.\n\n"
        "[Full proteins/windows FASTA](proteins.fasta), [SV ORF nucleotides](sv-orfs.fna), "
        "[all source-specific support rows](support.tsv), [complete event/placement ledger](events.json.gz).\n\n" +
        table(["Event sheet", "Type", "DNA provenance", "Source hypotheses", "Sequence groups"], index))
    write_json(output / "files.json", {str(p.relative_to(output)): digest(p)
                                      for p in sorted(output.rglob("*")) if p.is_file() and p.name != "files.json"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--audit", type=Path, action="append", required=True)
    parser.add_argument("--checks", type=Path, required=True)
    args = parser.parse_args()
    ledger = build(args.screen, args.audit, args.checks)
    publish(ledger, args.output)
    print(ledger["counts"])


if __name__ == "__main__":
    main()
