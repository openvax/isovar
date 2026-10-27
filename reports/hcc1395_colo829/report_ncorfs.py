"""Reference ncORF frames reconstructed from the original mutation-support RNA."""

import gzip
import json

from reportlab.platypus import Flowable, PageBreak, Spacer

from .build_evidence_pdfs import REPORT, TEAL, GRAY, PALE, ORANGE, para, section, data_table
from .report_landscape import NeutralReferenceComparison


def results_for(sample):
    with gzip.open(REPORT / 'ncorf_rna_audit/results.json.gz', 'rt') as stream:
        return [r for r in json.load(stream) if r['sample_id'] == sample]


def choose(candidates, prefer_twenty_spans=False):
    candidates = [p for p in candidates if p['focal_amino_acid_change']]
    if prefer_twenty_spans:
        strong = [p for p in candidates if p['coding_support']['single_read_span_support']['fragments'] >= 20]
        if strong:
            return max(strong, key=lambda p: (len(p['protein_sequence']), p['rna_support']['fragments']))
    return max(candidates, key=lambda p: (p['rna_support']['fragments'], len(p['protein_sequence']))) if candidates else None


def changed_default_examples(sample):
    examples = []
    for result in results_for(sample):
        if result['pass_name'] != 'default_coverage':
            continue
        for model in result['models']:
            if model['reference']['kind'] != 'extra_orf':
                continue
            candidate = choose([p for p in result['protein_candidates'] if p['model_id'] == model['model_id']], True)
            if candidate:
                examples.append((result, model, candidate))
    return examples


class CodingCoverage(Flowable):
    """Coverage only where the displayed coding sequence was observed."""
    def __init__(self, candidate, reference_length):
        super().__init__()
        self.candidate, self.reference_length = candidate, reference_length
        self.width, self.height = 512, 92

    def draw(self):
        c, p = self.canv, self.candidate
        c.translate(0, 15)
        n = len(p['protein_sequence'])
        depths = p['coding_support']['fragment_depth'][:3 * (n + int(p['ends_with_stop_codon']))]
        per_codon = [min(depths[i:i + 3]) for i in range(0, len(depths), 3)]
        maximum = max(per_codon, default=1)
        c.setFillColor(GRAY)
        c.setFont('Helvetica', 7.5)
        c.drawString(0, 65, 'RNA fragment support per codon (minimum base depth; stop included when observed)')
        for i, depth in enumerate(per_codon):
            focal = p['mutation_interval'][0] <= i < p['mutation_interval'][1]
            focal_stop = p['ends_with_stop_codon'] and i == n and p['mutation_interval'] == [n, n]
            c.setFillColor(ORANGE if focal or focal_stop else TEAL)
            c.rect(53 + i * 5.4, 25, 4.7, 29 * depth / max(1, maximum), fill=1, stroke=0)
        c.setFillColor(GRAY)
        c.drawRightString(46, 25, '0')
        c.drawRightString(46, 49, str(maximum))
        c.setFillColor(PALE)
        c.rect(53, 8, 450, 7, fill=1, stroke=0)
        start = p['reference_protein_start_0based']
        c.setFillColor(TEAL)
        c.rect(53 + 450 * start / self.reference_length, 8,
               450 * min(n, self.reference_length - start) / self.reference_length, 7, fill=1, stroke=0)
        c.setFillColor(GRAY)
        c.drawString(0, 7, 'Extent')
        c.drawString(53, -2, f'RNA {n} aa from reference position {start + 1}; reference ORF {self.reference_length} aa')
        c.drawRightString(503, -2, 'No reference filling')


def ncorf_story(sample, data):
    results = results_for(sample)
    default = {r['call_id']: r for r in results if r['pass_name'] == 'default_coverage'}
    sensitivity = {r['call_id']: r for r in results if r['pass_name'] == 'single_read_sensitivity'}
    rows, strong, weak = [], [], []
    model_count = 0
    for ident, result in sorted(default.items(), key=lambda pair: -pair[1]['atlas_alt_fragments']):
        exploratory = sensitivity[ident]
        for model in result['models']:
            reference = model['reference']
            if reference['kind'] != 'extra_orf':
                continue
            model_count += 1
            candidates = [p for p in result['protein_candidates'] if p['model_id'] == model['model_id']]
            low = [p for p in exploratory['protein_candidates'] if p['model_id'] == model['model_id']]
            best = choose(candidates, True)
            weak_best = choose(low)
            if best:
                strong.append((result, model, best))
                state = f"Changed local protein; {best['rna_support']['fragments']} assembled fragments"
            elif candidates:
                state = 'Focal allele is synonymous in this ORF'
            elif model['reason'] == 'reference_start_affected':
                state = 'Reference start codon affected; no new initiation inferred'
            elif weak_best:
                weak.append((exploratory, model, weak_best))
                state = f"Sensitivity only: {weak_best['rna_support']['fragments']}-fragment candidate"
            else:
                state = ('RNA path incompatible with this isoform' if model['reason'] == 'no_compatible_rna_reads'
                         else 'No retained RNA/frame match under either budget')
            rows.append((f"<b>{reference['source_gene_name']}</b><br/>{result['igv_locus']}",
                         reference['orf_id'], str(result['atlas_alt_fragments']), state))
    story = section('ncORFs and other antigen-source categories')
    story.append(para(f'<b>{len(default)} ALT-supported calls → {model_count} validated ncORF models → '
                      f'{len(strong)} amino-acid-changing models at default coverage.</b> '
                      f'The separate single-read sensitivity pass adds {len(weak)} low-support model(s). '
                      'All input models are assessed; ordinary CDS competitors remain in the exported results. '
                      'Alternative ORFs on one call reuse its molecules and are not separate mutations.'))
    story.append(data_table(['Gene / IGV (GRCh38)', 'Reference ORF', 'ALT', 'RNA result'], rows,
                            [186, 100, 34, 192]))
    sources = para('Sources: ORFdata 7,264-model catalogue 2026-10-11.1; complete ORF paths/proteins validated '
                   'against Ensembl 87. COLO829 uses the whole-run selected-site BAM. HCC1395 uses all 36 '
                   'original ALT witnesses for these loci, recovered from SRA and realigned genome-wide with '
                   'minimap2 2.26; all five ALT counts reproduce. Default minimum read coverage is 2; the '
                   'sensitivity pass uses 1 and reuses the same molecules. Full sequences, model alternatives, '
                   'original witnesses, source hashes and budgets: '
                   '<link href="https://github.com/openvax/isovar/blob/master/reports/hcc1395_colo829/'
                   'ncorf_rna_audit/manifest.json">ncORF audit manifest</link>.', 'FootSmall')
    limits = para('Assessment limits: RNA and reference Ribo-seq/MS annotations do not measure mutant translation '
                  'or HLA presentation. General cryptic initiation and retained-intron ORFs remain unresolved '
                  'beyond the focused DDX17 audit. Viral RNA/ORFs and cancer/testis isoform surveys are not assessed. '
                  'Germline/co-somatic colors require source-qualified genotype and explicit phasing.', 'FootSmall')
    story.append(Spacer(1, 8))
    if sample == 'HCC1395':
        story += section('What prevents stronger protein calls')
        story.append(para('The 22 NELFB ALT fragments change a reference start codon. They support the RNA allele, '
                          'but cannot establish an alternative initiation site. The ordinary coding-frame '
                          'competitor remains in the export. For UQCC1 and ST3GAL1, relaxing coverage alone '
                          'does not produce a retained RNA/frame match. MIR762HG supports one catalogue '
                          'isoform path; the other path is incompatible with these reads.'))
        story.append(para('The two sensitivity examples on the next page each depend on one original ONT '
                          'molecule. Their nonfocal differences can include sequencing errors. They identify '
                          'follow-up targets, not high-confidence new proteins. The recovered input is '
                          'ALT-witness selected, so it supplies no reference-depth or allele-frequency denominator.'))
        story.extend([sources, limits, PageBreak()])
        story += section('Low-support ncORF sensitivity examples')
        story.append(para('Independent confirmation is required for both candidates. Reference residues remain '
                          'neutral; only the focal RNA change is orange.'))
    examples = strong if strong else sorted(weak, key=lambda x: -len(x[2]['protein_sequence']))
    for result, model, candidate in examples[:2]:
        reference = model['reference']
        extent = 'complete RNA coding candidate' if candidate['starts_at_reference_start'] and candidate['ends_with_stop_codon'] else 'partial RNA coding candidate'
        support = candidate['coding_support']
        story.append(para(f"<b>{reference['source_gene_name']} | {reference['orf_id']} | {extent}</b>"))
        story.append(para(f"{candidate['rna_support']['fragments']} assembled fragments; "
                          f"{support['single_read_span_support']['fragments']} fragments each span the displayed coding interval; "
                          f"minimum base depth {support['minimum_fragment_depth']}. "
                          f"{result['source']}; {result['pass_name'].replace('_', ' ')}.", 'FootSmall'))
        sequence = candidate['protein_sequence'] + ('*' if candidate['ends_with_stop_codon'] else '')
        first, last = candidate['mutation_interval']
        if candidate['ends_with_stop_codon'] and first == last == len(candidate['protein_sequence']):
            last += 1  # An observed focal stop has a symbol, but no amino-acid residue.
        story.append(NeutralReferenceComparison(candidate['reference_protein_window'], sequence, first, last))
        story.append(CodingCoverage(candidate, len(reference['protein_sequence'])))
        if sample == 'COLO829':
            story.append(para('Complete reference ORF (RNA above is not extended from this sequence):', 'FootSmall'))
            sequence = reference['protein_sequence']
            story.append(para('<br/>'.join(sequence[i:i + 70] for i in range(0, len(sequence), 70)), 'TinyMono'))
            story.append(para('The observed start-to-residue-15 window has G5E relative to this 85-aa downstream ORF. '
                              'A longer 23-aa candidate retains 28 assembled fragments but only 11 single-read spans '
                              '(minimum base depth 15). The overlapping c17riboseqorf123 frame and PITPNM1 uORF '
                              'are synonymous for their focal alleles; they are not counted as mutant proteins.', 'FootSmall'))
        else:
            story.append(para('Complete reference ORF: ' + reference['protein_sequence'], 'TinyMono'))
            story.append(para('Single-read hypothesis. Nonfocal differences remain unclassified; a premature RNA stop '
                              'or an apparent complete coding interval does not establish a translated product. '
                              'Independent original-RNA confirmation is still needed.', 'FootSmall'))
        story.append(Spacer(1, 4))
    story.append(para('<font color="#ba5a29">&#9632;</font> Focal change &nbsp; '
                      '<font color="#526373">&#9632;</font> Unclassified nonfocal difference &nbsp; '
                      '<font color="#0f706d">&#9632;</font> RNA-covered interval. '
                      '* Observed stop. Reference residues stay neutral. Coverage outside the displayed RNA is unassessed.', 'FootSmall'))
    if sample == 'COLO829':
        story.extend([sources, limits])
    story.append(PageBreak())
    return story
