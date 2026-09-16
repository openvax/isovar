"""Rendering must preserve unknown coding status and distinct coordinates."""
from dataclasses import asdict
import json

import pytest

from isovar.cli import commands
from isovar.fusion import reconstruct_fusion
from isovar.fusion_visualization import fusion_figures, save_fusion_figures
from examples.osteosarc_context_figures import load_context, verify_context, context_panels
from .test_fusion import example


def test_fusion_panels_use_actual_boundaries_and_do_not_invent_proteins(tmp_path):
    pytest.importorskip("matplotlib")
    fusion,refs,reads=example(insert="GCT")
    result=reconstruct_fusion(fusion,refs,reads,peptide_lengths=(2,3))
    panels=dict(fusion_figures(result,refs,{"donor.1":"Donor-name"}))
    assert set(panels)=={"junction","donor-context","acceptor-context","protein-1","junction-peptides-1"}
    assert panels["junction"].get_facecolor()==(1,1,1,1)
    assert any("Donor-name" in t.get_text() for t in panels["donor-context"].axes[0].texts)
    assert len(panels["donor-context"].axes[0].patches)==1
    result=reconstruct_fusion(fusion,(),reads)
    assert result["status"]=="unresolved_frame"
    assert set(dict(fusion_figures(result,refs)))=={"junction","donor-context","acceptor-context"}
    directory=save_fusion_figures(result,tmp_path,iter(refs),dpi=80)
    assert {p.stem for p in directory.glob("*.png")}=={"junction","donor-context","acceptor-context"}
    assert len(list(directory.glob("*.svg")))==3
    assert (directory/"all-figures.pdf").read_bytes().startswith(b"%PDF")
    assert json.loads((directory/"evidence.json").read_text())["translations"]==[]
    assert len(json.loads((directory/"reference-models.json").read_text())["references"])==2
    with pytest.raises(FileExistsError):
        save_fusion_figures(result,tmp_path,refs,dpi=80)
    with pytest.raises(ValueError,match="dpi"):
        save_fusion_figures(result,tmp_path/"invalid",refs,dpi=0)


def test_fusion_cli_figures_are_optional_and_timestamped(tmp_path):
    pytest.importorskip("matplotlib")
    fusion,refs,reads=example()
    source=tmp_path/"input.json"
    source.write_text(json.dumps(dict(fusion=asdict(fusion),references=[asdict(r) for r in refs],reads=[asdict(r) for r in reads])))
    commands.run(["fusion","--input",str(source),"--output",str(tmp_path/"result.json"),
                  "--plot-dir",str(tmp_path/"figures"),"--dpi","80"])
    assert len(list((tmp_path/"figures").glob("20*Z/*/all-figures.pdf")))==1


def test_context_counts_and_protein_are_reproduced_from_original_records():
    data=load_context()
    protein=verify_context(data)
    assert len(protein["compound"])==74
    assert protein["mutation_offsets"]==[2,71]
    assert protein["single"][2]=="T" and protein["compound"][2]=="M"
    assert protein["reference"][71]=="R" and protein["compound"][71]=="M"
    assert [r["minimum_quality"] for r in protein["witnesses"]]==[32,23,23]
    pytest.importorskip("matplotlib")
    panels=list(context_panels(data,protein))
    assert len(panels)==6
    assert all(fig.get_facecolor()==(1,1,1,1) for _,_,fig in panels)
