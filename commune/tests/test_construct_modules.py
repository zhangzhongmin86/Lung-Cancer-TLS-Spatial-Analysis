from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


def test_construct_modules_uses_adjusted_cooccurrence_table(tmp_path: Path) -> None:
    modules = pd.DataFrame(
        {
            "interaction_key": ["A|||B|||L1|||R1", "A|||C|||L2|||R2"],
            "module": [1, 1],
            "weight": [0.8, 0.7],
            "source": ["A", "A"],
            "target": ["B", "C"],
            "ligand": ["L1", "L2"],
            "receptor": ["R1", "R2"],
        }
    )
    metadata = pd.DataFrame(
        {
            "celltype_2": ["Parent1", "Parent2", "Parent3"],
            "celltype_3": ["A", "B", "C"],
        },
        index=["cell1", "cell2", "cell3"],
    )
    cooccurrence = pd.DataFrame(
        {
            "Cell_Type_1": ["B"],
            "Cell_Type_2": ["A"],
            "Correlation": [0.6],
            "P_value": [0.001],
            "N": [20],
            "Adjusted_p_value": [0.01],
            "Correlation_type": ["Moderate positive"],
            "group": ["Tumor_metastasis"],
        }
    )

    modules_file = tmp_path / "modules.csv"
    metadata_file = tmp_path / "metadata.csv"
    cooccurrence_dir = tmp_path / "cooccurrence"
    output_dir = tmp_path / "output"
    cooccurrence_dir.mkdir()
    modules.to_csv(modules_file, index=False)
    metadata.to_csv(metadata_file)
    cooccurrence.to_csv(cooccurrence_dir / "cooccurrence_all_groups.csv", index=False)

    script = Path(__file__).parents[1] / "scripts" / "construct_modules.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--modules",
            str(modules_file),
            "--cooccurrence-dir",
            str(cooccurrence_dir),
            "--metadata",
            str(metadata_file),
            "--output-dir",
            str(output_dir),
        ],
        check=True,
    )

    constrained = pd.read_csv(output_dir / "constrained_module_interactions.csv")
    assert constrained["interaction_key"].tolist() == ["A|||B|||L1|||R1"]
    assert constrained["supported_groups"].tolist() == ["Tumor_metastasis"]
    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["retained_modules"] == 1
    assert manifest["retained_interactions"] == 1
