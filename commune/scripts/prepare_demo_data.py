#!/usr/bin/env python3
"""Create a reproducible, small COMMUNE input dataset from a backed h5ad."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


DEFAULT_IMMUNE_TYPES = {
    "T cells",
    "Myeloid cells",
    "B cells",
    "NK cells",
    "NKT cells",
    "Plasma cells",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Source h5ad file")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--groups", nargs="+", default=["Normal", "Tumor"])
    parser.add_argument("--samples-per-group", type=int, default=20)
    parser.add_argument("--cells-per-sample", type=int, default=100)
    parser.add_argument("--min-cells-per-type", type=int, default=3)
    parser.add_argument("--max-celltypes-per-sample", type=int, default=12)
    parser.add_argument("--exclude-tissue", nargs="*", default=["Blood"])
    parser.add_argument("--seed", type=int, default=20250713)
    parser.add_argument(
        "--no-h5ad",
        action="store_true",
        help="Do not write the convenience h5ad copy",
    )
    return parser.parse_args()


def choose_samples(
    obs: pd.DataFrame,
    groups: list[str],
    samples_per_group: int,
    cells_per_sample: int,
    exclude_tissue: list[str],
    rng: np.random.Generator,
) -> pd.DataFrame:
    required = {"sample", "group", "celltype_1", "celltype_2", "celltype_3"}
    missing = required.difference(obs.columns)
    if missing:
        raise ValueError(f"Missing required obs columns: {sorted(missing)}")

    work = obs.loc[obs["group"].astype(str).isin(groups)].copy()
    if exclude_tissue and "tissue_organ" in work:
        work = work.loc[~work["tissue_organ"].astype(str).isin(exclude_tissue)]
    work = work.loc[
        work["sample"].notna()
        & work["celltype_3"].notna()
        & work["celltype_3"].astype(str).ne("")
        & work["celltype_3"].astype(str).ne("nan")
    ].copy()
    work["is_immune"] = work["celltype_1"].astype(str).isin(DEFAULT_IMMUNE_TYPES)

    stats = (
        work.groupby(["group", "sample"], observed=True)
        .agg(
            total_cells=("celltype_3", "size"),
            n_celltypes=("celltype_3", "nunique"),
            immune_cells=("is_immune", "sum"),
        )
        .reset_index()
    )
    stats["nonimmune_cells"] = stats["total_cells"] - stats["immune_cells"]
    eligible = stats.loc[
        (stats["total_cells"] >= cells_per_sample)
        & (stats["immune_cells"] > 0)
        & (stats["nonimmune_cells"] > 0)
        & (stats["n_celltypes"] >= 3)
    ].copy()

    selected = []
    for group in groups:
        candidates = eligible.loc[eligible["group"].astype(str).eq(group)].copy()
        if len(candidates) < samples_per_group:
            raise ValueError(
                f"Group {group!r} has only {len(candidates)} eligible samples; "
                f"requested {samples_per_group}."
            )
        # Prefer diverse samples, but randomize ties reproducibly.
        candidates["tie_break"] = rng.random(len(candidates))
        candidates = candidates.sort_values(
            ["n_celltypes", "total_cells", "tie_break"],
            ascending=[False, False, True],
        )
        selected.append(candidates.head(samples_per_group))
    return pd.concat(selected, ignore_index=True)


def sample_cells(
    obs: pd.DataFrame,
    selected_samples: pd.DataFrame,
    cells_per_sample: int,
    min_cells_per_type: int,
    max_celltypes_per_sample: int,
    rng: np.random.Generator,
) -> np.ndarray:
    chosen: list[int] = []
    for row in selected_samples.itertuples(index=False):
        mask = obs["group"].astype(str).eq(str(row.group)) & obs["sample"].astype(
            str
        ).eq(str(row.sample))
        sample_obs = obs.loc[mask]
        celltype_counts = sample_obs["celltype_3"].astype(str).value_counts()
        ordered_celltypes = celltype_counts.index.tolist()
        max_types = min(
            max_celltypes_per_sample,
            max(1, cells_per_sample // max(1, min_cells_per_type)),
        )
        immune_mask = sample_obs["celltype_1"].astype(str).isin(DEFAULT_IMMUNE_TYPES)
        immune_types = set(sample_obs.loc[immune_mask, "celltype_3"].astype(str))
        nonimmune_types = set(sample_obs.loc[~immune_mask, "celltype_3"].astype(str))
        required_types = []
        for candidates in (immune_types, nonimmune_types):
            required = next((ct for ct in ordered_celltypes if ct in candidates), None)
            if required is not None:
                required_types.append(required)
        celltypes = list(dict.fromkeys(required_types + ordered_celltypes))[:max_types]

        sample_chosen: list[int] = []
        for celltype in celltypes:
            indices = sample_obs.index[
                sample_obs["celltype_3"].astype(str).eq(celltype)
            ].to_numpy(dtype=np.int64)
            take = min(min_cells_per_type, len(indices))
            sample_chosen.extend(rng.choice(indices, size=take, replace=False).tolist())

        sample_chosen = list(dict.fromkeys(sample_chosen))
        remaining_n = cells_per_sample - len(sample_chosen)
        if remaining_n > 0:
            eligible_remaining = sample_obs.index[
                sample_obs["celltype_3"].astype(str).isin(celltypes)
            ].to_numpy(dtype=np.int64)
            remaining = np.setdiff1d(
                eligible_remaining,
                np.asarray(sample_chosen, dtype=np.int64),
                assume_unique=False,
            )
            take = min(remaining_n, len(remaining))
            sample_chosen.extend(rng.choice(remaining, size=take, replace=False).tolist())
        chosen.extend(sample_chosen[:cells_per_sample])
    return np.asarray(sorted(chosen), dtype=np.int64)


def read_csr_rows(h5ad_path: Path, rows: np.ndarray) -> sparse.csr_matrix:
    with h5py.File(h5ad_path, "r") as handle:
        x = handle["X"]
        if not isinstance(x, h5py.Group) or not {"data", "indices", "indptr"}.issubset(x):
            raise ValueError("This utility currently requires CSR-encoded h5ad X data.")
        shape = tuple(int(v) for v in x.attrs.get("shape", ()))
        n_vars = shape[1] if len(shape) == 2 else int(handle["var"]["_index"].shape[0])
        source_indptr = x["indptr"]
        data_parts = []
        index_parts = []
        output_indptr = np.zeros(len(rows) + 1, dtype=np.int64)
        for output_row, source_row in enumerate(rows, start=1):
            start = int(source_indptr[source_row])
            end = int(source_indptr[source_row + 1])
            data_parts.append(x["data"][start:end])
            index_parts.append(x["indices"][start:end])
            output_indptr[output_row] = output_indptr[output_row - 1] + (end - start)
        data = np.concatenate(data_parts) if data_parts else np.array([], dtype=np.int64)
        indices = (
            np.concatenate(index_parts) if index_parts else np.array([], dtype=np.int64)
        )
    return sparse.csr_matrix((data, indices, output_indptr), shape=(len(rows), n_vars))


def write_gzipped_mtx(matrix: sparse.spmatrix, output: Path) -> None:
    plain = output.with_suffix("") if output.suffix == ".gz" else output
    mmwrite(plain, matrix)
    if output.suffix == ".gz":
        with plain.open("rb") as source, gzip.open(output, "wb", compresslevel=6) as target:
            shutil.copyfileobj(source, target)
        plain.unlink()


def main() -> None:
    args = parse_args()
    # anndata 0.10 writes NumPy-backed strings, whereas pandas 3 defaults to
    # Arrow-backed strings that this anndata release cannot serialize.
    pd.options.future.infer_string = False
    input_path = Path(args.input).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    backed = ad.read_h5ad(input_path, backed="r")
    obs = backed.obs.copy()
    obs.index = np.arange(len(obs), dtype=np.int64)
    var_names = backed.var_names.astype(str).to_numpy()
    original_cell_ids = backed.obs_names.astype(str).to_numpy()
    backed.file.close()

    selected_samples = choose_samples(
        obs,
        args.groups,
        args.samples_per_group,
        args.cells_per_sample,
        args.exclude_tissue,
        rng,
    )
    selected_rows = sample_cells(
        obs,
        selected_samples,
        args.cells_per_sample,
        args.min_cells_per_type,
        args.max_celltypes_per_sample,
        rng,
    )
    matrix = read_csr_rows(input_path, selected_rows)

    selected_obs = obs.loc[selected_rows].copy()
    selected_obs.index = pd.Index(
        np.asarray(original_cell_ids[selected_rows], dtype=object), dtype=object
    )
    selected_obs.index.name = "cell_id"
    selected_obs = selected_obs.drop(columns=["is_immune"], errors="ignore")
    for column in selected_obs.columns:
        if isinstance(selected_obs[column].dtype, pd.CategoricalDtype):
            selected_obs[column] = selected_obs[column].map(str).astype(object)
        elif pd.api.types.is_string_dtype(selected_obs[column].dtype):
            selected_obs[column] = selected_obs[column].map(str).astype(object)

    selected_samples.to_csv(output_dir / "selected_samples.csv", index=False)
    selected_obs.to_csv(output_dir / "metadata.csv")
    pd.Series(var_names, name="gene").to_csv(
        output_dir / "genes.tsv", sep="\t", index=False, header=False
    )
    pd.Series(selected_obs.index, name="cell_id").to_csv(
        output_dir / "cells.tsv", sep="\t", index=False, header=False
    )
    write_gzipped_mtx(matrix.transpose().tocsc(), output_dir / "counts.mtx.gz")

    if not args.no_h5ad:
        demo = ad.AnnData(
            X=matrix,
            obs=selected_obs,
            var=pd.DataFrame(index=pd.Index(var_names, name="gene")),
        )
        demo.write_h5ad(output_dir / "commune_demo.h5ad", compression="gzip")

    manifest = {
        "source": str(input_path),
        "groups": args.groups,
        "samples_per_group": args.samples_per_group,
        "cells_per_sample": args.cells_per_sample,
        "max_celltypes_per_sample": args.max_celltypes_per_sample,
        "seed": args.seed,
        "shape_cells_by_genes": list(matrix.shape),
        "nonzero_entries": int(matrix.nnz),
        "selected_samples": int(selected_obs["sample"].nunique()),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
