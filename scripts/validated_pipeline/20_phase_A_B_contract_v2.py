#!/usr/bin/env python
"""Phase A-B v2 audit for the PC/NC snRNA-seq input.

This script is deliberately read-only for raw data and protected submission
assets. It establishes the v2 contract before any v2 pseudobulk analysis.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


_context_root = os.environ.get("IUA_CONTEXT_ROOT")
if not _context_root:
    raise RuntimeError("Set IUA_CONTEXT_ROOT to the audited context-validation workspace before running this provenance script.")
ROOT = Path(_context_root).expanduser().resolve()
PROJECT = ROOT.parents[1]
H5AD = ROOT / "raw" / "PC_NC_JCI2026" / "natproendometrium.h5ad"
VERSION = "IUA_CONTEXT_REPAIR_V2"
OUTCOME_PRIMARY = ["PGR", "frozen_candidate25_score"]
OUTCOME_SECONDARY = ["ESR1", "stable_epithelial_core_score"]
FROZEN_PATH = PROJECT / "results" / "reference_WOI_program_epithelial.tsv"
CORE_PATH = PROJECT / "results" / "reference_epithelial_WOI_program_stability_core_cell_specificity.tsv"
MARKERS = {
    "epithelial_identity": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"],
    "glandular_identity": ["PAX8", "MUC1", "WFDC2", "SCGB2A1"],
    "stromal_exclusion": ["DCN", "LUM", "COL1A1", "COL3A1"],
    "immune_exclusion": ["PTPRC", "LYZ", "CD3D", "MS4A1"],
    "endothelial_exclusion": ["PECAM1", "VWF", "KDR", "EMCN"],
}


def stamp() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path, block: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_obs(h5: h5py.File, key: str) -> np.ndarray:
    obj = h5[f"obs/{key}"]
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        cats = [decode(x) for x in obj["categories"][...]]
        codes = obj["codes"][...]
        return np.asarray([cats[int(i)] if int(i) >= 0 else "NA" for i in codes], dtype=object)
    return np.asarray([decode(x) for x in obj[...]], dtype=object)


def read_var(h5: h5py.File, key: str) -> np.ndarray:
    obj = h5[f"raw/var/{key}"]
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        cats = [decode(x) for x in obj["categories"][...]]
        return np.asarray([cats[int(i)] if int(i) >= 0 else "NA" for i in obj["codes"][...]], dtype=object)
    return np.asarray([decode(x) for x in obj[...]], dtype=object)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> None:
    for name in ("metadata", "interim", "results", "reports", "logs", "checkpoints", "figures", "scripts"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    if not H5AD.exists():
        raise FileNotFoundError(H5AD)
    if not FROZEN_PATH.exists() or not CORE_PATH.exists():
        raise FileNotFoundError("Frozen outcome definitions are missing from the preserved project.")

    # Exact protected-asset inventory. These are verification records, not copies.
    protected_roots = [PROJECT / "投稿", PROJECT / "湖南省科技创新计划项目申报书（宫腔粘连）.pdf"]
    protected = []
    for root in protected_roots:
        files = list(root.rglob("*")) if root.is_dir() else [root]
        for item in files:
            if item.is_file():
                protected.append({
                    "asset_role": "submission_asset" if "投稿" in str(item) else "administrative_source",
                    "path": str(item), "size_bytes": item.stat().st_size,
                    "sha256": sha256(item), "verified_at": stamp(), "mutation_policy": "read_only",
                })
    write_tsv(pd.DataFrame(protected), ROOT / "metadata" / "protected_assets_v2.tsv")

    # v1 outputs are retained, but not blindly elevated to v2 conclusions.
    v1_files = [
        ROOT / "results" / "PC_NC_broad_epithelial_receptor_effects.tsv",
        ROOT / "results" / "PC_NC_broad_epithelial_programme_tests.tsv",
        ROOT / "results" / "PC_NC_broad_epithelial_edgeR_QL_all_genes.tsv",
        ROOT / "interim" / "PC_NC_broad_epithelial_pseudobulk_counts.tsv.gz",
        PROJECT / "results" / "GSE215968_disjoint_gate_receptor_effects.tsv",
        PROJECT / "results" / "GSE215968_disjoint_gate_programme_tests.tsv",
    ]
    reuse = []
    for item in v1_files:
        if item.exists():
            reuse.append({
                "path": str(item), "sha256": sha256(item), "size_bytes": item.stat().st_size,
                "v2_role": "historical_audit_only" if "PC_NC_broad" in item.name else "preserved_cross_context_source",
                "reason": "Broad PC/NC epithelial compartment and original BH families require v2 correction." if "PC_NC_broad" in item.name else "No new donor independence is implied by reuse.",
            })
    write_tsv(pd.DataFrame(reuse), ROOT / "metadata" / "v1_reusable_artifacts_v2.tsv")

    frozen = pd.read_csv(FROZEN_PATH, sep="\t")
    core = pd.read_csv(CORE_PATH, sep="\t")
    frozen_genes = frozen["gene"].astype(str).tolist()
    core_flag = "epithelial_enriched_2fold_reference_rule"
    if core_flag not in core.columns:
        raise ValueError(f"Missing frozen core flag: {core_flag}")
    core_genes = core.loc[core[core_flag].astype(bool), "gene"].astype(str).tolist()
    tested_set = set(["ESR1", "PGR", *frozen_genes, *core_genes])
    marker_rows = []
    for panel, genes in MARKERS.items():
        for gene in genes:
            marker_rows.append({"panel": panel, "gene": gene, "overlaps_frozen_outcome_set": gene in tested_set})
    marker_frame = pd.DataFrame(marker_rows)
    write_tsv(marker_frame, ROOT / "metadata" / "PC_NC_outcome_marker_set_overlap_v2.tsv")
    if marker_frame["overlaps_frozen_outcome_set"].any():
        raise ValueError("Annotation marker panel overlaps a tested endpoint gene; revise the fixed panel.")

    with h5py.File(H5AD, "r") as h5:
        donor = read_obs(h5, "donor_id")
        cycle = read_obs(h5, "CycleType")
        cell_type = read_obs(h5, "cell_type")
        genes = read_var(h5, "gene_symbols")
        raw = h5["raw/X"]
        indptr = raw["indptr"]
        indices = raw["indices"]
        values = raw["data"]
        n_cells, n_genes = h5["raw/X"].attrs["shape"] if "shape" in h5["raw/X"].attrs else (len(indptr) - 1, len(genes))
        if int(n_cells) != len(donor) or int(n_genes) != len(genes):
            raise ValueError("AnnData axes, donor metadata, and raw feature axis are inconsistent.")
        if len(set(genes)) != len(genes):
            raise ValueError("Gene symbols are not unique; no symbol-level pseudobulk should proceed.")
        if set(np.unique(cycle)) != {"Controlled", "Natural"}:
            raise ValueError(f"Unexpected CycleType values: {sorted(set(cycle))}")
        donor_cycle = pd.DataFrame({"donor_id": donor, "CycleType": cycle}).drop_duplicates()
        if donor_cycle["donor_id"].duplicated().any():
            raise ValueError("A donor maps to more than one CycleType.")
        donor_cycle = donor_cycle.sort_values("donor_id").reset_index(drop=True)
        if donor_cycle["CycleType"].value_counts().to_dict() != {"Natural": 9, "Controlled": 7}:
            raise ValueError("Expected 7 Controlled and 9 Natural donors were not observed.")

        # Full stored-value integer semantic check without dense materialization.
        non_integer, negative, nnz = 0, 0, 0
        for start in range(0, len(values), 5_000_000):
            block = values[start:min(start + 5_000_000, len(values))]
            non_integer += int(np.count_nonzero(np.abs(block - np.rint(block)) > 1e-8))
            negative += int(np.count_nonzero(block < 0))
            nnz += len(block)
        if non_integer or negative:
            raise ValueError(f"raw/X is not a non-negative integer count matrix: noninteger={non_integer}, negative={negative}")

        # Cell composition and protocol ledger reflect available H5AD metadata only.
        comp = (pd.DataFrame({"donor_id": donor, "CycleType": cycle, "cell_type": cell_type})
                  .groupby(["donor_id", "CycleType", "cell_type"], as_index=False).size()
                  .rename(columns={"size": "nuclei"}))
        donor_totals = comp.groupby(["donor_id", "CycleType"], as_index=False)["nuclei"].sum().rename(columns={"nuclei": "all_nuclei"})
        comp = comp.merge(donor_totals, on=["donor_id", "CycleType"], validate="many_to_one")
        comp["fraction_of_donor_nuclei"] = comp["nuclei"] / comp["all_nuclei"]
        write_tsv(comp.sort_values(["donor_id", "cell_type"]), ROOT / "metadata" / "PC_NC_cell_composition_by_donor_v2.tsv")
        protocol = donor_cycle.copy()
        protocol["nuclei_total"] = protocol["donor_id"].map(donor_totals.set_index("donor_id")["all_nuclei"])
        protocol["sampling_day_or_WOI_anchor"] = "NA_not_encoded_in_H5AD"
        protocol["hormone_dose_or_regimen"] = "NA_not_encoded_in_H5AD"
        protocol["age"] = "NA_not_encoded_in_H5AD"
        protocol["batch_or_library"] = "NA_not_encoded_in_H5AD"
        protocol["source_metadata_scope"] = "H5AD obs only; protocol requires source-paper verification"
        write_tsv(protocol, ROOT / "metadata" / "PC_NC_donor_protocol_ledger_v2.tsv")

        # Audit the publisher's fine labels with a fixed, outcome-disjoint marker panel.
        marker_genes = marker_frame["gene"].tolist()
        marker_idx = {g: int(np.where(genes == g)[0][0]) for g in marker_genes if np.any(genes == g)}
        label_levels = sorted(pd.unique(cell_type).tolist())
        label_n = pd.Series(cell_type).value_counts().to_dict()
        sums = defaultdict(float); detects = defaultdict(int)
        # Scan raw CSR in row chunks; only the fixed marker features are retained.
        marker_lookup = {v: k for k, v in marker_idx.items()}
        for row_start in range(0, len(donor), 5000):
            row_end = min(row_start + 5000, len(donor))
            data_start, data_end = int(indptr[row_start]), int(indptr[row_end])
            chunk_idx = indices[data_start:data_end]
            chunk_dat = values[data_start:data_end]
            if len(chunk_idx) == 0:
                continue
            widths = np.diff(indptr[row_start:row_end + 1]).astype(np.int64)
            chunk_labels = np.repeat(cell_type[row_start:row_end], widths)
            use = np.isin(chunk_idx, list(marker_lookup))
            for idx, value, label in zip(chunk_idx[use], chunk_dat[use], chunk_labels[use]):
                gene = marker_lookup[int(idx)]
                sums[(label, gene)] += float(value)
                detects[(label, gene)] += 1
        annotation = []
        for label in label_levels:
            for _, row in marker_frame.iterrows():
                gene = row["gene"]
                present = gene in marker_idx
                total = sums[(label, gene)] if present else np.nan
                detected = detects[(label, gene)] if present else np.nan
                annotation.append({
                    "cell_type": label, "panel": row["panel"], "gene": gene, "gene_in_h5ad": present,
                    "nuclei_in_label": int(label_n[label]), "mean_raw_count_per_nucleus": total / label_n[label] if present else np.nan,
                    "detection_fraction": detected / label_n[label] if present else np.nan,
                    "annotation_role": "diagnostic_only_no_relabeling",
                })
        write_tsv(pd.DataFrame(annotation), ROOT / "metadata" / "PC_NC_annotation_marker_audit_v2.tsv")

        gene_audit = []
        for gene in ["ESR1", "PGR", *frozen_genes, *core_genes]:
            gene_audit.append({"gene": gene, "in_raw_gene_symbols": gene in set(genes),
                               "endpoint_membership": "primary" if gene in ["PGR", *frozen_genes] else "secondary"})
        gene_audit = pd.DataFrame(gene_audit).drop_duplicates("gene")
        write_tsv(gene_audit, ROOT / "metadata" / "PC_NC_gene_identifier_audit_v2.tsv")
        if not gene_audit.loc[gene_audit["gene"].isin(["ESR1", "PGR"]), "in_raw_gene_symbols"].all():
            raise ValueError("A receptor endpoint is absent from H5AD raw gene symbols.")

        matrices = pd.DataFrame([{
            "layer": "raw/X", "orientation": "cells_by_genes", "n_cells": int(n_cells), "n_genes": int(n_genes),
            "stored_nnz": int(nnz), "value_dtype": str(values.dtype), "all_stored_values_integer": non_integer == 0,
            "nonnegative": negative == 0, "semantic_decision": "count_compatible_for_donor_pseudobulk_edgeR",
        }])
        write_tsv(matrices, ROOT / "metadata" / "PC_NC_matrix_layer_audit_v2.tsv")
        schema = {
            "analysis_id": VERSION, "h5ad": str(H5AD), "h5ad_sha256": sha256(H5AD), "verified_at": stamp(),
            "raw_layer": "raw/X", "n_cells": int(n_cells), "n_genes": int(n_genes), "raw_nnz": int(nnz),
            "observed_cell_type_levels": label_levels, "donors": int(len(donor_cycle)),
            "group_counts": donor_cycle["CycleType"].value_counts().to_dict(),
            "cell_level_qc_reprocessing": "not_run; dataset uses public publisher annotations and v2 performs an audit without re-filtering/relabeling",
            "annotation_marker_overlap_pass": True,
        }
        (ROOT / "metadata" / "PC_NC_h5ad_schema_audit_v2.json").write_text(json.dumps(schema, indent=2), encoding="utf-8")

    manifest = pd.DataFrame([
        {"input_id": "PC_NC_JCI2026_H5AD", "role": "v2_primary_reference", "path": str(H5AD),
         "source": "CELLxGENE collection e80f5450-0a0f-4312-a728-431d50aeb8fe", "accession": "fb90c70d-5917-4d03-920b-c78b230b51e5",
         "release": "local_download_preexisting", "download_date": "2026-09-19", "genome_build": "NA_in_H5AD",
         "organism": "Homo sapiens", "ancestry": "NA", "tissue": "endometrium", "assay": "snRNA-seq", "sample_size": "16 donors",
         "sha256": sha256(H5AD), "notes": "raw/X verified count-compatible; 7 Controlled and 9 Natural donors"},
        {"input_id": "frozen_WOI25", "role": "primary_program_definition", "path": str(FROZEN_PATH),
         "source": "preserved_project_reference", "accession": "NA", "release": "frozen_pre_v2", "download_date": "NA", "genome_build": "NA",
         "organism": "Homo sapiens", "ancestry": "NA", "tissue": "endometrial epithelium", "assay": "gene_set", "sample_size": str(len(frozen_genes)),
         "sha256": sha256(FROZEN_PATH), "notes": "Frozen before v2 analysis"},
        {"input_id": "stable_epithelial_core", "role": "secondary_program_definition", "path": str(CORE_PATH),
         "source": "preserved_project_reference", "accession": "NA", "release": "frozen_pre_v2", "download_date": "NA", "genome_build": "NA",
         "organism": "Homo sapiens", "ancestry": "NA", "tissue": "endometrial epithelium", "assay": "gene_set", "sample_size": str(len(core_genes)),
         "sha256": sha256(CORE_PATH), "notes": "Frozen before v2 analysis"},
    ])
    write_tsv(manifest, ROOT / "metadata" / "input_manifest_v2.tsv")
    software = pd.DataFrame([{
        "tool": "Python", "version": sys.version.replace("\n", " "), "source": "local_miniconda", "environment": sys.executable,
        "executable": sys.executable, "checked_at": stamp(), "notes": f"platform={platform.platform()}"},
        {"tool": "h5py", "version": h5py.__version__, "source": "Python package", "environment": sys.executable,
         "executable": sys.executable, "checked_at": stamp(), "notes": "HDF5 read-only audit"},
        {"tool": "pandas", "version": pd.__version__, "source": "Python package", "environment": sys.executable,
         "executable": sys.executable, "checked_at": stamp(), "notes": "tabular audit outputs"},
    ])
    write_tsv(software, ROOT / "metadata" / "software_manifest_v2.tsv")

    report = "\n".join([
        "# Phase A-B v2 data-contract audit", "", f"- Analysis ID: `{VERSION}`", f"- Verified at: `{stamp()}`",
        "- Raw H5AD was inspected read-only; `raw/X` stored values are non-negative integers and retain the cells-by-genes orientation.",
        "- Donor mapping was rechecked: 7 Controlled and 9 Natural donors; a donor maps to a single CycleType.",
        "- Public fine cell-type labels were audited with a fixed marker panel that is disjoint from ESR1, PGR, the frozen WOI set and the stable core; no labels were changed.",
        "- The H5AD does not encode individual age, dose/regimen, or timing anchor. Those fields remain unavailable rather than inferred.",
        "- Existing broad-epithelial v1 outputs are preserved for audit only. v2 must use a nonciliated epithelial main compartment.",
    ])
    (ROOT / "reports" / "phase_A_B_contract_v2.md").write_text(report + "\n", encoding="utf-8")
    stages = pd.DataFrame([
        {"stage": "A_baseline_and_protection", "status": "complete", "completed_at": stamp(), "inputs": "protected submission assets; v1 outputs", "outputs": "metadata/protected_assets_v2.tsv;metadata/v1_reusable_artifacts_v2.tsv", "qa": "hash inventory created; no protected asset mutation", "notes": "v1 PASS labels were not inherited"},
        {"stage": "B_semantic_protocol_annotation_audit", "status": "complete", "completed_at": stamp(), "inputs": "raw PC/NC H5AD and frozen definitions", "outputs": "metadata/PC_NC_*_v2;reports/phase_A_B_contract_v2.md", "qa": "count semantics, donor mapping, marker overlap and gene coverage checked", "notes": "Protocol fields not encoded in H5AD remain NA"},
        {"stage": "C_primary_correction", "status": "pending", "completed_at": "NA", "inputs": "pending validated pseudobulk", "outputs": "NA", "qa": "NA", "notes": "must use nonciliated epithelium and corrected families"},
        {"stage": "D_robustness", "status": "pending", "completed_at": "NA", "inputs": "pending C", "outputs": "NA", "qa": "NA", "notes": "LODO, thresholds, compartment and LOGO"},
        {"stage": "E_cross_context", "status": "pending", "completed_at": "NA", "inputs": "pending C/D and IUA branch verification", "outputs": "NA", "qa": "NA", "notes": "conditional only"},
        {"stage": "F_resource_preflight", "status": "pending", "completed_at": "NA", "inputs": "public metadata", "outputs": "NA", "qa": "NA", "notes": "access and independence boundary"},
        {"stage": "G_figures", "status": "pending", "completed_at": "NA", "inputs": "validated result tables", "outputs": "NA", "qa": "NA", "notes": "v2 figures only if evidence warrants"},
        {"stage": "H_decision_report", "status": "pending", "completed_at": "NA", "inputs": "all prior phases", "outputs": "NA", "qa": "NA", "notes": "no manuscript drafting"},
    ])
    write_tsv(stages, ROOT / "checkpoints" / "stage_status_v2.tsv")
    print(json.dumps({"analysis_id": VERSION, "status": "A_B_complete", "h5ad_sha256": sha256(H5AD)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
