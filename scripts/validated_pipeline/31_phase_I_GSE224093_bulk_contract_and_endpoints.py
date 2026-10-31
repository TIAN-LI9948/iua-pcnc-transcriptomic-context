"""Phase I: contract-first, tissue-level external IUA endpoint analysis.

This script deliberately does not treat GSE224093 as a menstrual-phase,
epithelial, or WOI replication.  It analyses four previously frozen endpoints
in the provider's bulk endometrial FPKM matrix as a separate exploratory
family.  It never alters submission assets, primary PC/NC results, or raw
files, and refuses to overwrite an existing Phase I delivery.
"""
from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import os
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


_context_root = os.environ.get("IUA_CONTEXT_ROOT")
if not _context_root:
    raise RuntimeError("Set IUA_CONTEXT_ROOT to the audited context-validation workspace before running this provenance script.")
ROOT = Path(_context_root).expanduser().resolve()
PROJECT = ROOT.parents[1]
RAW = ROOT / "raw" / "GSE224093"
METADATA = ROOT / "metadata"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"
CHECKPOINTS = ROOT / "checkpoints"
SOFT = RAW / "GSE224093_family.soft.gz"
FPKM = RAW / "GSE224093_gene_FPKM.txt.gz"
WOI25 = PROJECT / "results" / "reference_WOI_program_epithelial.tsv"
STABLE = PROJECT / "results" / "reference_epithelial_WOI_program_stability_core_cell_specificity.tsv"
SEED = 202609192
BOOTSTRAPS = 20_000

OUTPUTS = [
    METADATA / "GSE224093_input_contract_v1.tsv",
    METADATA / "GSE224093_sample_ledger_v1.tsv",
    METADATA / "GSE224093_software_manifest_v1.tsv",
    METADATA / "GSE224093_phase_I_qa_v1.json",
    RESULTS / "GSE224093_gene_mapping_v1.tsv",
    RESULTS / "GSE224093_endpoint_scores_bulk_v1.tsv",
    RESULTS / "GSE224093_endpoint_tests_bulk_v1.tsv",
    REPORTS / "GSE224093_bulk_context_report_v1.md",
    CHECKPOINTS / "phase_I_GSE224093_bulk_v1.tsv",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: list[float]) -> list[float]:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    ranked = arr[order] * len(arr) / np.arange(1, len(arr) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out.tolist()


def parse_soft(path: Path) -> pd.DataFrame:
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    records: list[dict[str, str]] = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        gsm = lines[0].strip()
        description = next((x.split(" = ", 1)[1] for x in lines if x.startswith("!Sample_description = ")), "")
        title = next((x.split(" = ", 1)[1] for x in lines if x.startswith("!Sample_title = ")), "")
        source = next((x.split(" = ", 1)[1] for x in lines if x.startswith("!Sample_source_name_ch1 = ")), "")
        match = re.match(r"^(Con|Ba)([1-7])\s+from", description)
        if not match:
            raise ValueError(f"Cannot derive locked matrix ID/group from SOFT description for {gsm}: {description!r}")
        matrix_id = f"{match.group(1)}{match.group(2)}"
        group = "normal_control" if match.group(1) == "Con" else "severe_IUA"
        records.append(
            {
                "GSM": gsm,
                "matrix_sample_id": matrix_id,
                "group": group,
                "sample_title": title,
                "sample_description": description,
                "source_name": source,
                "phase_or_cycle_timing": "not_reported_in_SOFT",
                "hormone_exposure": "not_reported_in_SOFT",
                "sample_unit": "bulk_endometrium",
            }
        )
    meta = pd.DataFrame(records).sort_values("matrix_sample_id", key=lambda x: x.str.replace(r"(Con|Ba)", lambda z: "0" if z.group(1) == "Con" else "1", regex=True)).reset_index(drop=True)
    if len(meta) != 14 or meta["GSM"].nunique() != 14 or meta["matrix_sample_id"].nunique() != 14:
        raise ValueError("SOFT metadata do not establish 14 unique samples.")
    if meta.groupby("group").size().to_dict() != {"normal_control": 7, "severe_IUA": 7}:
        raise ValueError("SOFT metadata do not establish the prespecified 7-versus-7 contrast.")
    return meta


def exact_two_sided(values: np.ndarray, severe_index: np.ndarray) -> tuple[float, int]:
    n = len(values)
    k = len(severe_index)
    obs = values[severe_index].mean() - np.delete(values, severe_index).mean()
    all_stats = []
    for selected in itertools.combinations(range(n), k):
        selected = np.asarray(selected, dtype=int)
        all_stats.append(values[selected].mean() - np.delete(values, selected).mean())
    all_stats = np.asarray(all_stats, dtype=float)
    p = np.mean(np.abs(all_stats) >= (abs(obs) - 1e-12))
    return float(p), int(len(all_stats))


def bootstrap_difference(values: np.ndarray, severe_index: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    severe = values[severe_index]
    control = np.delete(values, severe_index)
    draws = (
        rng.choice(severe, size=(BOOTSTRAPS, len(severe)), replace=True).mean(axis=1)
        - rng.choice(control, size=(BOOTSTRAPS, len(control)), replace=True).mean(axis=1)
    )
    return tuple(np.quantile(draws, [0.025, 0.975]).tolist())


def main() -> None:
    if any(path.exists() for path in OUTPUTS):
        present = [str(path) for path in OUTPUTS if path.exists()]
        raise RuntimeError(f"Refusing to overwrite existing Phase I outputs: {present}")
    for path in [SOFT, FPKM, WOI25, STABLE]:
        if not path.exists():
            raise FileNotFoundError(path)

    sample_ledger = parse_soft(SOFT)
    raw = pd.read_csv(FPKM, sep="\t", compression="gzip", dtype={"GeneName": "string"})
    if raw.columns[0] != "GeneName" or raw.shape != (34043, 15):
        raise ValueError(f"Unexpected FPKM matrix schema: {raw.shape}, {raw.columns[:3].tolist()}")
    matrix_samples = raw.columns[1:].tolist()
    if set(matrix_samples) != set(sample_ledger["matrix_sample_id"]):
        raise ValueError("FPKM sample columns and SOFT-derived sample ledger do not match exactly.")
    if raw["GeneName"].isna().any():
        raise ValueError("Missing GeneName values are not permitted.")
    numeric = raw.iloc[:, 1:].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy()).all() or (numeric.to_numpy() < 0).any():
        raise ValueError("FPKM matrix contains non-finite or negative values.")

    # Deterministic duplicate handling declared before endpoint calculation.
    raw["GeneName"] = raw["GeneName"].str.strip()
    duplicate_rows = int(raw["GeneName"].duplicated().sum())
    collapsed = raw.groupby("GeneName", sort=True, as_index=True)[matrix_samples].mean()
    log_fpkm = np.log2(collapsed + 1.0)

    woi = pd.read_csv(WOI25, sep="\t")
    stable = pd.read_csv(STABLE, sep="\t")
    frozen_woi = woi["gene"].astype(str).str.strip().tolist()
    stable_core = stable.loc[stable["epithelial_enriched_2fold_reference_rule"].astype(bool), "gene"].astype(str).str.strip().tolist()
    if len(frozen_woi) != 25 or len(stable_core) != 11:
        raise ValueError("Frozen endpoint definitions have changed from their contract sizes.")

    endpoint_definitions = [
        {"endpoint": "PGR_log2_FPKM", "family": "exploratory_external_4_endpoint", "kind": "gene", "genes": ["PGR"]},
        {"endpoint": "frozen_WOI25_zscore", "family": "exploratory_external_4_endpoint", "kind": "score", "genes": frozen_woi},
        {"endpoint": "ESR1_log2_FPKM", "family": "exploratory_external_4_endpoint", "kind": "gene", "genes": ["ESR1"]},
        {"endpoint": "stable_epithelial_core_zscore", "family": "exploratory_external_4_endpoint", "kind": "score", "genes": stable_core},
    ]
    gene_universe = set(log_fpkm.index)
    mapping_rows = []
    score_values: dict[str, np.ndarray] = {}
    for definition in endpoint_definitions:
        mapped = [gene for gene in definition["genes"] if gene in gene_universe]
        for gene in definition["genes"]:
            mapping_rows.append(
                {
                    "endpoint": definition["endpoint"],
                    "gene": gene,
                    "mapped": gene in gene_universe,
                    "aggregation": "mean_FPKM_across_duplicate_symbols_then_log2_FPKM_plus_1",
                    "endpoint_role": definition["kind"],
                }
            )
        if definition["kind"] == "gene":
            if mapped != definition["genes"]:
                raise ValueError(f"Required receptor missing: {definition['endpoint']}")
            score_values[definition["endpoint"]] = log_fpkm.loc[mapped[0], matrix_samples].to_numpy(dtype=float)
        else:
            if not mapped:
                raise ValueError(f"No genes mapped for {definition['endpoint']}")
            values = log_fpkm.loc[mapped, matrix_samples].to_numpy(dtype=float)
            z_values = (values - values.mean(axis=1, keepdims=True)) / values.std(axis=1, ddof=1, keepdims=True)
            if not np.isfinite(z_values).all():
                raise ValueError(f"Non-finite z-score in {definition['endpoint']}")
            score_values[definition["endpoint"]] = z_values.mean(axis=0)

    mapping = pd.DataFrame(mapping_rows)
    ordered_ledger = sample_ledger.set_index("matrix_sample_id").loc[matrix_samples].reset_index()
    severe_index = np.flatnonzero(ordered_ledger["group"].eq("severe_IUA").to_numpy())
    rng = np.random.default_rng(SEED)
    tests = []
    long_scores = []
    for definition in endpoint_definitions:
        endpoint = definition["endpoint"]
        values = score_values[endpoint]
        severe_mean = float(values[severe_index].mean())
        control_mean = float(np.delete(values, severe_index).mean())
        low, high = bootstrap_difference(values, severe_index, rng)
        p_value, assignments = exact_two_sided(values, severe_index)
        mapped_n = int(mapping.loc[(mapping.endpoint == endpoint) & mapping.mapped, "gene"].nunique())
        tests.append(
            {
                "dataset": "GSE224093",
                "analysis_scope": "exploratory_bulk_endometrium_not_WOI_or_epithelial_replication",
                "contrast": "severe_IUA_minus_normal_control",
                "endpoint": endpoint,
                "endpoint_family": definition["family"],
                "severe_IUA_n": 7,
                "normal_control_n": 7,
                "severe_IUA_mean": severe_mean,
                "normal_control_mean": control_mean,
                "difference_severe_IUA_minus_control": severe_mean - control_mean,
                "bootstrap_ci95_low": low,
                "bootstrap_ci95_high": high,
                "exact_permutation_p_two_sided": p_value,
                "exact_BH_q_within_external_family": np.nan,
                "permutation_assignments": assignments,
                "bootstrap_replicates": BOOTSTRAPS,
                "transform": "log2(FPKM+1); scores are mean per-gene z-scores across all 14 samples",
                "frozen_genes": len(definition["genes"]) if definition["kind"] == "score" else np.nan,
                "mapped_genes": mapped_n if definition["kind"] == "score" else np.nan,
                "cycle_phase": "not_reported",
                "cellular_compartment": "bulk_endometrium",
                "inferential_boundary": "tissue-level association; composition and timing confounded",
            }
        )
        for sample_id, group, value in zip(matrix_samples, ordered_ledger["group"], values):
            long_scores.append(
                {
                    "dataset": "GSE224093",
                    "matrix_sample_id": sample_id,
                    "group": group,
                    "endpoint": endpoint,
                    "value": float(value),
                    "value_scale": "log2_FPKM_plus_1" if definition["kind"] == "gene" else "mean_per_gene_zscore_log2_FPKM_plus_1",
                }
            )
    tests_df = pd.DataFrame(tests)
    tests_df["exact_BH_q_within_external_family"] = bh(tests_df["exact_permutation_p_two_sided"].tolist())

    qa = {
        "analysis_id": "GSE224093_BULK_ENDPOINTS_V1",
        "completed_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "status": "PASS_WITH_TISSUE_LEVEL_BOUNDARY",
        "input_hashes": {"SOFT": sha256(SOFT), "FPKM": sha256(FPKM), "frozen_WOI25": sha256(WOI25), "stable_epithelial_core": sha256(STABLE)},
        "matrix_shape_raw": [int(raw.shape[0]), int(raw.shape[1])],
        "matrix_shape_collapsed": [int(log_fpkm.shape[0]), int(log_fpkm.shape[1])],
        "duplicate_gene_rows_collapsed_by_mean_FPKM": duplicate_rows,
        "sample_group_counts": ordered_ledger.groupby("group").size().to_dict(),
        "phase_reported": False,
        "hormone_exposure_reported": False,
        "cell_type_resolved": False,
        "exact_assignments": 3432,
        "external_endpoint_family_size": 4,
        "formal_manuscript_written": False,
        "submission_assets_modified": False,
        "eligible_for_WOI_replication_claim": False,
        "eligible_for_epithelial_replication_claim": False,
    }

    input_contract = pd.DataFrame([
        {
            "input_id": "GSE224093_SOFT",
            "role": "external_bulk_metadata",
            "path": str(SOFT),
            "source": "NCBI GEO FTP",
            "accession": "GSE224093",
            "release": "public_2023-06-12",
            "download_date": datetime.now(timezone.utc).astimezone().date().isoformat(),
            "genome_build": "GRCh38",
            "organism": "Homo sapiens",
            "tissue": "bulk endometrium",
            "assay": "bulk RNA-seq metadata",
            "sample_size": "14 (7 severe IUA; 7 normal control)",
            "sha256": qa["input_hashes"]["SOFT"],
            "limitations": "cycle phase and hormone exposure not reported; tissue composition unresolved",
        },
        {
            "input_id": "GSE224093_FPKM",
            "role": "external_bulk_processed_expression",
            "path": str(FPKM),
            "source": "NCBI GEO supplementary file",
            "accession": "GSE224093",
            "release": "public_2023-06-12",
            "download_date": datetime.now(timezone.utc).astimezone().date().isoformat(),
            "genome_build": "GRCh38",
            "organism": "Homo sapiens",
            "tissue": "bulk endometrium",
            "assay": "RSEM FPKM",
            "sample_size": "14 (7 severe IUA; 7 normal control)",
            "sha256": qa["input_hashes"]["FPKM"],
            "limitations": "provider-normalized FPKM; not count-compatible; no cell-type resolution",
        },
    ])
    software = pd.DataFrame([
        {
            "tool": "Python",
            "version": sys.version.split()[0],
            "source": "C:/Users/zh/miniconda3/python.exe",
            "environment": platform.platform(),
            "purpose": "contract validation and exact endpoint permutation analysis",
        },
        {"tool": "numpy", "version": np.__version__, "source": "Python package", "environment": "active", "purpose": "numerical analysis"},
        {"tool": "pandas", "version": pd.__version__, "source": "Python package", "environment": "active", "purpose": "tabular input/output"},
    ])

    # Write only after all calculations and QA have succeeded.
    input_contract.to_csv(OUTPUTS[0], sep="\t", index=False)
    sample_ledger.to_csv(OUTPUTS[1], sep="\t", index=False)
    software.to_csv(OUTPUTS[2], sep="\t", index=False)
    mapping.to_csv(OUTPUTS[4], sep="\t", index=False)
    pd.DataFrame(long_scores).to_csv(OUTPUTS[5], sep="\t", index=False)
    tests_df.to_csv(OUTPUTS[6], sep="\t", index=False)
    OUTPUTS[3].write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    report = f"""# GSE224093 外部 bulk 队列端点审计（v1）

## 数据合同

- GEO SOFT 明确记录 7 名 severe IUA 与 7 名 normal control 的 bulk endometrium 样本；样本列与官方 FPKM 文件逐一匹配。
- 官方处理文件是 RSEM FPKM（GRCh38），不是原始计数；所有分析仅在 `log2(FPKM+1)` 尺度进行。
- `GeneName` 中 {duplicate_rows} 行重复符号按预注册的确定性规则合并：先在 FPKM 尺度取均值，后取对数。未按端点结果选择行。
- 月经周期相位、激素暴露及细胞组成没有在 SOFT 中报告。因此该队列不构成 WOI、上皮、受体活性、治疗反应或临床容受性的验证。

## 统计合同

- 只检验既有的 PGR、冻结 WOI25 程序、ESR1、稳定上皮核心四个端点；它们构成独立的 4 项探索性 BH 家族。
- 每个端点采用 7-vs-7 的完整标签枚举（3,432 种分配）得到双侧精确 P 值，20,000 次分层 donor/sample bootstrap 给出均值差 CI。
- 程序分数为每个基因在全部 14 个 bulk 样本内的 z 分数均值。该标准化仅支持同一 bulk 队列内的组间关联，不能与 snRNA-seq logCPM 效应量直接合并。

## 科学解释边界

结果文件中的差值均为 severe IUA 减 normal control。无论显著性如何，最多可描述为“独立 bulk 组织层面的疾病相关端点模式”；不能作为 IUA 与 HRT 分离、WOI 相位复制或非纤毛上皮复制的证据。
"""
    OUTPUTS[7].write_text(report, encoding="utf-8")
    checkpoint = pd.DataFrame([
        {
            "stage": "I_GSE224093_external_bulk_contract_and_endpoints",
            "status": "complete",
            "completed_at": qa["completed_at"],
            "input_fingerprint": qa["input_hashes"]["FPKM"],
            "parameter_fingerprint": "log2(FPKM+1); duplicate mean-FPKM; exact 7v7; 20k bootstrap; frozen four endpoint BH family",
            "outputs": ";".join(str(path.relative_to(ROOT)) for path in OUTPUTS[:-1]),
            "qa": qa["status"],
            "notes": "Exploratory bulk-tissue association only; not a WOI/epithelial replication.",
        }
    ])
    checkpoint.to_csv(OUTPUTS[8], sep="\t", index=False)
    print(json.dumps({"status": qa["status"], "endpoint_rows": len(tests_df), "sample_rows": len(sample_ledger), "duplicate_rows": duplicate_rows}, ensure_ascii=False))


if __name__ == "__main__":
    main()
