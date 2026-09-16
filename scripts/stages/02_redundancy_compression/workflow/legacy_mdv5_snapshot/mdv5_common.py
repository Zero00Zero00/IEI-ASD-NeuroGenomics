#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import importlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml


class MDV5Error(RuntimeError):
    pass


def load_config(path: str | os.PathLike[str]) -> dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        raise MDV5Error(f"Config missing: {p}")
    with p.open("rt", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise MDV5Error("YAML root must be a mapping")
    return cfg


def root_from_cfg(cfg: dict[str, Any]) -> Path:
    return Path(os.environ.get("ROOT", cfg["project_root"])).expanduser().resolve()


def resolve(root: Path, value: str | os.PathLike[str]) -> Path:
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (root / p).resolve()


def outdir(cfg: dict[str, Any], root: Path) -> Path:
    p = resolve(root, cfg["output_dir"])
    p.mkdir(parents=True, exist_ok=True)
    return p


def output_path(cfg: dict[str, Any], root: Path, key: str) -> Path:
    return outdir(cfg, root) / cfg["outputs"][key]


def read_table(path: Path, dtype: Any = None) -> pd.DataFrame:
    if not path.is_file():
        raise MDV5Error(f"Required table missing: {path}")
    compression = "gzip" if path.suffix == ".gz" else "infer"
    return pd.read_csv(path, sep="\t", dtype=dtype, compression=compression, low_memory=False)


def write_table(df: pd.DataFrame, path: Path, gzip_output: bool | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzip_output is None:
        gzip_output = path.suffix == ".gz"
    tmp = path.with_name(path.name + ".tmp")
    if gzip_output:
        df.to_csv(tmp, sep="\t", index=False, compression="gzip")
    else:
        df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmpname = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wt", encoding="utf-8") as handle:
            handle.write(text)
            if not text.endswith("\n"):
                handle.write("\n")
        os.replace(tmpname, path)
    finally:
        if os.path.exists(tmpname):
            os.unlink(tmpname)


def atomic_json(path: Path, obj: Any) -> None:
    atomic_text(path, json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False))


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def normalize_source(value: Any) -> str:
    s = str(value).strip()
    k = re.sub(r"[^a-z0-9]", "", s.lower())
    if k in {"go", "gobp", "geneontologybiologicalprocess", "biologicalprocess"}:
        return "GO_BP"
    if k in {"reactome", "react"}:
        return "Reactome"
    return s


def pick_col(df: pd.DataFrame, aliases: Iterable[str], label: str, required: bool = True) -> str | None:
    direct = {str(c): str(c) for c in df.columns}
    canon = {re.sub(r"[^a-z0-9]", "", str(c).lower()): str(c) for c in df.columns}
    for a in aliases:
        if a in direct:
            return direct[a]
        k = re.sub(r"[^a-z0-9]", "", a.lower())
        if k in canon:
            return canon[k]
    if required:
        raise MDV5Error(f"Cannot identify column for {label}; aliases={list(aliases)}; columns={list(df.columns)}")
    return None


def canonicalize_pathway_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    c_source = pick_col(out, ["source", "pathway_source", "database"], "source")
    c_id = pick_col(out, ["pathway_id", "term_id", "ID", "id"], "pathway_id")
    c_name = pick_col(out, ["pathway_name", "term_name", "name", "Pathway"], "pathway_name")
    rename = {c_source: "source", c_id: "pathway_id", c_name: "pathway_name"}
    out = out.rename(columns=rename)
    out["source"] = out["source"].map(normalize_source)
    out["pathway_id"] = out["pathway_id"].astype(str).str.strip()
    out["pathway_name"] = out["pathway_name"].astype(str).str.strip()
    out["pathway_key"] = out["source"] + "::" + out["pathway_id"]
    return out


def canonicalize_membership(df: pd.DataFrame) -> pd.DataFrame:
    out = canonicalize_pathway_columns(df)
    c_hid = pick_col(out, ["HGNC_id", "hgnc_id", "HGNC ID"], "HGNC_id")
    c_sym = pick_col(out, ["approved_symbol", "HGNC_symbol", "symbol", "gene"], "approved_symbol")
    out = out.rename(columns={c_hid: "HGNC_id", c_sym: "approved_symbol"})
    out["HGNC_id"] = out["HGNC_id"].astype(str).str.strip()
    out["approved_symbol"] = out["approved_symbol"].astype(str).str.strip()
    return out


def canonicalize_universe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    c_hid = pick_col(out, ["HGNC_id", "hgnc_id", "HGNC ID"], "HGNC_id")
    c_sym = pick_col(out, ["HGNC_symbol", "approved_symbol", "symbol", "gene"], "HGNC_symbol")
    c_core = pick_col(out, ["CoreSeed_flag", "coreseed_flag", "CoreSeed"], "CoreSeed_flag")
    out = out.rename(columns={c_hid: "HGNC_id", c_sym: "HGNC_symbol", c_core: "CoreSeed_flag"})
    out["HGNC_id"] = out["HGNC_id"].astype(str).str.strip()
    out["HGNC_symbol"] = out["HGNC_symbol"].astype(str).str.strip()
    return out


def as_bool(series: pd.Series) -> pd.Series:
    def conv(x: Any) -> bool:
        if pd.isna(x):
            return False
        s = str(x).strip().lower()
        return s in {"1", "true", "t", "yes", "y", "pass", "coreseed"}
    return series.map(conv)


def ensure_numeric(df: pd.DataFrame, aliases: Iterable[str], target: str, required: bool = True) -> pd.DataFrame:
    c = pick_col(df, aliases, target, required=required)
    if c is None:
        return df
    if c != target:
        df = df.rename(columns={c: target})
    df[target] = pd.to_numeric(df[target], errors="coerce")
    if required and df[target].isna().any():
        bad = int(df[target].isna().sum())
        raise MDV5Error(f"Column {target} contains {bad} non-numeric/missing values")
    return df


def sentinel_passes(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    txt = path.read_text(encoding="utf-8", errors="replace")
    pats = [
        r"(?im)^\s*(?:status|technical_status|data_provenance_status)\s*[:=]\s*PASS\s*$",
        r"(?im)^\s*STATUS\s*:\s*PASS\s*$",
        r"(?im)^\s*PASS\s*$",
        r"(?im)\bstatus\s*=\s*PASS\b",
    ]
    return any(re.search(p, txt) for p in pats) or "PASS" in txt.upper()


def parse_checksum_manifest(manifest: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with manifest.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"^([0-9a-fA-F]{64})\s+[* ]?(.+?)\s*$", line)
            if not m:
                raise MDV5Error(f"Unparseable checksum line in {manifest}: {line}")
            rows.append((m.group(1).lower(), m.group(2)))
    if not rows:
        raise MDV5Error(f"Checksum manifest empty: {manifest}")
    return rows


def resolve_checksum_target(root: Path, manifest: Path, rel: str) -> Path:
    p = Path(rel)
    candidates = []
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.extend([root / p, manifest.parent / p])
    for c in candidates:
        if c.is_file():
            return c.resolve()
    raise MDV5Error(f"Checksum target not found for '{rel}' from manifest {manifest}")


def verify_checksum_manifest(root: Path, manifest: Path) -> tuple[int, list[str]]:
    if not manifest.is_file():
        raise MDV5Error(f"Checksum manifest missing: {manifest}")
    failures: list[str] = []
    rows = parse_checksum_manifest(manifest)
    for expected, rel in rows:
        try:
            target = resolve_checksum_target(root, manifest, rel)
            actual = sha256_file(target)
            if actual.lower() != expected:
                failures.append(f"MISMATCH {rel} expected={expected} actual={actual}")
        except Exception as exc:
            failures.append(f"ERROR {rel}: {exc}")
    return len(rows), failures


def current_input_hashes(cfg: dict[str, Any], root: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for section in ("upstream", "optional_upstream"):
        for key, value in cfg.get(section, {}).items():
            p = resolve(root, value)
            if p.is_file():
                result[f"{section}.{key}"] = {
                    "path": str(p),
                    "size": p.stat().st_size,
                    "sha256": sha256_file(p),
                }
    return result


def dependency_versions(names: Iterable[str]) -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in names:
        module_name = "yaml" if name == "yaml" else name
        mod = importlib.import_module(module_name)
        versions[name] = str(getattr(mod, "__version__", "unknown"))
    return versions


def static_forbidden_scan(cfg: dict[str, Any], root: Path) -> list[str]:
    lock = cfg.get("gwas_blind_lock", {})
    literals = [str(x) for x in lock.get("static_forbidden_path_literals", [])]
    hits: list[str] = []
    for rel in lock.get("static_scan_files", []):
        p = resolve(root, rel)
        if not p.is_file():
            hits.append(f"MISSING_SCAN_TARGET:{rel}")
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        for lit in literals:
            if lit.lower() in txt.lower():
                hits.append(f"FORBIDDEN_LITERAL:{rel}:{lit}")
    return hits


def validate_configured_inputs_gwas_blind(cfg: dict[str, Any]) -> list[str]:
    lock = cfg.get("gwas_blind_lock", {})
    allowed = [str(x).lower() for x in lock.get("allowed_input_prefixes", [])]
    forbidden = [str(x).lower() for x in lock.get("forbidden_input_markers", [])]
    hits: list[str] = []
    for key, value in cfg.get("upstream", {}).items():
        v = str(value).replace("\\", "/")
        vl = v.lower()
        if allowed and not any(vl.startswith(a) for a in allowed):
            hits.append(f"INPUT_OUTSIDE_ALLOWLIST:{key}:{value}")
        if any(tok in vl for tok in forbidden):
            hits.append(f"FORBIDDEN_INPUT:{key}:{value}")
    return hits


def write_pass(path: Path, fields: dict[str, Any]) -> None:
    lines = [f"{k}={v}" for k, v in fields.items()]
    atomic_text(path, "\n".join(lines))


def path_key_sort(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["source", "pathway_id"], kind="mergesort").reset_index(drop=True)
