"""Shared path resolution for the VMH_NP pipeline.

The data root is taken from the environment variable VMH_NP_DATA if set, otherwise from
config/paths.json ("data_root"). Downloaded Allen files live in <data_root>/abc_atlas,
mirroring the relative paths in config/data_manifest.csv.
"""
import csv
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config"
RESULTS = REPO / "results"


def data_root() -> Path:
    env = os.environ.get("VMH_NP_DATA")
    if env:
        return Path(env).expanduser()
    with open(CONFIG / "paths.json") as fh:
        return Path(json.load(fh)["data_root"]).expanduser()


def abc_root() -> Path:
    return data_root() / "abc_atlas"


def manifest_rows():
    with open(CONFIG / "data_manifest.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def abc_file(directory: str, name: str) -> Path:
    """Local path of a manifest entry, e.g. abc_file('WMB-10Xv3', 'WMB-10Xv3-HY/log2')."""
    for r in manifest_rows():
        if r["directory"] == directory and r["name"] == name:
            return abc_root() / r["relative_path"]
    raise KeyError(f"{directory}/{name} is not in config/data_manifest.csv")
