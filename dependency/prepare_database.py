#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

DEP_DIR = Path.cwd()
DB_DIR = DEP_DIR / "databases"
TOOL_DIR = DEP_DIR / "tools"
MODEL_DIR = DEP_DIR / "models"

def run_cmd(cmd: List[str], cwd: Path | None = None) -> None:
    # Concise echo so users know what's happening
    print("  >", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def snapshot_dirs(root: Path) -> set:
    return {d for d in os.listdir(root) if (root / d).is_dir()}

def extract_and_get_new_dir(archive_path: Path, target_root: Path) -> List[Path]:
    before = snapshot_dirs(target_root)
    # Quiet extraction: no verbose file listing
    if archive_path.suffix in [".gz", ".tgz", ".tar"]:
        run_cmd(["tar", "-xzf", archive_path.name], cwd=target_root)
    elif archive_path.suffix == ".zip":
        run_cmd(["unzip", "-q", "-o", archive_path.name], cwd=target_root)
    else:
        raise ValueError(f"Unsupported archive format: {archive_path}")
    after = snapshot_dirs(target_root)
    return [target_root / d for d in sorted(list(after - before))]

def safe_rename_to(target_dir: Path, desired_name: str) -> Path:
    final_path = target_dir.parent / desired_name
    if final_path.exists():
        return final_path
    try:
        target_dir.rename(final_path)
        return final_path
    except Exception:
        shutil.move(str(target_dir), str(final_path))
        return final_path

def require_binaries(bins: List[str]) -> None:
    missing = [b for b in bins if shutil.which(b) is None]
    if missing:
        raise EnvironmentError(f"Missing system binaries: {missing}")

def assert_cwd_structure():
    ensure_dir(DB_DIR); ensure_dir(TOOL_DIR); ensure_dir(MODEL_DIR)

TASKS: Dict[str, Dict[str, str]] = {
    "VirSorter2": {"kind": "db", "url": "https://osf.io/v46sc/download", "archive": "VirSorter2.tgz", "target_dir": "VirSorter2_db"},
    "CheckV": {"kind": "db", "url": "https://zenodo.org/records/14033148/files/checkv1.5.tar.gz", "archive": "checkv1.5.tar.gz", "target_dir": "CheckV_db"},
    "geNomad": {"kind": "db", "url": "https://zenodo.org/records/14033148/files/genomad_db_v1.7.tar.gz", "archive": "genomad_db_v1.7.tar.gz", "target_dir": "geNomad_db"},
    "VIBRANT": {"kind": "db", "url": "https://zenodo.org/records/14033148/files/vibrant.tar.gz", "archive": "vibrant.tar.gz", "target_dir": "VIBRANT_db"},
    "Kaiju": {"kind": "db", "url": "https://zenodo.org/records/14033148/files/kaiju.tar.gz", "archive": "kaiju.tar.gz", "target_dir": "Kaiju_db"},
    "DeepVirFinder": {"kind": "tool", "url": "https://github.com/jessieren/DeepVirFinder/archive/refs/heads/master.zip", "archive": "DeepVirFinder.zip", "target_dir": "DeepVirFinder"},
    "VPAC": {"kind": "model", "url": "https://zenodo.org/records/17008739/files/VPAC_models.zip", "archive": "VPAC_models.zip", "target_dir": "VPAC_models"},
}

CATEGORY_DIR = {"db": DB_DIR, "tool": TOOL_DIR, "model": MODEL_DIR}

def aria2_cmd(url: str, out_name: str) -> List[str]:
    # aria2c concise progress; resume; parallel; no auto-renaming
    return [
        "aria2c",
        "-c",
        "-x", "16",
        "-s", "16",
        "-j", "16",
        "-k", "100M",
        "--auto-file-renaming=false",
        "--summary-interval=0",
        "--console-log-level=warn",
        "-o", out_name,
        url,
    ]

def wget_cmd(url: str, out_name: str) -> List[str]:
    # Show only progress bar; suppress other noise
    return ["wget", "-q", "--show-progress", "--progress=bar:force:noscroll", "-c", url, "-O", out_name]

def pick_downloader(use_aria2: bool) -> callable:
    if use_aria2 and shutil.which("aria2c") is not None:
        def _cmd(url: str, out_name: str) -> List[str]:
            return aria2_cmd(url, out_name)
        return _cmd
    # fallback to wget
    def _cmd(url: str, out_name: str) -> List[str]:
        return wget_cmd(url, out_name)
    return _cmd

def prepare_one(name: str, clean: bool, dl_cmd_builder: callable) -> Tuple[str, Path]:
    task = TASKS[name]
    kind, url, archive, target_dirname = task["kind"], task["url"], task["archive"], task["target_dir"]
    root = CATEGORY_DIR[kind]
    ensure_dir(root)
    archive_path = root / archive
    final_dir = root / target_dirname

    print(f"\n=== Preparing {name} ({kind}) ===")

    if final_dir.exists():
        print(f"  Already exists: {final_dir}")
        return (name, final_dir.resolve())

    if name == "VirSorter2":
        tmp = root / "download"
        if tmp.exists():
            tmp.unlink()
        run_cmd(dl_cmd_builder(url, "download"), cwd=root)
        if archive_path.exists():
            archive_path.unlink()
        tmp.rename(archive_path)
    else:
        run_cmd(dl_cmd_builder(url, archive), cwd=root)

    new_dirs = extract_and_get_new_dir(archive_path, root)

    if name == "VirSorter2":
        src = root / "db"
        if not src.exists() and len(new_dirs) == 1:
            src = new_dirs[0]
        final_dir = safe_rename_to(src, target_dirname)
    else:
        if len(new_dirs) == 1:
            final_dir = safe_rename_to(new_dirs[0], target_dirname)
        else:
            ensure_dir(final_dir)
            for d in new_dirs:
                shutil.move(str(d), str(final_dir))

    if clean and archive_path.exists():
        try:
            archive_path.unlink()
        except:
            pass

    print(f"  Done -> {final_dir.resolve()}")
    return (name, final_dir.resolve())

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare dependencies required by VMP.\n"
            "This script must be executed inside the directory:  VMP/dependency\n\n"
            "Examples:\n"
            "  python prepare_database.py -db all\n"
            "  python prepare_database.py -db VirSorter2,CheckV,VPAC -clean\n"
            "  python prepare_database.py -db all --aria2\n"
            "  python prepare_database.py -db geNomad --no-aria2\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-db",
        metavar="NAMES",
        required=True,
        help=(
            "Comma-separated list or 'all'.\n"
            "Use a comma (no spaces required) to select multiple components, e.g. 'VirSorter2,CheckV,VPAC'.\n"
            "Available components:\n"
            "  VirSorter2, CheckV, geNomad, VIBRANT, Kaiju, DeepVirFinder, VPAC\n\n"
            "Examples:\n"
            "  -db all\n"
            "  -db VirSorter2,CheckV\n"
            "  -db DeepVirFinder,VPAC"
        )
    )

    parser.add_argument(
        "-clean",
        action="store_true",
        help=(
            "Remove downloaded archive files after extraction.\n"
            "Default: keep archive files."
        )
    )

    aria2_group = parser.add_mutually_exclusive_group()
    aria2_group.add_argument(
        "-aria2",
        dest="aria2",
        action="store_true",
        help="Prefer aria2c for downloading (falls back to wget if aria2c is not found)."
    )
    aria2_group.add_argument(
        "-no-aria2",
        dest="aria2",
        action="store_false",
        help="Do not use aria2c; use wget."
    )
    parser.set_defaults(aria2=None)  # Auto-detect if not specified

    args = parser.parse_args()

    assert_cwd_structure()
    # tar and unzip are always required; downloader depends on args
    require_binaries(["tar", "unzip"])
    if args.aria2 is None:
        # Auto: use aria2c if present
        use_aria2 = shutil.which("aria2c") is not None
    else:
        use_aria2 = args.aria2 and (shutil.which("aria2c") is not None)
        if args.aria2 and not use_aria2:
            print("  aria2c requested but not found; falling back to wget.")

    dl_cmd_builder = pick_downloader(use_aria2)
    if not use_aria2:
        require_binaries(["wget"])

    keys = list(TASKS.keys())
    if args.db.lower() == "all":
        sel = keys
    else:
        sel = []
        lut = {k.lower(): k for k in keys}
        for x in args.db.split(","):
            x = x.strip()
            if x.lower() in lut:
                sel.append(lut[x.lower()])
            else:
                print(f"Unknown component: {x}")
                sys.exit(1)

    results = []
    for name in sel:
        results.append(prepare_one(name, clean=args.clean, dl_cmd_builder=dl_cmd_builder))

    print("\n=== Summary ===")
    for nm, p in sorted(results):
        print(f"{nm}: '{p}'")

    print("\nPlease update VMP/config.yml accordingly.")

if __name__ == "__main__":
    main()
