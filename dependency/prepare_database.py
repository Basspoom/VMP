import argparse, os, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional

DEP_DIR = Path.cwd()
DB_DIR = DEP_DIR / "databases"
TOOL_DIR = DEP_DIR / "tools"
MODEL_DIR = DEP_DIR / "models"

def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> None:
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
    s = archive_path.suffix.lower()
    if s in [".gz", ".tgz", ".tar"]:
        run_cmd(["tar", "-xzf", archive_path.name], cwd=target_root)
    elif s == ".zip":
        run_cmd(["unzip", "-q", "-o", archive_path.name], cwd=target_root)
    else:
        raise ValueError(f"Unsupported archive format: {archive_path}")
    after = snapshot_dirs(target_root)
    return [target_root / d for d in sorted(list(after - before))]

def list_children(p: Path) -> List[Path]:
    return [c for c in p.iterdir() if c.name not in (".DS_Store",)]

def flatten_if_wrapped(final_dir: Path, max_depth: int = 6) -> None:
    depth = 0
    while depth < max_depth:
        depth += 1
        children = list_children(final_dir)
        if not children:
            return
        if any(ch.is_file() for ch in children):
            return
        if len(children) != 1 or not children[0].is_dir():
            return
        inner = children[0]
        for c in list_children(inner):
            _safe_move(c, final_dir)
        try:
            inner.rmdir()
        except Exception:
            return

def _safe_remove(p: Path) -> None:
    try:
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        else:
            p.unlink(missing_ok=True)
    except FileNotFoundError:
        pass

def _safe_move(src: Path, dst_dir: Path) -> None:
    """Move src into dst_dir, overwriting if target exists."""
    ensure_dir(dst_dir)
    dest = dst_dir / src.name
    if dest.exists() or dest.is_symlink():
        _safe_remove(dest)
    shutil.move(str(src), str(dst_dir))

def move_dirs_into(final_dir: Path, new_dirs: List[Path]) -> None:
    ensure_dir(final_dir)
    for d in new_dirs:
        if d.is_dir():
            for c in list_children(d):
                _safe_move(c, final_dir)
            try:
                d.rmdir()
            except Exception:
                pass
        else:
            _safe_move(d, final_dir)
    flatten_if_wrapped(final_dir)

def require_binaries(bins: List[str]) -> None:
    missing = [b for b in bins if shutil.which(b) is None]
    if missing:
        raise EnvironmentError(f"Missing system binaries: {missing}")

def assert_cwd_structure():
    ensure_dir(DB_DIR); ensure_dir(TOOL_DIR); ensure_dir(MODEL_DIR)

def aria2_cmd(url: str, out_name: str) -> List[str]:
    return [
        "aria2c","-c","-x","16","-s","16","-j","16","-k","100M",
        "--auto-file-renaming=false","--summary-interval=0","--console-log-level=warn",
        "-o", out_name, url,
    ]

def wget_cmd(url: str, out_name: str) -> List[str]:
    return ["wget","-q","--show-progress","--progress=bar:force:noscroll","-c",url,"-O",out_name]

def pick_downloader(use_aria2: bool) -> Callable[[str, str], List[str]]:
    if use_aria2 and shutil.which("aria2c") is not None:
        return lambda url,out: aria2_cmd(url,out)
    return lambda url,out: wget_cmd(url,out)

TASKS: Dict[str, Dict[str, str]] = {
    "VirSorter2": {"kind":"db","url":"https://osf.io/v46sc/download","archive":"VirSorter2.tgz","target_dir":"VirSorter2_db"},
    "CheckV": {"kind":"db","url":"https://zenodo.org/records/14033148/files/checkv1.5.tar.gz","archive":"checkv1.5.tar.gz","target_dir":"CheckV_db"},
    "geNomad": {"kind":"db","url":"https://portal.nersc.gov/genomad/__data__/genomad_db_v1.9.tar.gz","archive":"genomad_db_v1.7.tar.gz","target_dir":"geNomad_db"},
    # "geNomad": {"kind":"db","url":"https://zenodo.org/records/14033148/files/genomad_db_v1.7.tar.gz","archive":"genomad_db_v1.7.tar.gz","target_dir":"geNomad_db"},
    "VIBRANT": {"kind":"db","url":"https://zenodo.org/records/14033148/files/vibrant.tar.gz","archive":"vibrant.tar.gz","target_dir":"VIBRANT_db"},
    "Kaiju": {"kind":"db","url":"https://zenodo.org/records/14033148/files/kaiju.tar.gz","archive":"kaiju.tar.gz","target_dir":"Kaiju_db"},
    "DeepVirFinder": {"kind":"tool","url":"https://github.com/jessieren/DeepVirFinder/archive/refs/heads/master.zip","archive":"DeepVirFinder.zip","target_dir":"DeepVirFinder"},
    "VPAC": {"kind":"model","url":"https://zenodo.org/records/17008739/files/VPAC_models.zip","archive":"VPAC_models.zip","target_dir":"VPAC_models"},

    "Homo_sapiens": {"kind":"db","url":"http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz","archive":"Homo_sapiens_UCSC_hg38.tar.gz","target_dir":"Homo_sapiens"},
    "checkM_db": {"kind":"db","url":"https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz","archive":"checkm_data_2015_01_16.tar.gz","target_dir":"CheckM_db"},
    "genomescope2.0": {"kind":"tool","url":"https://github.com/tbenavi1/genomescope2.0/archive/refs/heads/master.zip","archive":"genomescope2.0.zip","target_dir":"genomescope2.0"},
    "Evo": {"kind":"model","url":"https://zenodo.org/records/14033148/files/evo-model.zip","archive":"evo-model.zip","target_dir":"Evo"},
}

CATEGORY_DIR = {"db": DB_DIR, "tool": TOOL_DIR, "model": MODEL_DIR}

def prepare_hg38(root: Path, archive_path: Path, final_dir: Path) -> None:
    tmp_dirs = extract_and_get_new_dir(archive_path, root)
    base = None
    for d in tmp_dirs:
        cand = d / "UCSC" / "hg38" / "Sequence" / "Chromosomes"
        if cand.is_dir():
            base = cand; break
        for dd in d.rglob("Chromosomes"):
            if dd.is_dir():
                base = dd; break
        if base: break
    if not base:
        raise RuntimeError("hg38 Chromosomes dir not found after extraction")

    chr_names = [f"chr{i}.fa" for i in range(1,23)] + ["chrX.fa","chrY.fa"]
    ensure_dir(final_dir)
    out_fa = final_dir / "hg38.fa"
    if out_fa.exists() or out_fa.is_symlink():
        _safe_remove(out_fa)
    with out_fa.open("wb") as w:
        for nm in chr_names:
            p = base / nm
            if not p.exists():
                raise RuntimeError(f"Missing {nm} in hg38 package")
            with p.open("rb") as r:
                shutil.copyfileobj(r, w)
    for d in tmp_dirs:
        try:
            shutil.rmtree(d)
        except Exception:
            pass

def prepare_one(name: str, clean: bool, dl_cmd_builder: Callable[[str, str], List[str]]) -> Tuple[str, Path]:
    t = TASKS[name]
    kind, url, archive, target_dirname = t["kind"], t["url"], t["archive"], t["target_dir"]
    root = CATEGORY_DIR[kind]
    ensure_dir(root)
    final_dir = root / target_dirname

    if name == "Evo":
        archive_path = final_dir / archive
    else:
        archive_path = root / archive

    print(f"\n=== Preparing {name} ({kind}) ===")
    if final_dir.exists():
        print(f"  Overwrite: {final_dir}")

    if name == "VirSorter2":
        tmp = root / "download"
        if tmp.exists():
            tmp.unlink()
        run_cmd(wget_cmd(url, "download"), cwd=root)
        if archive_path.exists():
            archive_path.unlink()
        tmp.rename(archive_path)
    elif name == "Evo":
        ensure_dir(final_dir)
        run_cmd(dl_cmd_builder(url, archive), cwd=final_dir)
    else:
        run_cmd(dl_cmd_builder(url, archive), cwd=root)

    if name == "Homo_sapiens":
        prepare_hg38(root, archive_path, final_dir)
    elif name == "Evo":
        extract_and_get_new_dir(archive_path, final_dir)
        flatten_if_wrapped(final_dir)
    else:
        new_dirs = extract_and_get_new_dir(archive_path, root)
        if name == "VirSorter2":
            src = root / "db"
            if not src.exists() and len(new_dirs) == 1:
                src = new_dirs[0]
            if src.exists() and src.is_dir():
                ensure_dir(final_dir)
                for c in list_children(src):
                    _safe_move(c, final_dir)
                try:
                    src.rmdir()
                except Exception:
                    pass
            else:
                move_dirs_into(final_dir, new_dirs)
            flatten_if_wrapped(final_dir)
        elif name in {"checkM_db", "genomescope2.0", "Evo", "DeepVirFinder", "VPAC", "geNomad", "VIBRANT", "Kaiju", "CheckV"}:
            move_dirs_into(final_dir, new_dirs)
        else:
            move_dirs_into(final_dir, new_dirs)

    if clean and archive_path.exists():
        try:
            archive_path.unlink()
        except Exception:
            pass

    print(f"  Done -> {final_dir.resolve()}")
    return (name, final_dir.resolve())

def main():
    parser = argparse.ArgumentParser(
        description="Prepare deps in VMP/dependency",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-db", metavar="NAMES", required=True,
        help=(
            "Comma list or 'all'. Available:\n"
            "  VirSorter2, CheckV, geNomad, VIBRANT, Kaiju, DeepVirFinder, VPAC,\n"
            "  Homo_sapiens, checkM_db, genomescope2.0, Evo"
        )
    )
    parser.add_argument("-clean", action="store_true", help="Remove archives after extraction.")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("-aria2", dest="aria2", action="store_true", help="Prefer aria2c")
    g.add_argument("-no-aria2", dest="aria2", action="store_false", help="Disable aria2c")
    parser.set_defaults(aria2=None)
    args = parser.parse_args()

    assert_cwd_structure()
    require_binaries(["tar","unzip"])

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
                print(f"Unknown component: {x}"); sys.exit(1)

    if args.aria2 is None:
        use_aria2 = shutil.which("aria2c") is not None
    else:
        use_aria2 = args.aria2 and (shutil.which("aria2c") is not None)
        if args.aria2 and not use_aria2:
            print("  aria2c requested but not found; fallback to wget.")

    dl_cmd_builder = pick_downloader(use_aria2)
    if not use_aria2:
        require_binaries(["wget"])
    if "VirSorter2" in sel or "Homo_sapiens" in sel:
        require_binaries(["wget"])

    results = []
    for name in sel:
        results.append(prepare_one(name, clean=args.clean, dl_cmd_builder=dl_cmd_builder))

    print("\n=== Summary ===")
    for nm, p in sorted(results):
        print(f"{nm}: '{p}'")
    print("\nUpdate VMP/config.yml if needed.")

if __name__ == "__main__":
    main()
