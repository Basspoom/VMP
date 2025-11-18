#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shlex
import argparse
import subprocess
from pathlib import Path
from typing import List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.markdown import Markdown

try:
    import yaml
except Exception:
    yaml = None

HERE = Path(__file__).resolve().parent

# Component scripts colocated with this file
QC_SCRIPT          = HERE / "QC.py"
ASM_SCRIPT         = HERE / "Assembly.py"
VPAC_SINGLE_SCRIPT = HERE / "VPAC_single.py"
VPAC_DUAL_SCRIPT   = HERE / "VPAC_dual.py"
CLUST_SCRIPT       = HERE / "Clustering.py"
BINNING_SCRIPT     = HERE / "Binning.py"

console = Console()

# Global log control (set in main from CLI flags)
QUIET: bool = True


def custom_help():
    console.print(Panel(
        "VMP End-to-End\n[bold]Pipeline:[/bold] QC → Assembly → VPAC → Clustering → Binning\n"
        "Config-driven (envs/tools/models) via a single config.yml.",
        border_style="blue", title="VMP End-to-End", title_align="left"
    ))

    console.print(Markdown(
        "**I/O overview**\n"
        "- [QC]         input: `-i` → output: `-o/Clean_reads/`\n"
        "- [Assembly]   input: `-o/Clean_reads/` → output: `-o/raw_contigs/`\n"
        "- [VPAC]       input: `-o/raw_contigs/` → output: `-o/Viral_contigs/<asm>-single/` (mVC, mVC_dual, mNVC, mNVC_dual)\n"
        "- [Clustering] input: VPAC viral sets → output: `-o/cluster_<tool>/` (unchanged)\n"
        "- [Binning]    input: [bold]VPAC non-viral[/bold] sets ONLY → "
        "`Viral_contigs/<asm>-single/(mNVC | mNVC_dual)` (dual preferred when both)\n"
        "                 output: `Viral_contigs/<asm>-single/Bin/`"
    ))

    console.print(Panel("Usage: VMP_end2end.py [OPTIONS] -i RAW_FASTQ_DIR -o OUT_DIR -cf CONFIG",
                        border_style="cyan", title="Global parameters", title_align="left"))

    g = Table(show_header=False, box=None, pad_edge=False)
    g.add_column("Flag", style="bold cyan", no_wrap=True)
    g.add_column("Description", style="white")
    g.add_row("-i, --input",  "Directory of raw paired FASTQ (one subdir per sample).")
    g.add_row("-o, --output", "Root output directory.")
    g.add_row("-cf, --config","Path to config.yml (envs/tools/models).")
    g.add_row("-t, --threads","Threads for tools that accept -t/--threads (default: 60).")
    g.add_row("--parallel",   "Parallel samples (default: 2).")
    g.add_row("--resume/--no-resume","Auto-resume (default: on).")
    g.add_row("--from-step",  "Force start from a specific step.")
    g.add_row("--quiet / --show-logs", "Hide sub-step logs (default: --quiet).")
    console.print(g)
    console.print()

    console.print(Panel("Assembly options", border_style="magenta"))
    a = Table(show_header=False, box=None, pad_edge=False)
    a.add_column("Flag", style="bold cyan", no_wrap=True)
    a.add_column("Description", style="white")
    a.add_row("--assembly", "spades | megahit (default: spades).")
    a.add_row("--spades-metaviral", "SPAdes metaviral mode.")
    console.print(a)
    console.print()

    console.print(Panel("VPAC options", border_style="green"))
    v = Table(show_header=False, box=None, pad_edge=False)
    v.add_column("Flag", style="bold cyan", no_wrap=True)
    v.add_column("Description", style="white")
    v.add_row("--vpac", "single | dual | both (default: both).")
    v.add_row("--device", "Device for VPAC-dual (e.g., cuda:0 / cpu).")
    v.add_row("--lmin/--lmax", "VPAC-single length bounds (defaults: 3000 / 100000).")
    console.print(v)
    console.print()

    console.print(Panel("Clustering options", border_style="yellow"))
    c = Table(show_header=False, box=None, pad_edge=False)
    c.add_column("Flag", style="bold cyan", no_wrap=True)
    c.add_column("Description", style="white")
    c.add_row("--cluster", "cd-hit_mmseq | skani_pyleiden (default: skani_pyleiden).")
    console.print(c)
    console.print()

    console.print(Panel("Binning options (non-viral only)", border_style="blue"))
    b = Table(show_header=False, box=None, pad_edge=False)
    b.add_column("Flag", style="bold cyan", no_wrap=True)
    b.add_column("Description", style="white")
    b.add_row("--no-binning", "Disable Binning stage (default: enabled).")
    b.add_row("--bin-minfasta-kb", "Minimum bin size (Kbp, default 200).")
    b.add_row("--bin-min-len", "Minimal contig length (default 500).")
    b.add_row("--bin-min-completeness", "Minimum completeness to keep a bin (default 50).")
    b.add_row("--bin-max-contamination","Maximum contamination to keep a bin (default 10).")
    console.print(b)
    console.print(Rule(style="dim"))


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def build_parser() -> argparse.ArgumentParser:
    p = CustomArgumentParser(
        prog="VMP_end2end.py",
        description="QC → Assembly → VPAC → Clustering → Binning (config-driven).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required
    p.add_argument("-i", "--input",  required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-cf", "--config", required=True)

    # Assembly
    p.add_argument("--assembly", choices=["spades", "megahit"], default="spades")
    p.add_argument("--spades-metaviral", action="store_true")

    # VPAC
    p.add_argument("--vpac", choices=["single", "dual", "both"], default="both")
    p.add_argument("--device", default="cuda:0")
    p.add_argument("--lmin", type=int, default=3000)
    p.add_argument("--lmax", type=int, default=100000)

    # Clustering
    p.add_argument("--cluster", choices=["cd-hit_mmseq", "skani_pyleiden"], default="skani_pyleiden")

    # Common
    p.add_argument("-t", "--threads", type=int, default=60)
    p.add_argument("--parallel", type=int, default=2)

    # Resume
    g = p.add_mutually_exclusive_group()
    g.add_argument("--resume", dest="resume", action="store_true")
    g.add_argument("--no-resume", dest="resume", action="store_false")
    p.set_defaults(resume=True)
    p.add_argument("--from-step", choices=["qc", "assembly", "vpac_single", "vpac_dual", "clustering", "binning"])

    # Global log control
    lg = p.add_mutually_exclusive_group()
    lg.add_argument("--quiet", dest="quiet", action="store_true", help="Hide sub-step logs.")
    lg.add_argument("--show-logs", dest="quiet", action="store_false", help="Show sub-step logs.")
    p.set_defaults(quiet=True)  # default: hide logs

    # Binning control (only requested flags)
    p.add_argument("--no-binning", dest="binning", action="store_false", help="Disable Binning stage.")
    p.set_defaults(binning=True)
    p.add_argument("--bin-minfasta-kb", type=int, default=200)
    p.add_argument("--bin-min-len", type=int, default=500)
    p.add_argument("--bin-min-completeness", type=float, default=50.0)
    p.add_argument("--bin-max-contamination", type=float, default=10.0)

    return p


# --------------------------- Utilities & Envs ---------------------------------
def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Run a command respecting global QUIET. Hide stdout/stderr when QUIET=True."""
    global QUIET
    pretty = " ".join(shlex.quote(x) for x in cmd)
    if not QUIET:
        console.print(f"[bold white]$[/bold white] {pretty}")
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=(subprocess.DEVNULL if QUIET else None),
            stderr=(subprocess.DEVNULL if QUIET else None),
        )
    except FileNotFoundError as e:
        console.print(f"[red]Command not found:[/red] {pretty}")
        raise
    if proc.returncode != 0:
        if QUIET:
            console.print(f"[red]Step failed (exit {proc.returncode}). Re-run with --show-logs for details.[/red]")
        raise SystemExit(proc.returncode)


def ensure_exists(p: Path, hint: str = "") -> None:
    if not p.exists():
        msg = f"Required path not found: {p}"
        if hint:
            msg += f"  ({hint})"
        raise FileNotFoundError(msg)


def is_nonempty_dir(p: Path) -> bool:
    return p.exists() and p.is_dir() and any(p.iterdir())


def any_fastas(p: Path) -> bool:
    if not p.exists() or not p.is_dir():
        return False
    for ext in ("*.fa", "*.fna", "*.fasta", "*.contigs.fa", "*.contigs.fasta"):
        if any(p.glob(ext)):
            return True
    for sub in p.iterdir():
        if sub.is_dir():
            for ext in ("*.fa", "*.fna", "*.fasta", "*.contigs.fa", "*.contigs.fasta"):
                if any(sub.glob(ext)):
                    return True
    return False


def any_bins(clean_root: Path) -> bool:
    if not clean_root.exists():
        return False
    return any(clean_root.rglob("*.fa"))


def _first_existing_key(d: dict, keys: list[str]) -> Optional[str]:
    for k in keys:
        if k in d and d[k]:
            return str(d[k])
    return None


def load_envs_from_config(cfg_path: Path) -> dict:
    envs = {"qc": None, "assembly": None, "vpac_single": None, "vpac_dual": None, "clustering": None, "binning": None}
    if yaml is None:
        console.print("[yellow]PyYAML not available; running all steps in the current Python env.[/yellow]")
        return envs
    try:
        cfg = yaml.safe_load(Path(cfg_path).read_text(encoding="utf-8"))
        if not isinstance(cfg, dict):
            return envs
        envs["qc"]          = _first_existing_key(cfg, ["VMP_env", "VMP_new", "QC_env"])
        envs["assembly"]    = _first_existing_key(cfg, ["VMP_env", "VMP_new", "Assembly_env"])
        envs["vpac_single"] = _first_existing_key(cfg, ["VPAC_single", "VPAC_single_env"])
        envs["vpac_dual"]   = _first_existing_key(cfg, ["VPAC_dual", "VMP_dual", "VPAC_dual_env"])
        envs["clustering"]  = _first_existing_key(cfg, ["VMP_env", "VMP_new", "Clustering_env"])
        envs["binning"]     = _first_existing_key(cfg, ["Binning", "Binning_env"])
    except Exception as e:
        console.print(f"[yellow]Failed to parse config.yml for envs: {e}. Using current env.[/yellow]")
    return envs


def load_checkm2_diamond_from_config(cfg_path: Path) -> Optional[str]:
    """
    Expect a key like 'checkM_db' or 'CheckM_db' pointing to a base folder.
    We will append 'uniref100.KO.1.dmnd'.
    """
    if yaml is None:
        return None
    try:
        cfg = yaml.safe_load(Path(cfg_path).read_text(encoding="utf-8"))
        if not isinstance(cfg, dict):
            return None
        base = _first_existing_key(cfg, ["checkM_db", "CheckM_db", "checkm_db"])
        if not base:
            return None
        dmnd = Path(str(base)).expanduser().resolve() / "uniref100.KO.1.dmnd"
        return str(dmnd)
    except Exception:
        return None


def run_in_env(env_prefix: Optional[str], script: Path, args_list: list[str]) -> None:
    if env_prefix:
        python_bin = Path(env_prefix) / "bin" / "python"
        if python_bin.exists():
            run([str(python_bin), str(script)] + args_list)
            return
        run(["conda", "run", "--no-capture-output", "--prefix", str(env_prefix), "python", str(script)] + args_list)
        return
    run([sys.executable, str(script)] + args_list)


def step_header(title: str, subtitle: str = "") -> None:
    text = f"[bold]{title}[/bold]"
    if subtitle:
        text += f"\n[dim]{subtitle}[/dim]"
    console.print(Panel.fit(text, border_style="cyan", padding=(1, 2)))


# ------------------ Inputs for Clustering and Binning -------------------------
def get_cluster_input_dir(args, vpac_out_dir: Path) -> Path:
    """
    Clustering keeps using VPAC viral outputs as before.
    """
    mvc = vpac_out_dir / "mVC"
    mvc_dual = vpac_out_dir / "mVC_dual"
    if args.vpac == "single":
        return mvc
    elif args.vpac == "dual":
        return mvc_dual
    else:
        return mvc_dual if is_nonempty_dir(mvc_dual) else mvc


def get_binning_nonviral_dir(args, vpac_out_dir: Path) -> Path:
    """
    Binning must target NON-viral sets from VPAC:
      - single: mNVC/
      - dual:   mNVC_dual/  (also accept mMVC_dual/ if that's your naming)
      - both:   prefer mNVC_dual/ (or mMVC_dual/), else mNVC/
    """
    mnvc_single = vpac_out_dir / "mNVC"
    mnvc_dual   = vpac_out_dir / "mNVC_dual"
    mmvc_dual   = vpac_out_dir / "mMVC_dual"  # compatibility if your folder is named like this

    if args.vpac == "single":
        return mnvc_single
    elif args.vpac == "dual":
        # prefer mNVC_dual; fall back to mMVC_dual if present
        return mnvc_dual if is_nonempty_dir(mnvc_dual) else mmvc_dual
    else:
        # both: prefer dual non-viral if present
        if is_nonempty_dir(mnvc_dual):
            return mnvc_dual
        if is_nonempty_dir(mmvc_dual):
            return mmvc_dual
        return mnvc_single


def print_plan(args):
    tbl = Table(show_header=True, header_style="bold magenta")
    tbl.add_column("Stage", style="cyan", no_wrap=True)
    tbl.add_column("Key Options / Paths", style="white")

    tbl.add_row("QC",
                f"raw FASTQ: {args.input}\n"
                f"out: {args.output}\nconfig: {args.config}")

    asm_line = (
        f"tool: {args.assembly} | threads: {args.threads}\n"
        f"→ raw_contigs: {Path(args.output) / 'raw_contigs'}"
    )
    if args.assembly == "spades" and args.spades_metaviral:
        asm_line += " | --metaviral"
    tbl.add_row("Assembly", asm_line)

    vpac_out_dir = Path(args.output) / "Viral_contigs" / f"{args.assembly}-single"
    tbl.add_row("VPAC",
                f"mode: {args.vpac} | device: {args.device}\n"
                f"in: {Path(args.output) / 'raw_contigs'}\n"
                f"out: {vpac_out_dir}")

    cluster_in_dir = get_cluster_input_dir(args, vpac_out_dir)
    tbl.add_row("Clustering",
                f"tool: {args.cluster} | parallel: {args.parallel}\n"
                f"in: {cluster_in_dir}\n"
                f"out: {Path(args.output) / ('cluster_' + args.cluster)}")

    bin_out_dir = vpac_out_dir / "Bin"
    tbl.add_row("Binning",
                f"enabled: {args.binning} | threads: {args.threads} | parallel: {args.parallel}\n"
                f"non-viral in: decided at runtime (mNVC / mNVC_dual)\n"
                f"out: {bin_out_dir}\n"
                f"filters: minfasta-kb={args.bin_minfasta_kb}, min-len={args.bin_min_len}, "
                f"comp≥{args.bin_min_completeness}, cont≤{args.bin_max_contamination}\n"
                f"CheckM2 dmnd: from config 'checkM_db' + 'uniref100.KO.1.dmnd'")
    console.print(tbl)
    console.print(Rule())


def detect_resume_step(args) -> str:
    out_root        = Path(args.output).resolve()
    clean_reads_dir = out_root / "Clean_reads"
    raw_contigs_dir = out_root / "raw_contigs"
    vpac_out_dir    = out_root / "Viral_contigs" / f"{args.assembly}-single"
    mvc_dir         = vpac_out_dir / "mVC"
    mvc_dual_dir    = vpac_out_dir / "mVC_dual"
    cluster_out_dir = out_root / f"cluster_{args.cluster}"
    bin_out_dir     = vpac_out_dir / "Bin"

    if not is_nonempty_dir(clean_reads_dir):
        return "qc"
    if not any_fastas(raw_contigs_dir):
        return "assembly"
    if args.vpac == "single":
        if not is_nonempty_dir(mvc_dir):
            return "vpac_single"
    elif args.vpac == "dual":
        if not is_nonempty_dir(mvc_dual_dir):
            return "vpac_dual"
    else:
        if not is_nonempty_dir(mvc_dir):
            return "vpac_single"
        if not is_nonempty_dir(mvc_dual_dir):
            return "vpac_dual"
    if not is_nonempty_dir(cluster_out_dir):
        return "clustering"
    if args.binning:
        clean_root = bin_out_dir / "clean_bin"
        if not any_bins(clean_root):
            return "binning"
    return "done"


def run_qc(args, raw_fastq_dir: Path, out_root: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    qc_env = envs.get('qc')
    step_header("Step 1 — QC", "fastp / host removal / summaries")
    run_in_env(
        qc_env,
        QC_SCRIPT,
        ["-r", str(raw_fastq_dir), "-out", str(out_root), "-cf", str(cfg_path)] +
        ([] if args.threads is None else ["-p", str(args.threads)])
    )


def run_assembly(args, clean_reads_dir: Path, out_root: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    asm_env = envs.get('assembly')
    step_header("Step 2 — Assembly", f"{args.assembly.upper()} → raw_contigs/")
    asm_args = [
        "-r", str(clean_reads_dir),
        "-out", str(out_root),
        "-cf", str(cfg_path),
        "--tool", args.assembly,
        "--parallel", str(args.parallel),
        "-p", str(args.threads),
    ]
    if args.assembly == "spades" and args.spades_metaviral:
        asm_args.append("--metaviral")
    run_in_env(asm_env, ASM_SCRIPT, asm_args)


def run_vpac_single(args, raw_contigs_dir: Path, vpac_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    vs_env = envs.get('vpac_single')
    step_header("Step 3 — VPAC-single", "FNN classifier")
    run_in_env(
        vs_env,
        VPAC_SINGLE_SCRIPT,
        [
            "-i", str(raw_contigs_dir),
            "-o", str(vpac_out_dir),
            "-lmin", str(args.lmin),
            "-lmax", str(args.lmax),
            "-t", str(args.threads),
            "--parallel", "1",
            "-m", "FNN",
            "-cf", str(cfg_path),
        ],
    )


def run_vpac_dual(args, raw_contigs_dir: Path, vpac_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    vd_env = envs.get('vpac_dual')
    step_header("Step 3 — VPAC-dual", "dual-path classifier")
    run_in_env(
        vd_env,
        VPAC_DUAL_SCRIPT,
        ["-i", str(raw_contigs_dir), "-o", str(vpac_out_dir), "-d", args.device, "-cf", str(cfg_path)],
    )


def run_clustering(args, vpac_out_dir: Path, cluster_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    clu_env = envs.get('clustering')
    cluster_in_dir = get_cluster_input_dir(args, vpac_out_dir)
    if not is_nonempty_dir(cluster_in_dir):
        raise FileNotFoundError(
            f"Clustering input is empty or missing: {cluster_in_dir}\n"
            f"Expected VPAC viral outputs under 'mVC' or 'mVC_dual'."
        )
    step_header("Step 4 — Clustering", f"tool = {args.cluster}")
    run_in_env(
        clu_env,
        CLUST_SCRIPT,
        ["-i", str(cluster_in_dir), "-o", str(cluster_out_dir), "--tool", args.cluster,
         "--parallel", str(args.parallel), "-cf", str(cfg_path)],
    )


def run_binning(args, out_root: Path, vpac_out_dir: Path, cfg_path: Path):
    if not args.binning:
        console.print("[yellow]Binning stage disabled (--no-binning).[/yellow]")
        return
    envs = load_envs_from_config(cfg_path)
    bin_env = envs.get('binning')

    clean_reads_dir = out_root / "Clean_reads"
    bin_out_dir     = vpac_out_dir / "Bin"
    bin_out_dir.mkdir(parents=True, exist_ok=True)

    nonviral_dir = get_binning_nonviral_dir(args, vpac_out_dir)

    if not is_nonempty_dir(nonviral_dir):
        raise FileNotFoundError(
            f"Non-viral directory for binning is missing or empty: {nonviral_dir}\n"
            f"Expected VPAC non-viral outputs under 'mNVC' or 'mNVC_dual' (also accept 'mMVC_dual')."
        )

    dmnd = load_checkm2_diamond_from_config(cfg_path)
    if not dmnd:
        raise FileNotFoundError("Failed to resolve CheckM2 diamond db from config (checkM_db + 'uniref100.KO.1.dmnd').")
    if not Path(dmnd).exists():
        raise FileNotFoundError(f"CheckM2 diamond db not found: {dmnd}")

    step_header("Step 5 — Binning (non-viral)", "bowtie2/samtools → SemiBin2 → CheckM2 → dRep")

    bin_args = [
        "-cf", str(cfg_path),
        "-r",  str(clean_reads_dir),
        "-i",  str(nonviral_dir),
        "-out", str(bin_out_dir),
        "-t",  str(args.threads),
        "-p",  str(args.parallel),
        "--minfasta-kb", str(args.bin_minfasta_kb),
        "--min-len", str(args.bin_min_len),
        "--min-completeness", str(args.bin_min_completeness),
        "--max-contamination", str(args.bin_max_contamination),
        "--checkm2-db", dmnd,  # auto-resolved; always provided
    ]

    run_in_env(bin_env, BINNING_SCRIPT, bin_args)


def orchestrate(args):
    for script in (QC_SCRIPT, ASM_SCRIPT, VPAC_SINGLE_SCRIPT, VPAC_DUAL_SCRIPT, CLUST_SCRIPT, BINNING_SCRIPT):
        ensure_exists(script, "Component script missing in the same folder.")

    raw_fastq_dir   = Path(args.input).resolve()
    out_root        = Path(args.output).resolve()
    cfg_path        = Path(args.config).resolve()

    ensure_exists(raw_fastq_dir, "Input FASTQ dir")
    ensure_exists(cfg_path, "config.yml")
    out_root.mkdir(parents=True, exist_ok=True)

    clean_reads_dir = out_root / "Clean_reads"
    raw_contigs_dir = out_root / "raw_contigs"
    vpac_out_dir    = out_root / "Viral_contigs" / f"{args.assembly}-single"
    cluster_out_dir = out_root / f"cluster_{args.cluster}"
    bin_out_dir     = vpac_out_dir / "Bin"

    # Decide start step
    if args.from_step:
        start_step = args.from_step
    elif args.resume:
        start_step = detect_resume_step(args)
    else:
        start_step = "qc"

    # Show status
    status_tbl = Table(show_header=True, header_style="bold blue")
    status_tbl.add_column("Step", style="cyan", no_wrap=True)
    status_tbl.add_column("Path / Artifact", style="white")
    status_tbl.add_column("Status", style="white")

    def status_row(name, path, done):
        status_tbl.add_row(name, str(path), "[green]DONE[/green]" if done else "[yellow]PENDING[/yellow]")

    status_row("QC", clean_reads_dir, is_nonempty_dir(clean_reads_dir))
    status_row("Assembly", raw_contigs_dir, any_fastas(raw_contigs_dir))
    status_row("VPAC-single", vpac_out_dir / "mVC", is_nonempty_dir(vpac_out_dir / "mVC"))
    status_row("VPAC-dual", vpac_out_dir / "mVC_dual", is_nonempty_dir(vpac_out_dir / "mVC_dual"))
    status_row("Clustering", cluster_out_dir, is_nonempty_dir(cluster_out_dir))
    if args.binning:
        status_row("Binning", bin_out_dir / "clean_bin", any_bins(bin_out_dir / "clean_bin"))

    console.print(Panel(status_tbl, border_style="blue", title="Resume status", title_align="left"))
    console.print(Panel(f"Starting from: [bold]{start_step}[/bold]", border_style="magenta"))

    if start_step == "done":
        console.print(Panel.fit("[bold green]Everything already complete. Nothing to do.[/bold green]", border_style="green"))
        return

    # QC
    if start_step == "qc":
        run_qc(args, raw_fastq_dir, out_root, cfg_path)
    ensure_exists(clean_reads_dir, "QC should create Clean_reads/ under -o")

    # Assembly
    if start_step in ("qc", "assembly"):
        run_assembly(args, clean_reads_dir, out_root, cfg_path)
    ensure_exists(raw_contigs_dir, "Assembly should create raw_contigs/ under -o")

    # VPAC
    if args.vpac == "single":
        if start_step in ("qc", "assembly", "vpac_single"):
            run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
    elif args.vpac == "dual":
        if start_step in ("qc", "assembly", "vpac_dual"):
            run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
    else:  # both
        mvc = vpac_out_dir / "mVC"
        mvc_dual = vpac_out_dir / "mVC_dual"
        if start_step in ("qc", "assembly"):
            if not is_nonempty_dir(mvc):
                run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc_dual):
                run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
        elif start_step == "vpac_single":
            run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc_dual):
                run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
        elif start_step == "vpac_dual":
            run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc):
                run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)

    # Clustering
    if start_step in ("qc", "assembly", "vpac_single", "vpac_dual", "clustering"):
        run_clustering(args, vpac_out_dir, cluster_out_dir, cfg_path)

    # Binning (non-viral)
    if args.binning and start_step in ("qc", "assembly", "vpac_single", "vpac_dual", "clustering", "binning"):
        run_binning(args, out_root, vpac_out_dir, cfg_path)

    console.print(Panel.fit(
        f"[bold green]All done![/bold green]\n\n"
        f"[white]QC Clean reads:[/white] {clean_reads_dir}\n"
        f"[white]Assembly contigs:[/white] {raw_contigs_dir}\n"
        f"[white]VPAC output:[/white] {vpac_out_dir}\n"
        f"[white]Clustering output:[/white] {cluster_out_dir}\n"
        f"[white]Binning output:[/white] {bin_out_dir if args.binning else '(skipped)'}",
        border_style="green"
    ))


def main():
    global QUIET
    args = build_parser().parse_args()
    QUIET = bool(args.quiet)

    console.print(Panel.fit(
        Markdown("### VMP End-to-End\n**Pipeline:** QC → Assembly → VPAC → Clustering → Binning\nThis run is fully driven by your `config.yml` (envs, tools, models)."),
        border_style="blue"
    ))
    print_plan(args)
    orchestrate(args)


if __name__ == "__main__":
    main()
