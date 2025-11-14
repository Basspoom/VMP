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

QC_SCRIPT          = HERE / "QC.py"
ASM_SCRIPT         = HERE / "Assembly.py"
VPAC_SINGLE_SCRIPT = HERE / "VPAC_single.py"
VPAC_DUAL_SCRIPT   = HERE / "VPAC_dual.py"
CLUST_SCRIPT       = HERE / "Clustering.py"

console = Console()

def custom_help():
    intro = (
        "The 'VMP_end2end' orchestrator runs the complete pipeline:\n"
        "[bold]QC → Assembly → VPAC → Clustering[/bold]\n\n"
        "It wires the outputs/inputs across steps using your folder conventions and passes a single "
        "[bold]config.yml[/bold] to all components for environments, tools and models."
    )

    examples_md = Markdown(
        "\n**Examples:**\n"
        "```\n"
        "VMP_end2end.py \\\n"
        "  -i /data/raw_fastq \\\n"
        "  -o /data/VMP_run/out \\\n"
        "  -cf /data/VMP_run/config.yml \\\n"
        "  --assembly spades --spades-metaviral \\\n"
        "  --vpac both --device cuda:0 \\\n"
        "  --cluster skani_pyleiden \\\n"
        "  -t 60 --parallel 2\n"
        "\n"
        "VMP_end2end.py -i /raw -o /run -cf ./config.yml --assembly megahit --vpac single --cluster cd-hit_mmseq -t 48 --resume\n"
        "```\n"
    )

    io_md = Markdown(
        "**I/O conventions:**\n"
        "- [QC]        input: `-i` (raw paired FASTQ parent) → output: `-o/Clean_reads/`\n"
        "- [Assembly]  input: `-o/Clean_reads/` → output: `-o/raw_contigs/`\n"
        "- [VPAC]      input: `-o/raw_contigs/` → output: `-o/Viral_contigs/{spades|megahit}-single/` (mVC, mVC_dual)\n"
        "- [Clustering] input: VPAC output dir → output: `-o/cluster_{tool}/`\n"
        "\n"
        "[dim]*This orchestrator supports auto-resume and runs each step in its configured Conda env.*[/dim]"
    )

    console.print(Panel(intro, border_style="cyan", title="VMP End-to-End", title_align="left"))
    console.print(examples_md)
    console.print(io_md)
    console.print()

    out_tbl = Table(show_header=True, header_style="bold blue")
    out_tbl.add_column("Files/Dirs", style="cyan", no_wrap=False)
    out_tbl.add_column("Description", style="white")
    out_tbl.add_row("~/Clean_reads/", "QC-cleaned FASTQ pairs (per sample).")
    out_tbl.add_row("~/raw_contigs/", "Assembled contigs directory (shared for all samples).")
    out_tbl.add_row("~/Viral_contigs/<assembler>-single/mVC", "VPAC-single results (per-sample subfolders).")
    out_tbl.add_row("~/Viral_contigs/<assembler>-single/mVC_dual", "VPAC-dual results (if --vpac dual/both).")
    out_tbl.add_row("~/cluster_<tool>/", "Clustering results (cd-hit_mmseq or skani_pyleiden).")
    console.print(Panel(out_tbl, border_style="blue", title="Key Outputs", title_align="left"))
    console.print(Panel("Every component reads tool/env/model paths from the same config.yml.",
                        border_style="cyan"))
    console.print()

    console.print("[bold]Detailed parameters[/bold]\n")

    console.print(Panel("Usage: VMP_end2end.py [OPTIONS] -i RAW_FASTQ_DIR -o OUT_DIR -cf CONFIG",
                        border_style="cyan", title="Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-i, --input", "Directory of raw paired-end FASTQ (one subdir per sample).")
    g_tbl.add_row("-o, --output", "Root output directory for the whole run.")
    g_tbl.add_row("-cf, --config", "Path to config.yml (envs/tools/models for all steps).")
    g_tbl.add_row("-t, --threads", "Threads for tools that accept -t/--threads (default: 60).")
    g_tbl.add_row("--parallel", "Parallel samples where supported (default: 2).")
    g_tbl.add_row("--resume/--no-resume", "Auto-detect finished steps and continue (default: on).")
    g_tbl.add_row("--from-step", "Override resume and start from a specific step.")
    console.print(g_tbl)
    console.print()

    console.print(Panel("Assembly options", border_style="magenta"))
    a_tbl = Table(show_header=False, box=None, pad_edge=False)
    a_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    a_tbl.add_column("Description", style="white")
    a_tbl.add_row("--assembly", "Assembler: 'spades' or 'megahit' (default: spades).")
    a_tbl.add_row("--spades-metaviral", "If using SPAdes, enable metaviral mode.")
    console.print(a_tbl)
    console.print()

    console.print(Panel("VPAC options", border_style="green"))
    v_tbl = Table(show_header=False, box=None, pad_edge=False)
    v_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    v_tbl.add_column("Description", style="white")
    v_tbl.add_row("--vpac", "Run 'single', 'dual', or 'both' (default: both).")
    v_tbl.add_row("--device", "Device for VPAC-dual (e.g., cuda:0 / cpu).")
    v_tbl.add_row("--lmin/--lmax", "VPAC-single length bounds (defaults: 3000 / 100000).")
    console.print(v_tbl)
    console.print()

    console.print(Panel("Clustering options", border_style="yellow"))
    c_tbl = Table(show_header=False, box=None, pad_edge=False)
    c_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    c_tbl.add_column("Description", style="white")
    c_tbl.add_row("--cluster", "Toolchain: 'cd-hit_mmseq' or 'skani_pyleiden' (default: skani_pyleiden).")
    console.print(c_tbl)
    console.print(Rule(style="dim"))


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def build_parser() -> argparse.ArgumentParser:
    p = CustomArgumentParser(
        prog="VMP_end2end.py",
        description="End-to-end controller for VMP: QC → Assembly → VPAC → Clustering (config-driven).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required high-level
    p.add_argument("-i", "--input",  required=True, help="Directory containing raw paired FASTQ (one subdir per sample).")
    p.add_argument("-o", "--output", required=True, help="Root output directory for the whole run.")
    p.add_argument("-cf", "--config", required=True, help="Path to config.yml consumed by component scripts.")
    # Assembly choice
    p.add_argument("--assembly", choices=["spades", "megahit"], default="spades", help="Assembler to use.")
    p.add_argument("--spades-metaviral", action="store_true", help="If using SPAdes, add --metaviral.")
    # VPAC mode
    p.add_argument("--vpac", choices=["single", "dual", "both"], default="both",
                   help="Run only VPAC-single, only VPAC-dual, or both (same in/out).")
    p.add_argument("--device", default="cuda:0", help="Device for VPAC-dual (e.g., cuda:0 / cpu).")
    # Clustering
    p.add_argument("--cluster", choices=["cd-hit_mmseq", "skani_pyleiden"], default="skani_pyleiden",
                   help="Clustering / de-redundancy toolchain.")
    # Common resources
    p.add_argument("-t", "--threads", type=int, default=60, help="Threads for tools that accept -t/--threads.")
    p.add_argument("--parallel", type=int, default=2, help="Parallel samples (where supported).")
    # VPAC-single defaults
    p.add_argument("--lmin", type=int, default=3000, help="VPAC-single min contig length.")
    p.add_argument("--lmax", type=int, default=100000, help="VPAC-single max contig length.")
    # Resume & control
    g = p.add_mutually_exclusive_group()
    g.add_argument("--resume", dest="resume", action="store_true", help="Auto-detect finished steps and continue from the latest.")
    g.add_argument("--no-resume", dest="resume", action="store_false", help="Disable auto-resume; run all steps from the start.")
    p.set_defaults(resume=True)
    p.add_argument("--from-step", choices=["qc", "assembly", "vpac_single", "vpac_dual", "clustering"],
                   help="Start from a specific step (overrides --resume).")
    return p


# --------------------------- Utilities & Envs ---------------------------------
def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    pretty = " ".join(shlex.quote(x) for x in cmd)
    console.print(f"[bold white]$[/bold white] {pretty}")
    proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None)
    if proc.returncode != 0:
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

def _first_existing_key(d: dict, keys: list[str]) -> Optional[str]:
    for k in keys:
        if k in d and d[k]:
            return str(d[k])
    return None

def load_envs_from_config(cfg_path: Path) -> dict:
    envs = {"qc": None, "assembly": None, "vpac_single": None, "vpac_dual": None, "clustering": None}
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
    except Exception as e:
        console.print(f"[yellow]Failed to parse config.yml for envs: {e}. Falling back to current env.[/yellow]")
    return envs

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
    console.print(Panel.fit(text, border_style="cyan", padding=(1,2)))

def print_plan(args):
    tbl = Table(show_header=True, header_style="bold magenta")
    tbl.add_column("Stage", style="cyan", no_wrap=True)
    tbl.add_column("Key Options / Paths", style="white")

    tbl.add_row("QC",
        f"raw FASTQ in: {args.input}\n"
        f"out root: {args.output}\n"
        f"config: {args.config}"
    )
    asm_line = f"tool: {args.assembly} | threads: {args.threads}\n" \
               f"Assembly out base: {args.output}\n" \
               f"→ raw_contigs: {Path(args.output)/'raw_contigs'}"
    if args.assembly == "spades" and args.spades_metaviral:
        asm_line += " | --metaviral"
    tbl.add_row("Assembly", asm_line)

    tbl.add_row("VPAC",
        f"mode: {args.vpac} | device: {args.device}\n"
        f"input: {Path(args.output)/'raw_contigs'}\n"
        f"output: {Path(args.output)/'Viral_contigs'/f'{args.assembly}-single'}"
    )

    tbl.add_row("Clustering",
        f"tool: {args.cluster} | parallel: {args.parallel}\n"
        f"in: {Path(args.output)/'Viral_contigs'/f'{args.assembly}-single'}\n"
        f"out: {Path(args.output)/('cluster_'+args.cluster)}"
    )

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
    return "done"

def run_qc(args, raw_fastq_dir: Path, out_root: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    qc_env = envs.get('qc')
    step_header("Step 1 — QC", "fastp / host removal / summaries")
    run_in_env(qc_env, QC_SCRIPT,
               ["-r", str(raw_fastq_dir), "-out", str(out_root), "-cf", str(cfg_path)] +
               ([] if args.threads is None else ["-p", str(args.threads)]))

def run_assembly(args, clean_reads_dir: Path, out_root: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    asm_env = envs.get('assembly')
    step_header("Step 2 — Assembly", f"{args.assembly.upper()} → raw_contigs/")
    asm_args = ["-r", str(clean_reads_dir), "-out", str(out_root), "-cf", str(cfg_path),
                "--tool", args.assembly, "--parallel", str(args.parallel), "-p", str(args.threads)]
    if args.assembly == "spades" and args.spades_metaviral:
        asm_args.append("--metaviral")
    run_in_env(asm_env, ASM_SCRIPT, asm_args)

def run_vpac_single(args, raw_contigs_dir: Path, vpac_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    vs_env = envs.get('vpac_single')
    step_header("Step 3 — VPAC-single", "FNN classifier")
    run_in_env(vs_env, VPAC_SINGLE_SCRIPT,
               ["-i", str(raw_contigs_dir), "-o", str(vpac_out_dir),
                "-lmin", str(args.lmin), "-lmax", str(args.lmax),
                "-t", str(args.threads), "--parallel", "1", "-m", "FNN", "-cf", str(cfg_path)])

def run_vpac_dual(args, raw_contigs_dir: Path, vpac_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    vd_env = envs.get('vpac_dual')
    step_header("Step 3 — VPAC-dual", "dual-path classifier")
    run_in_env(vd_env, VPAC_DUAL_SCRIPT,
               ["-i", str(raw_contigs_dir), "-o", str(vpac_out_dir), "-d", args.device, "-cf", str(cfg_path)])

def run_clustering(args, vpac_out_dir: Path, cluster_out_dir: Path, cfg_path: Path):
    envs = load_envs_from_config(cfg_path)
    clu_env = envs.get('clustering')
    step_header("Step 4 — Clustering", f"tool = {args.cluster}")
    run_in_env(clu_env, CLUST_SCRIPT,
               ["-i", str(vpac_out_dir), "-o", str(cluster_out_dir), "--tool", args.cluster,
                "--parallel", str(args.parallel), "-cf", str(cfg_path)])

def orchestrate(args):
    for script in (QC_SCRIPT, ASM_SCRIPT, VPAC_SINGLE_SCRIPT, VPAC_DUAL_SCRIPT, CLUST_SCRIPT):
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

    start_step = None
    if args.from_step:
        start_step = args.from_step
    elif args.resume:
        start_step = detect_resume_step(args)
    else:
        start_step = "qc"

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
    console.print(Panel(status_tbl, border_style="blue", title="Resume status", title_align="left"))
    console.print(Panel(f"Starting from: [bold]{start_step}[/bold]", border_style="magenta"))

    if start_step == "done":
        console.print(Panel.fit("[bold green]Everything already complete. Nothing to do.[/bold green]", border_style="green"))
        return

    if start_step == "qc":
        run_qc(args, raw_fastq_dir, out_root, cfg_path)

    ensure_exists(clean_reads_dir, "QC should create Clean_reads/ under -o")

    if start_step in ("qc", "assembly"):
        run_assembly(args, clean_reads_dir, out_root, cfg_path)

    ensure_exists(raw_contigs_dir, "Assembly should create raw_contigs/ under -o")

    if args.vpac == "single":
        if start_step in ("qc", "assembly", "vpac_single"):
            run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
    elif args.vpac == "dual":
        if start_step in ("qc", "assembly", "vpac_dual"):
            run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
    else:
        mvc_dir = vpac_out_dir / "mVC"
        mvc_dual_dir = vpac_out_dir / "mVC_dual"
        if start_step in ("qc", "assembly"):
            if not is_nonempty_dir(mvc_dir):
                run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc_dual_dir):
                run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
        elif start_step == "vpac_single":
            run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc_dual_dir):
                run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
        elif start_step == "vpac_dual":
            run_vpac_dual(args, raw_contigs_dir, vpac_out_dir, cfg_path)
            if not is_nonempty_dir(mvc_dir):
                run_vpac_single(args, raw_contigs_dir, vpac_out_dir, cfg_path)

    if start_step in ("qc", "assembly", "vpac_single", "vpac_dual", "clustering"):
        run_clustering(args, vpac_out_dir, cluster_out_dir, cfg_path)

    console.print(Panel.fit(
        f"[bold green]All done![/bold green]\n\n"
        f"[white]QC Clean reads:[/white] {clean_reads_dir}\n"
        f"[white]Assembly contigs:[/white] {raw_contigs_dir}\n"
        f"[white]VPAC output:[/white] {vpac_out_dir}\n"
        f"[white]Clustering output:[/white] {cluster_out_dir}",
        border_style="green"))

def main():
    args = build_parser().parse_args()
    console.print(Panel.fit(
        Markdown("### VMP End-to-End\n**Pipeline:** QC → Assembly → VPAC → Clustering  \nThis run is fully driven by your `config.yml` (envs, tools, models)."),
        border_style="blue"))
    print_plan(args)
    orchestrate(args)

if __name__ == "__main__":
    main()
