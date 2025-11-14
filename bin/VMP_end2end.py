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
        "VMP_end2end.py -i /raw -o /run -cf ./config.yml --assembly megahit --vpac single --cluster cd-hit_mmseq -t 48\n"
        "```\n"
    )

    io_md = Markdown(
        "**I/O conventions:**\n"
        "- [QC]      input: `-i` (raw paired FASTQ parent) → output: `-o/Clean_reads/`\n"
        "- [Assembly] input: `-o/Clean_reads/` → output: `-o/raw_contigs/`\n"
        "- [VPAC]     input: `-o/raw_contigs/` → output: `-o/Viral_contigs/{spades|megahit}-single/` (mVC, mVC_dual)\n"
        "- [Clustering] input: VPAC output dir → output: `-o/cluster_{tool}/`\n"
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
    p.add_argument("-i", "--input",  required=True, help="Directory containing raw paired FASTQ (one subdir per sample).")
    p.add_argument("-o", "--output", required=True, help="Root output directory for the whole run.")
    p.add_argument("-cf", "--config", required=True, help="Path to config.yml consumed by component scripts.")
    p.add_argument("--assembly", choices=["spades", "megahit"], default="spades", help="Assembler to use.")
    p.add_argument("--spades-metaviral", action="store_true", help="If using SPAdes, add --metaviral.")
    p.add_argument("--vpac", choices=["single", "dual", "both"], default="both",
                   help="Run only VPAC-single, only VPAC-dual, or both (same in/out).")
    p.add_argument("--device", default="cuda:0", help="Device for VPAC-dual (e.g., cuda:0 / cpu).")
    p.add_argument("--cluster", choices=["cd-hit_mmseq", "skani_pyleiden"], default="skani_pyleiden",
                   help="Clustering / de-redundancy toolchain.")
    p.add_argument("-t", "--threads", type=int, default=60, help="Threads for tools that accept -t/--threads.")
    p.add_argument("--parallel", type=int, default=2, help="Parallel samples (where supported).")
    p.add_argument("--lmin", type=int, default=3000, help="VPAC-single min contig length.")
    p.add_argument("--lmax", type=int, default=100000, help="VPAC-single max contig length.")
    return p

def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Run a command; raise on non-zero exit, while echoing the command nicely."""
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

    step_header("Step 1 — QC", "fastp / host removal / summaries (see QC script)")
    run([sys.executable, str(QC_SCRIPT),
         "-r", str(raw_fastq_dir),
         "-o", str(out_root),
         "-cf", str(cfg_path),
    ] + ([] if args.threads is None else ["-t", str(args.threads)]))

    ensure_exists(clean_reads_dir, "QC should create Clean_reads/ under -o")

    step_header("Step 2 — Assembly", f"{args.assembly.upper()} → raw_contigs/")
    asm_cmd = [sys.executable, str(ASM_SCRIPT),
               "-r", str(clean_reads_dir),
               "-o", str(out_root),
               "--tool", args.assembly,
               "--parallel", str(args.parallel),
               "-t", str(args.threads)]
    if args.assembly == "spades" and args.spades_metaviral:
        asm_cmd.append("--metaviral")
    run(asm_cmd)

    ensure_exists(raw_contigs_dir, "Assembly should create raw_contigs/ under -o")

    step_header("Step 3 — VPAC classification", f"mode = {args.vpac}")
    if args.vpac in ("single", "both"):
        run([sys.executable, str(VPAC_SINGLE_SCRIPT),
             "-i", str(raw_contigs_dir),
             "-o", str(vpac_out_dir),
             "-lmin", str(args.lmin),
             "-lmax", str(args.lmax),
             "-t", str(args.threads),
             "--parallel", "1",
             "-m", "FNN",
             "-cf", str(cfg_path),
        ])

    if args.vpac in ("dual", "both"):
        run([sys.executable, str(VPAC_DUAL_SCRIPT),
             "-i", str(raw_contigs_dir),
             "-o", str(vpac_out_dir),
             "-d", args.device,
             "-cf", str(cfg_path),
        ])

    step_header("Step 4 — Clustering", f"tool = {args.cluster}")
    run([sys.executable, str(CLUST_SCRIPT),
         "-i", str(vpac_out_dir),
         "-o", str(cluster_out_dir),
         "--tool", args.cluster,
         "--parallel", str(args.parallel),
         "-cf", str(cfg_path),
    ])

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
