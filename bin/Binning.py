import os
import sys
import glob
import shutil
import argparse
import subprocess
import yaml
import pandas as pd
import multiprocessing as mp

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.markdown import Markdown

console = Console()


def run_cmd(cmd, desc=None):
    if isinstance(cmd, str):
        shell = True
        cmd_str = cmd
    else:
        shell = False
        cmd_str = " ".join(cmd)

    if desc:
        console.print(f"[bold cyan]> {desc}[/bold cyan]")
    console.print(f"[dim]{cmd_str}[/dim]")
    proc = subprocess.run(cmd, shell=shell)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit {proc.returncode}): {cmd_str}")


def get_sample_list(raw_dir):
    samples = []
    for name in os.listdir(raw_dir):
        p = os.path.join(raw_dir, name)
        if os.path.isdir(p):
            samples.append(name)
    samples = sorted(samples)
    console.print(f"[bold green]Sample list:[/bold green] {samples}")
    return samples


def find_contig_fasta(contig_dir, sample):
    patterns = [
        os.path.join(contig_dir, f"{sample}*.fa"),
        os.path.join(contig_dir, f"{sample}*.fasta"),
        os.path.join(contig_dir, f"{sample}*.fna"),
        os.path.join(contig_dir, f"{sample}*.fa.gz"),
        os.path.join(contig_dir, f"{sample}*.fasta.gz"),
        os.path.join(contig_dir, f"{sample}*.fna.gz"),
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))

    candidates = sorted(set(candidates))
    if not candidates:
        raise FileNotFoundError(f"No contig fasta found for sample {sample} under {contig_dir}")
    if len(candidates) > 1:
        console.print(
            f"[yellow]Warning: multiple contig files for {sample}, use first one:[/yellow]\n"
            + "\n".join(candidates)
        )
    return candidates[0]


def filter_checkm2_bins(
    quality_tsv,
    bins_dir,
    clean_dir,
    min_completeness,
    max_contamination,
    ext="fa",
):
    os.makedirs(clean_dir, exist_ok=True)
    if not os.path.exists(quality_tsv):
        console.print(f"[red]CheckM2 report not found: {quality_tsv}[/red]")
        return 0

    df = pd.read_csv(quality_tsv, sep="\t")
    required = {"Name", "Completeness", "Contamination"}
    if not required.issubset(df.columns):
        raise ValueError(f"quality_report.tsv missing columns: {required - set(df.columns)}")

    filtered = df[
        (df["Completeness"] >= min_completeness)
        & (df["Contamination"] <= max_contamination)
    ].copy()

    console.print(
        f"  Total bins: {len(df)}, "
        f"after filter (Comp>={min_completeness}, Cont<={max_contamination}): "
        f"[bold]{len(filtered)}[/bold]"
    )

    passed = []
    for name in filtered["Name"]:
        src = os.path.join(bins_dir, f"{name}.{ext}")
        if not os.path.exists(src):
            console.print(f"[yellow]  Skip: bin fasta not found: {src}[/yellow]")
            continue
        dst = os.path.join(clean_dir, f"{name}.{ext}")
        shutil.copy2(src, dst)
        passed.append(name)

    filtered = filtered[filtered["Name"].isin(passed)]
    if len(filtered) > 0:
        out_qc = os.path.join(clean_dir, "quality_report.filtered.tsv")
        filtered.to_csv(out_qc, sep="\t", index=False)
        console.print(f"  Filtered quality report saved: {out_qc}")
    else:
        console.print("  No bins copied after filtering.")
    return len(filtered)


def process_sample(
    sample,
    raw_dir,
    contig_dir,
    mapping_root,
    semibin_root,
    checkm_root,
    clean_root,
    threads,
    minfasta_kb,
    min_len,
    min_completeness,
    max_contamination,
    checkm2_db,
    keep_temp,
):

    try:
        console.print(Rule(title=f"[bold]Sample {sample}[/bold]"))

        reads1 = os.path.join(raw_dir, sample, f"{sample}_host_removed_1.fastq")
        reads2 = os.path.join(raw_dir, sample, f"{sample}_host_removed_2.fastq")

        if not (os.path.exists(reads1) and os.path.exists(reads2)):
            console.print(
                f"[yellow]Skip {sample}: host_removed FASTQ not found under {os.path.join(raw_dir, sample)}[/yellow]"
            )
            return

        try:
            contig_fa = find_contig_fasta(contig_dir, sample)
        except FileNotFoundError as e:
            console.print(f"[yellow]{e}[/yellow]")
            return

        console.print(f"  Contigs: {contig_fa}")
        console.print(f"  Reads:   {reads1}, {reads2}")

        sample_map_dir = os.path.join(mapping_root, sample)
        os.makedirs(sample_map_dir, exist_ok=True)

        index_prefix = os.path.join(sample_map_dir, f"{sample}.contigs")
        sam = os.path.join(sample_map_dir, f"{sample}.sam")
        bam = os.path.join(sample_map_dir, f"{sample}.bam")
        mapped_bam = os.path.join(sample_map_dir, f"{sample}.mapped.bam")
        sorted_bam = os.path.join(sample_map_dir, f"{sample}.mapped.sorted.bam")

        run_cmd(
            f"bowtie2-build --threads {threads} -f {contig_fa} {index_prefix}",
            desc=f"{sample}: bowtie2-build"
        )

        run_cmd(
            f"bowtie2 -q --fr -x {index_prefix} -1 {reads1} -2 {reads2} "
            f"-S {sam} -p {threads}",
            desc=f"{sample}: bowtie2 mapping"
        )

        run_cmd(
            f"samtools view -@ {threads} -h -b -S {sam} -o {bam}",
            desc=f"{sample}: samtools view"
        )
        run_cmd(
            f"samtools view -@ {threads} -b -F 4 {bam} -o {mapped_bam}",
            desc=f"{sample}: samtools filter mapped"
        )
        run_cmd(
            f"samtools sort -@ {threads} {mapped_bam} -o {sorted_bam}",
            desc=f"{sample}: samtools sort"
        )
        run_cmd(
            f"samtools index {sorted_bam}",
            desc=f"{sample}: samtools index"
        )

        if not keep_temp:
            for f in [sam, bam, mapped_bam]:
                if os.path.exists(f):
                    os.remove(f)

        sample_semibin_dir = os.path.join(semibin_root, sample)
        os.makedirs(sample_semibin_dir, exist_ok=True)

        semibin_cmd = (
            f"SemiBin2 single_easy_bin "
            f"-i {contig_fa} "
            f"-b {sorted_bam} "
            f"-o {sample_semibin_dir} "
            f"--environment global "
            f"--minfasta-kb {minfasta_kb} "
            f"--random-seed 66 "
            f"--max-edges 200 "
            f"--min-len {min_len} "
            f"--sequencing-type short_read "
            f"--verbose "
            # f"--no-recluster "
            f"-t 1 "
            f"--compression none"
        )
        run_cmd(semibin_cmd, desc=f"{sample}: SemiBin2 single_easy_bin")

        bins_dir = os.path.join(sample_semibin_dir, "output_bins")
        if not os.path.isdir(bins_dir):
            console.print(f"[red]{sample}: SemiBin2 output_bins not found: {bins_dir}[/red]")
            return

        sample_checkm_dir = os.path.join(checkm_root, sample)
        os.makedirs(sample_checkm_dir, exist_ok=True)

        checkm_cmd = (
            f"checkm2 predict "
            f"-x fa "
            f"--database_path {checkm2_db} "
            f"--threads {threads} "
            f"--input {bins_dir} "
            f"--output-directory {sample_checkm_dir}"
        )
        run_cmd(checkm_cmd, desc=f"{sample}: CheckM2 predict")

        quality_tsv = os.path.join(sample_checkm_dir, "quality_report.tsv")
        sample_clean_dir = os.path.join(clean_root, sample)

        console.print(
            f"  Filtering bins by completeness>={min_completeness}, "
            f"contamination<={max_contamination}"
        )
        n_kept = filter_checkm2_bins(
            quality_tsv=quality_tsv,
            bins_dir=bins_dir,
            clean_dir=sample_clean_dir,
            min_completeness=min_completeness,
            max_contamination=max_contamination,
            ext="fa",
        )
        console.print(f"[bold green]Sample {sample}: kept {n_kept} bins after filtering.[/bold green]")

    except Exception as e:
        console.print(f"[red]Error in sample {sample}: {e}[/red]")


def run_drep(clean_root, drep_out_dir, threads):
    console.print(Rule(title="[bold]dRep dereplication across all samples[/bold]"))

    if not os.path.isdir(clean_root):
        console.print(f"[yellow]clean_bin directory not found: {clean_root}, skip dRep.[/yellow]")
        return

    genomes = []
    for sample in sorted(os.listdir(clean_root)):
        sample_dir = os.path.join(clean_root, sample)
        if not os.path.isdir(sample_dir):
            continue
        genomes.extend(sorted(glob.glob(os.path.join(sample_dir, "*.fa"))))

    if not genomes:
        console.print("[yellow]No .fa bins found under clean_bin; skip dRep.[/yellow]")
        return

    os.makedirs(drep_out_dir, exist_ok=True)
    genomes_dir = os.path.join(drep_out_dir, "genomes")
    os.makedirs(genomes_dir, exist_ok=True)

    copied = []
    for g in genomes:
        dst = os.path.join(genomes_dir, os.path.basename(g))
        if not os.path.exists(dst):
            shutil.copy2(g, dst)
        copied.append(dst)

    drep_cmd = ["dRep", "dereplicate", drep_out_dir, "-p", str(threads), "-g"]
    drep_cmd.extend(copied)

    run_cmd(drep_cmd, desc="dRep dereplicate")

    console.print(
        f"[bold green]dRep dereplication finished. Non-redundant bins and reports are in: {drep_out_dir}/dereplicated_genomes[/bold green]"
    )


def custom_help():
    console = Console()

    intro = (
        "The 'Binning' script performs contig binning for each sample:\n"
        "  1) Map host-removed reads to contigs and build sorted BAM;\n"
        "  2) Run SemiBin2 to generate bins;\n"
        "  3) Run CheckM2 for quality assessment and filter bins by completeness/contamination;\n"
        "  4) Optionally run dRep across all clean bins to obtain a non-redundant MAG set."
    )
    console.print(Panel(intro, border_style="cyan", title="Binning", title_align="left"))

    examples = Markdown(
        "\n**Examples:**\n"
        "```bash\n"
        "Binning \\\n"
        "  -cf .../config.yml \\\n"
        "  -r  .../QC/Clean_reads \\\n"
        "  -i  .../Viral_contigs/megahit-single/mNVC_dual \\\n"
        "  -out .../Viral_contigs/megahit-single/Bin \\\n"
        "  -t 20 -p 4 \\\n"
        "  --minfasta-kb 25 --min-len 500 \\\n"
        "  --min-completeness 50 --max-contamination 10\n"
        "```"
    )
    console.print(examples)

    reads_md = Markdown(
        "**Input reads location (-r):**\n"
        "Directory of clean reads produced by QC, containing host-removed paired-end FASTQ files:\n"
        "```text\n"
        "/path/to/Clean_reads/A1/A1_host_removed_1.fastq\n"
        "/path/to/Clean_reads/A1/A1_host_removed_2.fastq\n"
        "/path/to/Clean_reads/A2/A2_host_removed_1.fastq\n"
        "/path/to/Clean_reads/A2/A2_host_removed_2.fastq\n"
        "...\n"
        "```\n"
    )
    contig_md = Markdown(
        "**Contig directory (-i):**\n"
        "Directory containing contigs for each sample, e.g.:\n"
        "```text\n"
        "/path/to/mNVC_dual/A1.mNVC.fasta\n"
        "/path/to/mNVC_dual/A2.mNVC.fasta\n"
        "```\n"
        "The script will search for files matching pattern `SAMPLE*.fa/fasta/fna`."
    )
    out_md = Markdown(
        "**Output directory (-out):**\n"
        "Will contain subdirectories:\n"
        "- `mapping/`  : Bowtie2 + samtools mapping results\n"
        "- `semibin2/`: SemiBin2 outputs (`output_bins/` per sample)\n"
        "- `checkm2/` : CheckM2 quality reports per sample\n"
        "- `clean_bin/`: filtered bins and `quality_report.filtered.tsv` per sample\n"
        "- `dRep/`    : non-redundant MAGs and dRep reports (if dRep is enabled)\n"
    )

    console.print(reads_md)
    console.print(contig_md)
    console.print(out_md)
    console.print(Rule())

    console.print("[bold]Global parameters[/bold]")
    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-cf, --config", "Path to config.yml for env/dbs.")
    g_tbl.add_row("-r, --raw_data", "Directory containing per-sample host_removed FASTQ subdirectories.")
    g_tbl.add_row("-i, --contig_dir", "Directory containing per-sample contig FASTA files.")
    g_tbl.add_row("-out, --out_dir", "Output directory for binning results.")
    g_tbl.add_row("-t, --thread", "Threads per sample.")
    g_tbl.add_row("-p, --parallel", "Number of samples to process in parallel (default: 1).")
    console.print(Panel(g_tbl, border_style="cyan", title="Global parameters", title_align="left"))

    console.print()
    console.print("[bold]Binning and filtering parameters[/bold]")
    b_tbl = Table(show_header=False, box=None, pad_edge=False)
    b_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    b_tbl.add_column("Description", style="white")
    b_tbl.add_row("--minfasta-kb", "Minimum bin size in Kbp (default: 200).")
    b_tbl.add_row("--min-len", "Minimal contig length for binning (default: 500).")
    b_tbl.add_row("--min-completeness", "Minimum completeness to keep a bin (default: 50).")
    b_tbl.add_row("--max-contamination", "Maximum contamination to keep a bin (default: 10).")
    b_tbl.add_row("--skip-drep", "Skip dRep dereplication across clean bins (default: run dRep).")
    b_tbl.add_row("--keep-temp", "Keep intermediate mapping files (sam/bam). Default: delete.")
    console.print(Panel(b_tbl, border_style="magenta", title="Binning and filtering", title_align="left"))

    console.print(Rule(style="dim"))


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def main():
    parser = CustomArgumentParser(description="Contig binning pipeline using SemiBin2, CheckM2, and optional dRep.")
    parser.add_argument("-cf", "--config", type=str, help="Path to config.yml for env/tools/dbs")
    parser.add_argument("-r", "--raw_data", type=str, required=True, help="Clean_reads directory from QC (contains per-sample host_removed fastq).")
    parser.add_argument("-i", "--contig_dir", type=str, required=True, help="Directory with per-sample contig FASTA files.")
    parser.add_argument("-out", "--out_dir", type=str, required=True, help="Output directory for binning results.")
    parser.add_argument("-t", "--thread", type=int, default=20, help="Threads per sample for bowtie2/samtools/SemiBin2/CheckM2 (default: 20).")
    parser.add_argument("-p", "--parallel", type=int, default=1, help="Number of samples to process in parallel (default: 1).")
    
    parser.add_argument("--minfasta-kb", type=int, default=200, help="Minimum bin size in Kbp (default: 200).")
    parser.add_argument("--min-len", type=int, default=500,  help="Minimal contig length used for binning (default: 500).")
    parser.add_argument("--min-completeness", type=float, default=50.0, help="Minimum completeness to keep a bin (default: 50).")
    parser.add_argument("--max-contamination", type=float, default=10.0,help="Maximum contamination to keep a bin (default: 10).")
    parser.add_argument("--checkm2-db", type=str, default=None, help="Optional: path to CheckM2 diamond db (.dmnd). If not set, use config.yml checkM_db/uniref100.KO.1.dmnd")
    parser.add_argument("--skip-drep", action="store_true", help="Skip dRep dereplication across clean bins (default: run dRep).")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate mapping files (sam/bam). Default: delete.")

    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config, "r") as f:
            raw = yaml.safe_load(f) or {}
            cfg.update(raw)

    bin_env = cfg.get("Binning")
    if bin_env:
        bin_env = str(bin_env).strip().strip("'\"")
        env_bin = os.path.join(bin_env, "bin")
        if os.path.isdir(env_bin):
            os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
            console.print(f"[green]Using Binning env bin:[/green] {env_bin}")
        else:
            console.print(f"[yellow]Warning: Binning env bin not found: {env_bin}[/yellow]")

    if args.checkm2_db:
        checkm2_db = args.checkm2_db
    else:
        checkm_root = cfg.get("checkM_db")
        if not checkm_root:
            console.print(
                "[red]Error: checkM_db not set in config.yml and --checkm2-db not provided.[/red]"
            )
            sys.exit(1)
        checkm_root = str(checkm_root).strip().strip("'\"")
        checkm2_db = os.path.join(checkm_root, "uniref100.KO.1.dmnd")

    console.print(f"[bold cyan]CheckM2 database:[/bold cyan] {checkm2_db}")

    raw_dir = args.raw_data
    contig_dir = args.contig_dir
    out_dir = args.out_dir
    threads = args.thread

    mapping_root = os.path.join(out_dir, "mapping")
    semibin_root = os.path.join(out_dir, "semibin2")
    checkm_root = os.path.join(out_dir, "checkm2")
    clean_root = os.path.join(out_dir, "clean_bin")
    drep_root = os.path.join(out_dir, "dRep")

    for d in [mapping_root, semibin_root, checkm_root, clean_root]:
        os.makedirs(d, exist_ok=True)

    samples = get_sample_list(raw_dir)

    jobs = []
    for sample in samples:
        jobs.append(
            (
                sample,
                raw_dir,
                contig_dir,
                mapping_root,
                semibin_root,
                checkm_root,
                clean_root,
                threads,
                args.minfasta_kb,
                args.min_len,
                args.min_completeness,
                args.max_contamination,
                checkm2_db,
                args.keep_temp,
            )
        )

    if args.parallel > 1:
        console.print(
            f"[bold cyan]Processing samples in parallel: {args.parallel} workers[/bold cyan]"
        )
        with mp.Pool(processes=args.parallel) as pool:
            pool.starmap(process_sample, jobs)
    else:
        for job in jobs:
            process_sample(*job)

    if args.skip_drep:
        console.print("[yellow]Skipping dRep dereplication (--skip-drep specified).[/yellow]")
    else:
        run_drep(clean_root, drep_root, threads)

    console.print(Rule())
    console.print("[bold cyan]Binning finished for all samples.[/bold cyan]")


if __name__ == "__main__":
    main()
