import os
import re
import glob
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional, Dict
import sys
import shlex
from pathlib import Path
from typing import Union
import pandas as pd
from tqdm import tqdm
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.markdown import Markdown

try:
    import yaml
except Exception as e:
    raise RuntimeError("PyYAML is required. Please `conda/pip install pyyaml`. Original error: %s" % e)

SCORE_DIR = Path(__file__).resolve().parent / "score"

def run_score(script_name: str, *args: Union[Path, str]) -> None:
    script_path = SCORE_DIR / script_name
    if not script_path.is_file():
        raise FileNotFoundError(f"[score] script not found: {script_path}")
    cmd = " ".join(
        [shlex.quote(sys.executable), shlex.quote(str(script_path))] +
        [shlex.quote(str(a)) for a in args]
    )
    sh(cmd)


def sh(cmd: str, cwd: Optional[Path] = None) -> None:
    print("  >", cmd)
    proc = subprocess.run(cmd, shell=True, cwd=str(cwd) if cwd else None)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")


def conda_run(prefix: Path, inner_cmd: str) -> str:
    return f"conda run --no-capture-output --prefix {prefix} {inner_cmd}"

def conda_run_py(prefix: Path, script_name: str, *args: Union[Path, str]) -> None:
    script_path = SCORE_DIR / script_name
    if not script_path.is_file():
        raise FileNotFoundError(f"[score] script not found: {script_path}")
    cmd = " ".join(
        ["conda run --no-capture-output --prefix", shlex.quote(str(prefix)), "python", shlex.quote(str(script_path))] +
        [shlex.quote(str(a)) for a in args]
    )
    sh(cmd)


def Getfasta_list(input_dir: str) -> List[str]:
    path = Path(input_dir)
    fasta_list = []
    for name in os.listdir(path):
        p = path / name
        if p.is_dir():
            continue
        if name.endswith(".contigs.fa"):
            fasta_list.append(name[: name.index('.')])
    print('Your sample list:', fasta_list)
    return fasta_list


def run_seqkit(fasta_list: List[str], input_dir: str, out_dir: str, min_length: int, max_length: int, thread: int, env_single: Path) -> None:
    print(f"Filering contig length during {min_length} ~ {max_length} kb ...")
    Path(out_dir, "filtered_contigs").mkdir(parents=True, exist_ok=True)
    for item in fasta_list:
        in_fa = Path(input_dir) / f"{item}.contigs.fa"
        out_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        cmd = conda_run(env_single, f"seqkit seq -m {min_length} -M {max_length} -w 0 -j {thread} {in_fa} > {out_fa}")
        sh(cmd)
    print("Done filer!")


def run_deepvirfinder_single(fasta_list: List[str], out_dir: str, thread: int, parallel: int, env_single: Path, dvf_dir: Path) -> None:
    print("Begin DeepVirFinder...")
    Path(out_dir, "dvf", "dvf_all").mkdir(parents=True, exist_ok=True)
    dvf_py = dvf_dir / "dvf.py"

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        work = Path(out_dir) / "dvf" / item
        work.mkdir(parents=True, exist_ok=True)
        cmd1 = conda_run(env_single, f"python {dvf_py} -i {in_fa} -o {work} -c {thread}")
        sh(cmd1)
        sh(f"cp {work}/{item}.*.txt {Path(out_dir)/'dvf'/'dvf_all'}")

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done DeepVirFinder!")


def init_virsorter2_db(env_single: Path, vs2_db: Path) -> None:
    print("Initializing VirSorter2 DB source...")
    cmd = conda_run(env_single, f"virsorter config --init-source --db-dir {vs2_db}")
    sh(cmd)


def run_virsorter2_single(fasta_list: List[str], out_dir: str, thread: int, parallel: int, env_single: Path) -> None:
    print("Begin VirSorter2...")
    Path(out_dir, "vs2", "vs2_all").mkdir(parents=True, exist_ok=True)

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        wdir = Path(out_dir) / "vs2" / f"{item}.out"
        cmd1 = conda_run(env_single, f"virsorter run -w {wdir} -i {in_fa} -j {thread} all --keep-original-seq")
        sh(cmd1)
        sh(f"cp {wdir}/final-viral-score.tsv {Path(out_dir)/'vs2'/'vs2_all'}/{item}-final-viral-score.tsv")

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done VirSorter2!")


def run_vibrant_single(fasta_list: List[str], out_dir: str, vibrant_db: Path, thread: int, parallel: int, env_single: Path) -> None:
    print("Begin VIBRANT...")
    Path(out_dir, "VIBRANT").mkdir(parents=True, exist_ok=True)

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        wdir = Path(out_dir) / "VIBRANT" / item
        wdir.mkdir(parents=True, exist_ok=True)
        cmd = conda_run(env_single, f"VIBRANT_run.py -t {thread} -d {vibrant_db} -f nucl -virome -i {in_fa} -folder {wdir}")
        sh(cmd)

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done VIBRANT!")


def do_checkV_single(fasta_list: List[str], out_dir: str, checkv_db: Path, thread: int, parallel: int, env_single: Path) -> None:
    print("Begin CheckV...")
    Path(out_dir, "CheckV").mkdir(parents=True, exist_ok=True)

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        wdir = Path(out_dir) / "CheckV" / item
        cmd = conda_run(env_single, f"checkv end_to_end -d {checkv_db} -t {thread} {in_fa} {wdir}")
        sh(cmd)

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done CheckV!")


def do_kaiju_single(fasta_list: List[str], out_dir: str, kaiju_db: Path, thread: int, parallel: int, env_single: Path) -> None:
    print("Begin kaiju...")
    Path(out_dir, "kaiju").mkdir(parents=True, exist_ok=True)
    nodes = kaiju_db / "nodes.dmp"
    names = kaiju_db / "names.dmp"
    fmi   = kaiju_db / "kaiju_db_refseq.fmi"

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        out_tax = Path(out_dir) / "kaiju" / f"{item}.kaiju.tax.out"
        out_name = Path(out_dir) / "kaiju" / f"{item}.kaiju.names.out"
        cmd1 = conda_run(env_single, f"kaiju -t {nodes} -f {fmi} -i {in_fa} -o {out_tax} -z {thread} -v")
        sh(cmd1)
        cmd2 = conda_run(env_single, f"kaiju-addTaxonNames -i {out_tax} -t {nodes} -n {names} -o {out_name} -r superkingdom")
        sh(cmd2)

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done kaiju!")


def run_genomad(fasta_list: List[str], genomad_db: Path, out_dir: str, thread: int, parallel: int, env_dual: Path) -> None:
    print("Begin geNomad...")
    Path(out_dir, "genomad").mkdir(parents=True, exist_ok=True)

    def process_item(item: str):
        in_fa = Path(out_dir) / "filtered_contigs" / f"{item}.contigs.fasta"
        cmd = conda_run(env_dual, f"genomad end-to-end {in_fa} {Path(out_dir)/'genomad'} {genomad_db} -t {thread}")
        sh(cmd)

    with ThreadPoolExecutor(max_workers=parallel) as ex:
        list(ex.map(process_item, fasta_list))
    print("Done geNomad!")


def deal_results(fasta_list: List[str], out_dir: str) -> None:
    dirs = [f"{out_dir}/results", f"{out_dir}/results/dvf", f"{out_dir}/results/vs2",
            f"{out_dir}/results/vibr", f"{out_dir}/results/checkV", f"{out_dir}/results/kaiju", f"{out_dir}/results/genomad"]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)

    print("Integrating key output files...")
    for item in fasta_list:
        sh(f"cp {Path(out_dir,'dvf','dvf_all')}/*  {Path(out_dir,'results','dvf')}")
        sh(f"cp {Path(out_dir,'vs2','vs2_all')}/*  {Path(out_dir,'results','vs2')}")
        sh(f"cp {Path(out_dir,'VIBRANT',item,'VIBRANT_'+item+'.contigs','VIBRANT_results_'+item+'.contigs','VIBRANT_genome_quality_'+item+'.contigs.tsv')}  {Path(out_dir,'results','vibr')}")
        sh(f"cp {Path(out_dir,'CheckV',item,'quality_summary.tsv')}  {Path(out_dir,'results','checkV',item+'.quality_summary.tsv')}")
        sh(f"cp {Path(out_dir,'kaiju',item+'.kaiju.names.out')}  {Path(out_dir,'results','kaiju')}")
        sh(f"cp {Path(out_dir,'genomad',item+'.contigs_marker_classification',item+'.contigs_features.tsv')}  {Path(out_dir,'results','genomad')}")
        sh(f"cp {Path(out_dir,'genomad',item+'.contigs_aggregated_classification',item+'.contigs_aggregated_classification.tsv')}  {Path(out_dir,'results','genomad')}")


def get_scores(out_dir: str) -> None:
    """
    Run 3 rounds of scoring by calling scripts in ./score/ with the current interpreter.
    This does NOT depend on PATH and will use SCORE_DIR to locate the scripts.
    """
    print("Perform three rounds of scoring: single tool, tunning addition and tunning remove.")

    print("Step1: Single tool")
    run_score("VPAC_part2_dvf.py",     Path(out_dir, "results", "dvf"))
    run_score("VPAC_part2_vs2.py",     Path(out_dir, "results", "vs2"))
    run_score("VPAC_part2_vibr.py",    Path(out_dir, "results", "vibr"))
    run_score("VPAC_part2_genomad.py", Path(out_dir, "results", "genomad"))

    print("Step2: Tunning addition")
    run_score(
        "VPAC_part2_rmv.py",
        "-c", Path(out_dir, "results", "checkV"),
        "-v", Path(out_dir, "results", "vs2"),
        "-g", Path(out_dir, "results", "genomad"),
    )

    print("Step3: Tunning remove")
    run_score(
        "VPAC_part2_add.py",
        "-k", Path(out_dir, "results", "kaiju"),
        "-c", Path(out_dir, "results", "checkV"),
        "-v", Path(out_dir, "results", "vs2"),
        "-g", Path(out_dir, "results", "genomad"),
    )


def score_summary(out_dir: str, input_dir: str, fasta_list: List[str]) -> None:
    dirs = [f"{out_dir}/score", f"{out_dir}/summary",
            f"{out_dir}/score/part1", f"{out_dir}/score/part2", f"{out_dir}/score/part3"]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)

    print("Integrating the results of three rounds of scoring...")
    sh(f"cp {Path(out_dir,'results','dvf')}/*.tsv  {Path(out_dir,'score','part1')}")
    sh(f"cp {Path(out_dir,'results','vs2')}/*_vs2.tsv  {Path(out_dir,'score','part1')}")
    sh(f"cp {Path(out_dir,'results','vibr')}/VIBRANT_genome_quality_*_vibrant.tsv  {Path(out_dir,'score','part1')}")
    sh(f"cp {Path(out_dir,'results','genomad')}/*_genomad.tsv  {Path(out_dir,'score','part1')}")

    sh(f"cp {Path(out_dir,'results','kaiju')}/*scores.tsv  {Path(out_dir,'score','part2')}")
    sh(f"cp {Path(out_dir,'results','checkV')}/*_checkV_add.tsv  {Path(out_dir,'score','part2')}")
    sh(f"cp {Path(out_dir,'results','vs2')}/*_vs2_add.tsv  {Path(out_dir,'score','part2')}")
    sh(f"cp {Path(out_dir,'results','genomad')}/*_genomad_add.tsv  {Path(out_dir,'score','part2')}")

    sh(f"cp {Path(out_dir,'results','checkV')}/*_checkV_rmv.tsv  {Path(out_dir,'score','part3')}")
    sh(f"cp {Path(out_dir,'results','vs2')}/*_vs2_rmv.tsv  {Path(out_dir,'score','part3')}")
    sh(f"cp {Path(out_dir,'results','genomad')}/*_genomad_rmv.tsv  {Path(out_dir,'score','part3')}")

    contigs_folder = Path(input_dir)
    score_folders = [Path(out_dir, 'score', 'part1'), Path(out_dir, 'score', 'part2'), Path(out_dir, 'score', 'part3')]

    score_files_info = {
        "part1": ["*.contigs.fasta_gt1bp_dvfpred_dvf.tsv", "*-final-viral-score_vs2.tsv", "VIBRANT_genome_quality_*.contigs_vibrant.tsv", "*.contigs_aggregated_classification_genomad.tsv"],
        "part2": ["*.kaiju.names_kaiju_scores.tsv", "*.quality_summary_checkV_add.tsv", "*-final-viral-score_vs2_add.tsv", "*.contigs_features_genomad_add.tsv"],
        "part3": ["*.quality_summary_checkV_rmv.tsv", "*-final-viral-score_vs2_rmv.tsv", "*.contigs_features_genomad_rmv.tsv"]
    }

    def custom_sort(seq_id: str) -> int:
        m = re.match(r"k\d+_(\d+)", seq_id)
        return int(m.group(1)) if m else 0

    for sample_name in tqdm(fasta_list, desc="Processing samples"):
        contigs_file = contigs_folder / f"{sample_name}.contigs.fa"
        sequence_ids = {line.strip().split()[0][1:] for line in open(contigs_file) if line.startswith('>')}
        sequence_ids = sorted(sequence_ids, key=custom_sort)
        summary_df = pd.DataFrame({"Sequence ID": sequence_ids})

        for step, score_folder in enumerate(tqdm(score_folders, desc=f"Processing steps for {sample_name}", leave=False), start=1):
            score_files_pattern = score_files_info.get(f"part{step}")
            if score_files_pattern is None:
                print("quit")
                continue
            for pattern in score_files_pattern:
                pattern_with_sample = pattern.replace("*", sample_name)
                score_files = [fn for fn in os.listdir(score_folder) if fn.endswith(pattern_with_sample)]
                for score_file in score_files:
                    score_file_path = Path(score_folder, score_file)
                    df = pd.read_csv(score_file_path, sep='\t')
                    df[df.columns[1:]] = df[df.columns[1:]].fillna(0).astype(int)
                    df.columns = [f"{step}_{score_file.split('.')[0]}_score{col}" if idx > 0 else col for idx, col in enumerate(df.columns)]
                    summary_df = summary_df.merge(df, how='left', left_on='Sequence ID', right_on=df.columns[0]).fillna(0)

        summary_df.iloc[:, 1:] = summary_df.iloc[:, 1:].astype(int)
        out_file = Path(out_dir, 'summary', f"{sample_name}_summary.tsv")
        summary_df.to_csv(out_file, index=False, sep='\t')


MODEL_MAP = {
    "CNN": "VPAC-single-CNN.pth",
    "FNN": "VPAC-single-FNN.h5",
    "GB":  "VPAC-single-GB.pkl",
    "KAN": "VPAC-single-KAN.pt",
    "RF":  "VPAC-single-RF.joblib",
    "SVC": "VPAC-single-SVC.joblib",
    "VAE": "VPAC-single-VAE.pth",
    "AE": "VPAC-single-AE.pth",
}


def get_viralseq(input_dir: str, out_dir: str, fasta_list: List[str], model_key: str, models_dir: Path, env_single: Path) -> None:
    Path(out_dir, 'viral_score').mkdir(parents=True, exist_ok=True)
    Path(out_dir, 'labels').mkdir(parents=True, exist_ok=True)
    Path(out_dir, 'mVC').mkdir(parents=True, exist_ok=True)
    Path(out_dir, 'mNVC').mkdir(parents=True, exist_ok=True)

    if model_key not in MODEL_MAP:
        raise ValueError(f"Unknown model '{model_key}'. Valid: {list(MODEL_MAP)}")

    model_path = models_dir / MODEL_MAP[model_key]
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")

    print(f"Predict viral and non-viral contigs by {model_key} model...")
    for item in fasta_list:
        in_summary = Path(out_dir, 'summary', f"{item}_summary.tsv")
        out_score  = Path(out_dir, 'viral_score', f"{item}_score.txt")
        conda_run_py(env_single, "predict.py", in_summary, model_path, out_score)
        all_headers = Path(out_dir, 'viral_score', 'temp_all_headers.txt')
        sh(f"grep '>' {Path(input_dir)/ (item+'.contigs.fa')} > {all_headers}")
        sh(f"awk -F ' ' '{{print $1}}' {all_headers} > {Path(out_dir,'viral_score','temp_ids.txt')}")
        sh(f"sed 's/>//' {Path(out_dir,'viral_score','temp_ids.txt')} > {Path(out_dir,'labels', item+'_all_labels.txt')}")
        sh(f"awk 'NR>1 && $3==\"Virus\"{{print $1}}' {out_score} > {Path(out_dir,'labels', item+'_viral_labels.txt')}")
        sh(conda_run(env_single, f"seqkit grep -w 0 -f {Path(out_dir,'labels', item+'_viral_labels.txt')}  {Path(input_dir, item+'.contigs.fa')} > {Path(out_dir,'mVC', item+'.mVC.fasta')}"))
        sh(conda_run(env_single, f"seqkit grep -w 0 -f {Path(out_dir,'labels', item+'_viral_labels.txt')} -v {Path(input_dir, item+'.contigs.fa')} > {Path(out_dir,'mNVC', item+'.mNVC.fasta')}"))

    print(f"Your viral contigs are in path '{Path(out_dir,'mVC')}/' ")
    print(f"Your non-viral contigs are in path '{Path(out_dir,'mNVC')}/'")


def check_viralseq(out_dir: str) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    from io import StringIO

    for mode in ['mVC', 'mNVC']:
        contigs_dir = Path(out_dir) / mode
        stat_dir = contigs_dir / 'contigs_stat'
        stat_dir.mkdir(parents=True, exist_ok=True)
        fa_files = list(contigs_dir.glob("*.fasta"))
        if not fa_files:
            print(f"No contig files found in {contigs_dir} for statistics.")
            continue
        try:
            files_str = " ".join(map(str, fa_files))
            stats_output = subprocess.check_output(
                f"seqkit stat -a {files_str}", shell=True, universal_newlines=True, stderr=subprocess.STDOUT
            )
            stats_df = pd.read_csv(StringIO(stats_output), sep='\t')
            stats_summary_path = stat_dir / 'contig_stats_summary.tsv'
            stats_df.to_csv(stats_summary_path, sep='\t', index=False)
            print(f"[{mode}] Contig statistics saved to {stats_summary_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error running seqkit stat for {mode}: {e.output}")
            continue
        except Exception as e:
            print(f"Error processing stats for {mode}: {str(e)}")
            continue

        all_samples_data = []
        for file_path in fa_files:
            sample_name = file_path.name.replace('.contigs.fa', '')
            print(f"[{mode}] Processing length distribution for {sample_name}...")
            try:
                cmd = f"seqkit fx2tab -n -l {file_path} | cut -f 2"
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
                lengths = [int(line.strip()) for line in process.stdout if line.strip()]
                bins = list(range(0, 100001, 100))
                counts, bin_edges = np.histogram(lengths, bins=bins)
                all_samples_data.append((sample_name, counts, bin_edges))
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                continue

        if not all_samples_data:
            print(f"[{mode}] No data available for plotting.")
            continue

        n_samples = len(all_samples_data)
        n_cols = 4
        n_rows = (n_samples + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        fig.suptitle(f'Contig Length Distribution ({mode})', fontsize=14)

        for idx, (sample_name, counts, bin_edges) in enumerate(all_samples_data):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col] if n_rows > 1 else axes[col] if n_cols > 1 else axes
            x = bin_edges[:-1]
            ax.bar(x, counts, width=100, align='edge', alpha=0.7, edgecolor='k', linewidth=0.3)
            ax.set_title(sample_name, fontsize=10)
            ax.set_xlabel('Length (bp)', fontsize=8)
            ax.set_ylabel('Count', fontsize=8)
            ax.set_xlim(0, 100000)
            ax.set_xticks(np.arange(0, 100001, 50000))
            ax.tick_params(axis='both', labelsize=8)
        for idx in range(n_samples, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            if n_rows > 1:
                axes[row, col].axis('off')
            else:
                (axes[col] if n_cols > 1 else axes).axis('off')
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        plot_pdf = stat_dir / 'contig_length_distribution.pdf'
        plot_png = stat_dir / 'contig_length_distribution.png'
        fig.savefig(plot_pdf, bbox_inches='tight')
        fig.savefig(plot_png, dpi=500, bbox_inches='tight')
        print(f"[{mode}] Plots saved to:\n  - {plot_pdf}\n  - {plot_png}")
        print("All processes have been done!")


def custom_help():
    console = Console()
    intro = (
        "This pipeline integrates DeepVirFinder, VirSorter2, and VIBRANT to identify viral contigs, "
        "then runs Kaiju and CheckV for taxonomy/quality, and geNomad for discovery.\n"
        "All tool/DB/env paths now come from a YAML config (see --config)."
    )
    examples_md = Markdown(
        "\n**Examples:**\n"
        "```\n"
        "VPAC_predict \\\n"
        "  -i /data/input -o /data/output \\\n"
        "  --min_length 3000 --max_length 50000 \\\n"
        "  --thread 40 --parallel 2 \\\n"
        "  --config /path/to/config.yml \\\n"
        "  --model FNN\n"
        "```\n"
    )
    reads_md = Markdown(
        "**Input FASTA directory (-i):** contains per-sample `*.contigs.fa`."
    )
    out_md = Markdown("**Output directory (-o):** all intermediates and finals will be placed here.")

    console.print(Panel(intro, border_style="cyan", title="VPAC_predict", title_align="left"))
    console.print(examples_md)
    console.print(reads_md)
    console.print(out_md)
    console.print()

    out_tbl = Table(show_header=True, header_style="bold blue")
    out_tbl.add_column("Files/Dirs", style="cyan")
    out_tbl.add_column("Description", style="white")

    out_tbl.add_row("~/outdir/mNVC/", "VPAC single-path classifier predicted non-viral contigs.")
    out_tbl.add_row("~/outdir/mVC/",  "VPAC single-path classifier predicted viral contigs.")

    console.print(Panel(out_tbl, border_style="blue", title="Outputs (overview)", title_align="left"))
    console.print(Panel("Modules run with sensible defaults; tune thresholds via config or CLI.", border_style="cyan"))
    console.print()


    console.print("[bold]Detailed parameters[/bold]\n")
    console.print(Panel("Usage: VPAC_predict [OPTIONS]", border_style="cyan", title="Usage & Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-i, --input_dir", "Directory of input FASTA files (required).")
    g_tbl.add_row("-o, --out_dir", "Directory for all outputs (required).")
    g_tbl.add_row("-lmin, --min_length", "Minimum contig length to analyze (default: 3000).")
    g_tbl.add_row("-lmax, --max_length", "Maximum contig length to analyze (default: 50000).")
    g_tbl.add_row("-t, --thread", "Number of parallel threads (default: 40).")
    g_tbl.add_row("--parallel", "Number of samples to process in parallel (default: 1).")
    g_tbl.add_row("-cf, --config", "Path to config.yml for env/tools/dbs (required).")
    console.print(g_tbl)
    console.print()

    console.print(Panel("Model selection", border_style="blue"))
    m_tbl = Table(show_header=False, box=None, pad_edge=False)
    m_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    m_tbl.add_column("Description", style="white")
    m_tbl.add_row("-m, --model", "Choose one of: CNN, FNN (default), GB, KAN, RF, SVC, VAE, AE")
    console.print(m_tbl)
    console.print()
    console.print(Rule(style="dim"))


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def main():
    parser = CustomArgumentParser(
        description=(
            "Integrates DeepVirFinder, VirSorter2, and VIBRANT to identify viral contigs; "
            "runs Kaiju and CheckV for taxonomy/quality; runs geNomad for viral discovery; "
            "merges scores and exports final viral sequences. Paths come from --config."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    parser.add_argument("-lmin", "--min_length", type=int, default=3000)
    parser.add_argument("-lmax", "--max_length", type=int, default=50000)
    parser.add_argument("-t", "--thread", type=int, default=40)
    parser.add_argument("--parallel", type=int, default=1)
    parser.add_argument("-cf", "--config", type=str, required=True, help="Path to config.yml for env/tools/dbs")
    parser.add_argument("-m", "--model", type=str, choices=list(MODEL_MAP.keys()), default="FNN")

    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    try:
        env_single = Path(cfg["VPAC_single"]).resolve()
        env_dual   = Path(cfg["VPAC_dual"]).resolve()
        vibrant_db = Path(cfg["VIBRANT_db"]).resolve()
        checkv_db  = Path(cfg["CheckV_db"]).resolve()
        kaiju_db   = Path(cfg["Kaiju_db"]).resolve()
        genomad_db = Path(cfg["geNomad_db"]).resolve()
        virsorter2_db = Path(cfg["VirSorter2_db"]).resolve()
        dvf_dir    = Path(cfg["DeepVirFinder"]).resolve()
        models_dir = Path(cfg["VPAC_models"]).resolve()
    except KeyError as e:
        raise KeyError(f"Missing required key in config.yml: {e}")

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    fasta_list = Getfasta_list(args.input_dir)
    run_seqkit(fasta_list, args.input_dir, args.out_dir, args.min_length, args.max_length, args.thread, env_single)

    init_virsorter2_db(env_single, virsorter2_db)

    run_deepvirfinder_single(fasta_list, args.out_dir, args.thread, args.parallel, env_single, dvf_dir)
    run_vibrant_single(fasta_list, args.out_dir, vibrant_db, args.thread, args.parallel, env_single)
    do_checkV_single(fasta_list, args.out_dir, checkv_db, args.thread, args.parallel, env_single)
    do_kaiju_single(fasta_list, args.out_dir, kaiju_db, args.thread, args.parallel, env_single)
    run_genomad(fasta_list, genomad_db, args.out_dir, args.thread, args.parallel, env_dual)
    run_virsorter2_single(fasta_list, args.out_dir, args.thread, args.parallel, env_single)

    deal_results(fasta_list, args.out_dir)
    get_scores(args.out_dir)
    score_summary(args.out_dir, args.input_dir, fasta_list)
    get_viralseq(args.input_dir, args.out_dir, fasta_list, args.model, models_dir, env_single)
    check_viralseq(args.out_dir)


if __name__ == "__main__":
    main()
