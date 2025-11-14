import os
import subprocess
import pandas as pd
import json
import argparse
import shlex
import shutil
import glob
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
import re
from pathlib import Path

def _load_kv_config(cfg_path: str) -> dict:
    """Load a simple key: value YAML (flat) into a dict, ignoring comments and blanks."""
    data = {}
    with open(cfg_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if ':' not in line:
                continue
            key, val = line.split(':', 1)
            key = key.strip()
            val = val.strip()
            if (val.startswith("'") and val.endswith("'")) or (val.startswith('"') and val.endswith('"')):
                val = val[1:-1]
            if ' #' in val:
                val = val.split(' #', 1)[0].strip()
            data[key] = val
    return data

ENV_PREFIX = None
ENV_BIN = None
ENV_RUN_ENVS = None 
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.markdown import Markdown
from rich.rule import Rule


def Getfastq_list(raw_data):
    path = raw_data
    fastq_list = []
    exclude_dirs = {'origin_quality', 'filtered_quality'}
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path) and item not in exclude_dirs:
            fastq_list.append(item)
    
    print('Filtered sample list:', fastq_list)
    return fastq_list

def process_sample_megahit(current_fasta, raw_data, assembly_dir, raw_contigs_dir, 
                          thread, k_min, k_max, k_step, min_length, min_count, 
                          k_list, memory, presets_parameters):
    print(f"Processing the sample: {current_fasta}")
    megahit_command = f"megahit -1 {os.path.join(raw_data, current_fasta, '*_1.fastq')} " \
                      f"-2 {os.path.join(raw_data, current_fasta, '*_2.fastq')} " \
                      f"-o {os.path.join(assembly_dir, current_fasta)} --out-prefix {current_fasta} -t {thread}"
    
    if k_min:
        megahit_command += f" --k-min {k_min}"
    if k_max:
        megahit_command += f" --k-max {k_max}"
    if k_step:
        megahit_command += f" --k-step {k_step}"
    if min_length:
        megahit_command += f" --min-contig-len {min_length}"
    if min_count:
        megahit_command += f" --min-count {min_count}"
    if k_list:
        megahit_command += f" --k-list {k_list}"
    if memory:
        megahit_command += f" --memory {memory}"
    if presets_parameters:
        megahit_command += f" --presets {presets_parameters}"

    try:
        exit_code = subprocess.call(megahit_command, shell=True, env=ENV_RUN_ENVS)
        if exit_code != 0:
            print(f"Error occurred during assembly of {current_fasta}. Command: {megahit_command}")
            return
    except Exception as e:
        print(f"Exception occurred for {current_fasta}: {str(e)}")
        return

    print(f"Sample {current_fasta} has completed assembly.")
    contig_file_path = os.path.join(assembly_dir, current_fasta, f"{current_fasta}.contigs.fa")
    if os.path.exists(contig_file_path):
        shutil.copy(contig_file_path, raw_contigs_dir)
    else:
        print(f"No contig file found for {current_fasta}. Skipping copy.")

def MEGAHIT(out_dir, raw_data, fastq_list, thread, k_min, k_max, k_step, min_length, min_count, k_list, memory, presets_parameters, parallel):
    path_data = os.path.dirname(out_dir)
    raw_contigs_dir = os.path.join(out_dir, "raw_contigs")
    
    if not os.path.exists(raw_contigs_dir):
        os.makedirs(raw_contigs_dir)
    
    print("Performing sequence assembly using MEGAHIT...")

    assembly_dir = os.path.join(out_dir, "Assembly")
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)

    process_func = partial(process_sample_megahit,
                          raw_data=raw_data,
                          assembly_dir=assembly_dir,
                          raw_contigs_dir=raw_contigs_dir,
                          thread=thread,
                          k_min=k_min,
                          k_max=k_max,
                          k_step=k_step,
                          min_length=min_length,
                          min_count=min_count,
                          k_list=k_list,
                          memory=memory,
                          presets_parameters=presets_parameters)

    with Pool(processes=parallel) as pool:
        pool.map(process_func, fastq_list)

    print(f"All your contigs are in '{raw_contigs_dir}'")

    fa_files = glob.glob(os.path.join(raw_contigs_dir, "*.fa"))
    if not fa_files:
        print("No contig files found for statistics.")
        return

    try:
        stats_output = subprocess.check_output(f"seqkit stat -a {' '.join(fa_files)}", shell=True, text=True, stderr=subprocess.STDOUT)
        stats_df = pd.read_csv(StringIO(stats_output), sep='\t')
        stats_summary_path = os.path.join(assembly_dir, 'contig_stats_summary.tsv')
        stats_df.to_csv(stats_summary_path, sep='\t', index=False)
        print(f"Contig statistics summary saved to {stats_summary_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit stat: {e.output}")
        return
    except Exception as e:
        print(f"Error processing stats: {str(e)}")
        return

    all_samples_data = []
    for file_path in fa_files:
        sample_name = os.path.basename(file_path).replace('.contigs.fa', '')
        print(f"Processing length distribution for {sample_name}...")
        try:
            cmd = f"seqkit fx2tab -n -l {file_path} | cut -f 2"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
            lengths = []
            for line in process.stdout:
                line = line.strip()
                if line:
                    lengths.append(int(line))
            bins = list(range(0, 20001, 100))
            counts, bin_edges = np.histogram(lengths, bins=bins)
            all_samples_data.append((sample_name, counts, bin_edges))
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue

    if not all_samples_data:
        print("No data available for plotting.")
        return

    n_samples = len(all_samples_data)
    n_cols = 4
    n_rows = (n_samples + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    fig.suptitle('Contig Length Distribution per Sample', fontsize=14)

    for idx, (sample_name, counts, bin_edges) in enumerate(all_samples_data):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows > 1:
            ax = axes[row, col]
        else:
            ax = axes[col] if n_cols > 1 else axes
        x = bin_edges[:-1] 
        ax.bar(x, counts, width=100, align='edge', edgecolor='k', alpha=0.7, linewidth=0.3)
        ax.set_title(sample_name, fontsize=10)
        ax.set_xlabel('Length (bp)', fontsize=8)
        ax.set_ylabel('Count', fontsize=8)
        ax.set_xlim(0, 20000)
        ax.set_xticks(np.arange(0, 20001, 2500))  
        ax.tick_params(axis='both', labelsize=8)
        ax.ticklabel_format(axis='y', style='plain')

    for idx in range(n_samples, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows > 1:
            axes[row, col].axis('off')
        else:
            if n_cols > 1:
                axes[col].axis('off')
            else:
                axes.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.96])  
    plot_pdf_path = os.path.join(assembly_dir, 'contig_length_distribution.pdf')
    plot_png_path = os.path.join(assembly_dir, 'contig_length_distribution.png')
    fig.savefig(plot_pdf_path, bbox_inches='tight')
    fig.savefig(plot_png_path, dpi=500, bbox_inches='tight')
    plt.close()

    print(f"Contig length distribution plots saved to:\n- {plot_pdf_path}\n- {plot_png_path}")


def process_sample_spades(current_fasta, raw_data, assembly_dir, raw_contigs_dir, 
                         thread, min_length, memory, k_sizes, cov_cutoff, 
                         meta, metaviral, rnaviral, additional_params):
    print(f"Processing the sample: {current_fasta}")
    spades_command = f"spades.py -1 {os.path.join(raw_data, current_fasta, '*_1.fastq')} " \
                    f"-2 {os.path.join(raw_data, current_fasta, '*_2.fastq')} " \
                    f"-o {os.path.join(assembly_dir, current_fasta)} --threads {thread}"

    if memory:
        spades_command += f" -m {memory}"
    if k_sizes:
        spades_command += f" -k {k_sizes}"
    if cov_cutoff:
        spades_command += f" --cov-cutoff {cov_cutoff}"
    if meta:
        spades_command += " --meta"
    if metaviral:
        spades_command += " --metaviral"
    if rnaviral:
        spades_command += " --rnaviral"
    if additional_params:
        spades_command += f" {additional_params}"

    try:
        exit_code = subprocess.call(spades_command, shell=True, env=ENV_RUN_ENVS)
        if exit_code != 0:
            print(f"Error occurred during assembly of {current_fasta}. Command: {spades_command}")
            return
    except Exception as e:
        print(f"Exception occurred for {current_fasta}: {str(e)}")
        return

    print(f"Sample {current_fasta} has completed assembly.")
    input_contig = os.path.join(assembly_dir, current_fasta, 'contigs.fasta')
    output_contig = os.path.join(assembly_dir, current_fasta, f'{current_fasta}.contigs.fa')
    if os.path.exists(input_contig):
        subprocess.call(f"seqkit seq {input_contig} -w 0 -u -m {min_length} -o {output_contig}", shell=True)
        if os.path.exists(output_contig):
            shutil.copy(output_contig, raw_contigs_dir)
        else:
            print(f"Contig processing failed for {current_fasta}")
    else:
        print(f"No contig file found for {current_fasta}")

def SPADES(out_dir, raw_data, fastq_list, thread, min_length, memory, k_sizes, cov_cutoff, meta, metaviral, rnaviral, additional_params, parallel):
    path_data = os.path.dirname(out_dir)
    raw_contigs_dir = os.path.join(out_dir, "raw_contigs")
    
    if not os.path.exists(raw_contigs_dir):
        os.makedirs(raw_contigs_dir)

    print("Performing sequence assembly using SPAdes...")

    assembly_dir = os.path.join(out_dir, "Assembly")
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    process_func = partial(process_sample_spades,
                          raw_data=raw_data,
                          assembly_dir=assembly_dir,
                          raw_contigs_dir=raw_contigs_dir,
                          thread=thread,
                          min_length=min_length,
                          memory=memory,
                          k_sizes=k_sizes,
                          cov_cutoff=cov_cutoff,
                          meta=meta,
                          metaviral=metaviral,
                          rnaviral=rnaviral,
                          additional_params=additional_params)

    with Pool(processes=parallel) as pool:
        pool.map(process_func, fastq_list)

    print(f"All your contigs are in '{raw_contigs_dir}'")

    fa_files = glob.glob(os.path.join(raw_contigs_dir, "*.fa"))
    if not fa_files:
        print("No contig files found for statistics.")
        return

    try:
        stats_output = subprocess.check_output(f"seqkit stat -a {' '.join(fa_files)}", shell=True, text=True, stderr=subprocess.STDOUT)
        stats_df = pd.read_csv(StringIO(stats_output), sep='\t')
        stats_summary_path = os.path.join(assembly_dir, 'contig_stats_summary.tsv')
        stats_df.to_csv(stats_summary_path, sep='\t', index=False)
        print(f"Contig statistics summary saved to {stats_summary_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit stat: {e.output}")
        return
    except Exception as e:
        print(f"Error processing stats: {str(e)}")
        return

    all_samples_data = []
    for file_path in fa_files:
        sample_name = os.path.basename(file_path).replace('.contigs.fa', '')
        print(f"Processing length distribution for {sample_name}...")
        try:
            cmd = f"seqkit fx2tab -n -l {file_path} | cut -f 2"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
            lengths = []
            for line in process.stdout:
                line = line.strip()
                if line:
                    lengths.append(int(line))

            bins = list(range(0, 20001, 100))
            counts, bin_edges = np.histogram(lengths, bins=bins)
            all_samples_data.append((sample_name, counts, bin_edges))
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue

    if not all_samples_data:
        print("No data available for plotting.")
        return

    n_samples = len(all_samples_data)
    n_cols = 4
    n_rows = (n_samples + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    fig.suptitle('Contig Length Distribution per Sample', fontsize=14)

    for idx, (sample_name, counts, bin_edges) in enumerate(all_samples_data):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows > 1:
            ax = axes[row, col]
        else:
            ax = axes[col] if n_cols > 1 else axes
        x = bin_edges[:-1]
        ax.bar(x, counts, width=100, align='edge', edgecolor='k', alpha=0.7, linewidth=0.3)
        ax.set_title(sample_name, fontsize=10)
        ax.set_xlabel('Length (bp)', fontsize=8)
        ax.set_ylabel('Count', fontsize=8)
        ax.set_xlim(0, 20000)
        ax.set_xticks(np.arange(0, 20001, 2500)) 
        ax.tick_params(axis='both', labelsize=8)
        ax.ticklabel_format(axis='y', style='plain')

    for idx in range(n_samples, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows > 1:
            axes[row, col].axis('off')
        else:
            if n_cols > 1:
                axes[col].axis('off')
            else:
                axes.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.96])  

    plot_pdf_path = os.path.join(assembly_dir, 'contig_length_distribution.pdf')
    plot_png_path = os.path.join(assembly_dir, 'contig_length_distribution.png')
    fig.savefig(plot_pdf_path, bbox_inches='tight')
    fig.savefig(plot_png_path, dpi=500, bbox_inches='tight')
    plt.close()

    print(f"Contig length distribution plots saved to:\n- {plot_pdf_path}\n- {plot_png_path}")


def rename_contig_headers(out_dir):
    raw_contigs_dir = os.path.join(out_dir, "raw_contigs")
    raw_contigs_new_dir = os.path.join(out_dir, "raw_contigs_new")
    
    os.makedirs(raw_contigs_new_dir, exist_ok=True)
    
    header_pattern = re.compile(r"^>NODE_(\d+).*_type_([a-z]+)")
    used_ids = set() 
    for fa_path in Path(raw_contigs_dir).glob("*.fa"):
        new_path = Path(raw_contigs_new_dir) / fa_path.name
        
        with open(fa_path, "r") as infile, open(new_path, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    match = header_pattern.match(line.strip())
                    if match:
                        node_num = match.group(1)
                        contig_type = match.group(2).lower()
                        suffix = "0001" if contig_type == "circular" else "0002"
                        new_header = f">k141_{node_num}{suffix}"
                        
                        count = 1
                        original_header = new_header
                        while new_header[1:] in used_ids:
                            new_header = f"{original_header[:-4]}{str(count).zfill(4)}"
                            count += 1
                        
                        used_ids.add(new_header[1:])
                        outfile.write(new_header + "\n")
                    else:
                        outfile.write(line) 
                else:
                    outfile.write(line)
    

def custom_help():
    console = Console()

    intro = (
        "The 'Assembly' script provides two assembly modules: MEGAHIT and SPAdes. "
        "It assembles quality-controlled FASTQ files into contigs and summarizes assembly statistics."
    )

    examples_md = Markdown(
        "\n**Examples:**\n"
        "```\n"
        "Assembly  -r /data  -out /results  --tool spades   -p 80  --meta  --spades-memory 1024  --cov-cutoff  auto\n"
        "Assembly  -r /data  -out /results  --tool megahit  -p 80   -kmi 21  -kma 141  -ks 10  --min-count 5  --memory 0.8\n"
        "```\n"
    )

    reads_md = Markdown(
        "**Input reads location (-r):** directory that contains per-sample subfolders with decompressed paired-end FASTQ files, e.g.\n"
        "```\n"
        "/xxxx/Clean_reads/A1/A1_1.fastq\n"
        "/xxxx/Clean_reads/A1/A1_2.fastq\n"
        "/xxxx/Clean_reads/A2/A2_1.fastq\n"
        "/xxxx/Clean_reads/A2/A2_2.fastq\n"
        "...\n"
        "```\n"
        "Therefore, `-r` should be set to `/xxxx/Clean_reads`.\n"
    )

    out_md = Markdown(
        "**Output directory (-out):** the root folder where assembled contigs and reports will be written."
    )

    console.print(Panel(intro, border_style="cyan", title="Assembly", title_align="left"))
    console.print(examples_md)
    console.print(reads_md)
    console.print(out_md)
    console.print()

    out_tbl = Table(show_header=True, header_style="bold blue")
    out_tbl.add_column("Files/Dirs", style="cyan", no_wrap=True)
    out_tbl.add_column("Description", style="white")
    out_tbl.add_row("~/raw_contigs/A1_contigs.fa & ~/raw_contigs/A1_contigs.fasta",
                    "Assembled contigs (per sample).")
    out_tbl.add_row("~/Assembly/contig_stats_summary.tsv",
                    "Summary statistics: num_seqs, sum_len, min_len, avg_len, max_len, GC(%).")
    out_tbl.add_row("~/Assembly/contig_length_distribution.pdf & ~/Assembly/contig_length_distribution.png",
                    "Histogram of contig length distribution.")
    console.print(Panel(out_tbl, border_style="blue", title="Outputs", title_align="left"))
    console.print(Panel("Both tools work out-of-the-box with sensible defaults; tune parameters as needed.",
                        border_style="cyan"))
    console.print()

    # ====================== Detailed parameters ===============================
    console.print("[bold]Detailed parameters[/bold]\n")

    # --- Usage & Global
    console.print(Panel("Usage: Assembly [OPTIONS] INPUT OUTPUT",
                        border_style="cyan", title="Usage & Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-cf, --config", "DPath to config.yml for env/tools/dbs.")
    g_tbl.add_row("-r, --raw_data", "Directory containing per-sample subdirectories of FASTQ files.")
    g_tbl.add_row("-out, --out_dir", "Output directory for assembly results.")
    g_tbl.add_row("-p, --thread", "Number of threads (default: 80).")
    g_tbl.add_row("-m, --min_length", "Minimum contig length to keep (default: 500).")
    g_tbl.add_row("--tool", "Assembly backend: 'megahit' or 'spades' (required).")
    g_tbl.add_row("--parallel", "Number of samples to process in parallel (default: 2).")
    console.print(g_tbl)
    console.print()

    # --- MEGAHIT --------------------------------------------------------------
    console.print(Panel("MEGAHIT Parameters", border_style="magenta"))
    m_tbl1 = Table(show_header=False, box=None, pad_edge=False)
    m_tbl1.add_column("Flag", style="bold cyan", no_wrap=True)
    m_tbl1.add_column("Description", style="white")
    console.print("[bold]K-mer schedule:[/bold]")
    m_tbl1.add_row("-kmi, --k_min", "Minimum k-mer size (default: 21).")
    m_tbl1.add_row("-kma, --k_max", "Maximum k-mer size (default: 141).")
    m_tbl1.add_row("-ks,  --k_step", "Increment for k-mer sizes (default: 10).")
    m_tbl1.add_row("--k-list", "Explicit comma-separated k list (e.g., 21,29,39).")
    console.print(m_tbl1)

    m_tbl2 = Table(show_header=False, box=None, pad_edge=False)
    m_tbl2.add_column("Flag", style="bold cyan", no_wrap=True)
    m_tbl2.add_column("Description", style="white")
    console.print()
    console.print("[bold]Graph / memory / presets:[/bold]")
    m_tbl2.add_row("--min-count", "Minimum k-mer multiplicity for filtering (default: 2).")
    m_tbl2.add_row("--memory", "Max memory fraction for construction (default: 0.9).")
    m_tbl2.add_row("--presets", "Override a group of parameters (MEGAHIT presets).")
    console.print(m_tbl2)
    console.print()

    # --- SPAdes ---------------------------------------------------------------
    console.print(Panel("SPAdes Parameters", border_style="green"))
    s_tbl1 = Table(show_header=False, box=None, pad_edge=False)
    s_tbl1.add_column("Flag", style="bold cyan", no_wrap=True)
    s_tbl1.add_column("Description", style="white")
    console.print("[bold]Core settings:[/bold]")
    s_tbl1.add_row("--spades-memory", "RAM limit in GB (default: 250).")
    s_tbl1.add_row("--k-sizes", "Odd k-mer sizes < 128 (comma-separated).")
    s_tbl1.add_row("--cov-cutoff", "Coverage cutoff (float), or 'auto' / 'off'.")
    console.print(s_tbl1)

    s_tbl2 = Table(show_header=False, box=None, pad_edge=False)
    s_tbl2.add_column("Flag", style="bold cyan", no_wrap=True)
    s_tbl2.add_column("Description", style="white")
    console.print()
    console.print("[bold]Pipelines / extras:[/bold]")
    s_tbl2.add_row("--meta", "Use metagenomic mode.")
    s_tbl2.add_row("--metaviral", "Run metaviralSPAdes for virus detection.")
    s_tbl2.add_row("--rnaviral", "RNA virus assembly from RNA-Seq data.")
    s_tbl2.add_row("--spades-params", "Additional SPAdes parameters (single string).")
    console.print(s_tbl2)

    console.print()
    console.print(Rule(style="dim"))



class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="Metagenomic assembly of paired-end reads into contigs using MEGAHIT or SPAdes, with summary statistics.")
    parser.add_argument('-cf','--config', type=str, help='Path to config.yml for env/tools/dbs')
    parser.add_argument("--parallel", type=int, default=2)
    parser.add_argument("-r", "--raw_data", type=str, required=True)
    parser.add_argument("-out", "--out_dir", type=str, required=True)
    parser.add_argument("-p", "--thread", type=int, default=80)
    parser.add_argument("-m", "--min_length", type=int, default=500)
    parser.add_argument("--tool", choices=['megahit', 'spades'], required=True)
    parser.add_argument("-kmi", "--k_min", type=int, default=21)
    parser.add_argument("-kma", "--k_max", type=int, default=141)
    parser.add_argument("-ks", "--k_step", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=2)
    parser.add_argument("--k-list", type=str)
    parser.add_argument("--memory", type=float, default=0.9)
    parser.add_argument("--presets", type=str)
    parser.add_argument("--spades-memory", type=int, default=250)
    parser.add_argument("--k-sizes", type=str)
    parser.add_argument("--cov-cutoff", type=float)
    parser.add_argument("--meta", action='store_true')
    parser.add_argument("--metaviral", action='store_true')
    parser.add_argument("--rnaviral", action='store_true')
    parser.add_argument("--spades-params", type=str)
    args = parser.parse_args()


    if args.config:
        cfg = _load_kv_config(args.config)
        vmp_env = cfg.get('VMP_env')
        if not vmp_env:
            raise SystemExit("[Assembly] config.yml missing 'VMP_env' key")

        global ENV_PREFIX, ENV_BIN, ENV_RUN_ENVS
        ENV_PREFIX = vmp_env.rstrip('/')
        ENV_BIN = os.path.join(ENV_PREFIX, 'bin')

        if not os.path.isdir(ENV_BIN):
            raise SystemExit(f"[Assembly] ENV bin not found: {ENV_BIN}")

        ENV_RUN_ENVS = os.environ.copy()
        ENV_RUN_ENVS['PATH'] = ENV_BIN + os.pathsep + ENV_RUN_ENVS.get('PATH', '')


    else:

        raise SystemExit("[Assembly] Please provide --config to point to config.yml containing VMP_env")

    fastq_list = Getfastq_list(args.raw_data)

    if args.tool == 'megahit':
        MEGAHIT(args.out_dir, args.raw_data, fastq_list, args.thread, args.k_min, args.k_max, args.k_step, 
                args.min_length, args.min_count, args.k_list, args.memory, args.presets, args.parallel)
    elif args.tool == 'spades':
        SPADES(args.out_dir, args.raw_data, fastq_list, args.thread, args.min_length, args.spades_memory, 
               args.k_sizes, args.cov_cutoff, args.meta, args.metaviral, args.rnaviral, args.spades_params, args.parallel)
        rename_contig_headers(args.out_dir)
        
if __name__ == "__main__":
    main()