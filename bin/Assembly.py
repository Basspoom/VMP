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
        exit_code = subprocess.call(megahit_command, shell=True)
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
        exit_code = subprocess.call(spades_command, shell=True)
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
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "The 'Assembly' script contains two assemble modules, namely MEGAHIT and SPAdes. This script is used for assembly. It can assemble the quality-controlled fastq files to obtain contigs." + RESET)

    print("\n" + PURPLE + "Examples:" + RESET)
    print("  Assembly  -r /data  -out /results  --tool spades -p 80  --meta  --spades-memory 1024  --cov-cutoff auto")
    print("  Assembly  -r /data  -out /results  --tool megahit  -p 80  -kmi 21  -kma 141  -ks 10  --min-count 5  --memory 0.8 ")

    print("\n" + GREEN + "The -r parameter points to the directory containing the paired-end sequencing FASTQ files (please decompress first), for example:" + RESET)
    print("  " + BLUE + "/xxxx/Clean_reads/A1/A1_1.fastq" + RESET)
    print("  " + BLUE + "/xxxx/Clean_reads/A1/A1_2.fastq" + RESET)
    print("  " + BLUE + "/xxxx/Clean_reads/A2/A2_1.fastq" + RESET)
    print("  " + BLUE + "/xxxx/Clean_reads/A2/A2_2.fastq" + RESET)
    print("  ......")
    print("Therefore, -r should be set to /xxxx/Clean_reads.")

    print("\n" + GREEN + "The -out parameter specifies the output directory for quality control and assembly." + RESET)

    print("\n" + BLUE + "  ~/raw_contigs/A1_contigs.fa" + RESET + " & " + BLUE + "~/raw_contigs/A1_contigs.fasta" + RESET + ": Contigs obtained from assembly.")
    print("  " + BLUE + "~/Assembly/contig_stats_summary.tsv" + RESET + ": Statistics of assembled contigs, including indicators such as num_seqs, sum_len, min_len, avg_len, max_len, and GC(%).")
    print("  " + BLUE + "~/Assembly/contig_length_distribution.pdf" + RESET + " & " + BLUE + "~/Assembly/contig_length_distribution.png" + RESET + ": Histogram of the length distribution of assembled contigs.")

    print("\n" + GREEN + "Both tools can be analyzed using preset parameters, and users can adjust parameters according to their data needs." + RESET)


    print("\n" + PURPLE + "Detailed parameters")
    print(" " + GREEN + "=" * 50 + " " + RESET + YELLOW + "Usage & Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)

    print(YELLOW + "Global parameters:" + RESET)
    print(f"  {'-r, --raw_data':<40} Directory containing subdirectories for single-sample data.")
    print(f"  {'-out, --out_dir':<40} Output directory for quality control and assembly.")
    print(f"  {'-p, --thread':<40} Number of parallel threads (default: 80).")
    print(f"  {'-m, --min_length':<40} Min length of assembly contigs (Default: 500).")
    print(f"  {'--tool':<40} Assembly tool to use: 'megahit' or 'spades'.")
    print(f"  {'--parallel':<40} Number of samples to process in parallel (default: 1).")

    print("\n" + GREEN + "=" * 51 + " " + RESET + YELLOW + "Assembly Tool Parameters" + RESET + " " + GREEN + "=" * 51 + RESET)

    print("\n" + YELLOW + "MEGAHIT Parameters:" + RESET)
    print(f"  {'-kmi, --k_min':<40} Minimum k-mer size (default: 21).")
    print(f"  {'-kma, --k_max':<40} Maximum k-mer size (default: 141).")
    print(f"  {'-ks, --k_step':<40} Increment of k-mer size (default: 10).")
    print(f"  {'--min-count':<40} Minimum multiplicity for filtering (default: 2).")
    print(f"  {'--k-list':<40} Comma-separated list of k-mer sizes (e.g. 21,29,39).")

    print(f"  {'--memory':<40} Max memory for construction (fraction of total memory, default: 0.9).")
    print(f"  {'--presets':<40} Override a group of parameters.")

    print("\n" + YELLOW + "SPAdes Parameters:" + RESET)
    print(f"  {'--spades-memory':<40} RAM limit for SPAdes in Gb (default: 250).")
    print(f"  {'--k-sizes':<40} List of k-mer sizes (must be odd and less than 128).")
    print(f"  {'--cov-cutoff':<40} Coverage cutoff value (positive float, 'auto', or 'off').")
    print(f"  {'--meta':<40} Use this flag for metagenomic data.")
    print(f"  {'--metaviral':<40} Run metaviralSPAdes pipeline for virus detection.")
    print(f"  {'--rnaviral':<40} Enable virus assembly module from RNA-Seq data.")
    print(f"  {'--spades-params':<40} Additional SPAdes parameters as a single string.")



class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="Upstream of metagenomics (quality control and contigs assembly).")
    parser.add_argument("--parallel", type=int, default=2, help="Number of samples to process in parallel (default: 1)")

    # Basic parameters
    parser.add_argument("-r", "--raw_data", type=str, required=True, help="Directory containing subdirectories for single-sample data.")
    parser.add_argument("-out", "--out_dir", type=str, required=True, help="Output directory for quality control and assembly.")
    parser.add_argument("-p", "--thread", type=int, default=80, help="Number of parallel threads (default: 80)")
    parser.add_argument("-m", "--min_length", type=int, default=500, help="Min length of assembly contigs (Default: 500)")
    parser.add_argument("--tool", choices=['megahit', 'spades'], required=True, help="Assembly tool to use: 'megahit' or 'spades'.")

    # megahit parameters
    parser.add_argument("-kmi", "--k_min", type=int, default=21, help="Minimum k-mer size (Default: 21)")
    parser.add_argument("-kma", "--k_max", type=int, default=141, help="Maximum k-mer size (Default: 141)")
    parser.add_argument("-ks", "--k_step", type=int, default=10, help="Increment of k-mer size (Default: 10)")
    parser.add_argument("--min-count", type=int, default=2, help="Minimum multiplicity for filtering (Default: 2)")
    parser.add_argument("--k-list", type=str, help="Comma-separated list of k-mer sizes (e.g. 21,29,39)")
    parser.add_argument("--memory", type=float, default=0.9, help="Max memory for construction (fraction of total memory, Default: 0.9)")
    parser.add_argument("--presets", type=str, help="Override a group of parameters.")

    # SPAdes parameters
    parser.add_argument("--spades-memory", type=int, default=250, help="RAM limit for SPAdes in Gb (default: 250)")
    parser.add_argument("--k-sizes", type=str, help="List of k-mer sizes (must be odd and less than 128, e.g. '21,33,55')")
    parser.add_argument("--cov-cutoff", type=float, help="Coverage cutoff value (a positive float number, or 'auto', or 'off')")
    parser.add_argument("--meta", action='store_true', help="Use this flag for metagenomic data.")
    parser.add_argument("--metaviral", action='store_true', help="Run metaviralSPAdes pipeline for virus detection.")
    parser.add_argument("--rnaviral", action='store_true', help="Enable virus assembly module from RNA-Seq data.")
    parser.add_argument("--spades-params", type=str, help="Additional SPAdes parameters as a single string (e.g. --careful --only-assembler).")


    args = parser.parse_args()

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