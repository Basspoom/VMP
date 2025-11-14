import os
import sys
import re
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.markdown import Markdown

def _load_vmp_env_prefix(cfg_path: str) -> str:
    import yaml
    with open(cfg_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    candidates = []
    for key in ("VMP_env",):
        if key in data:
            candidates.append(data.get(key))
    for lvl in ("env", "environments"):
        if isinstance(data.get(lvl), dict) and "VMP_env" in data[lvl]:
            candidates.append(data[lvl]["VMP_env"])

    for c in candidates:
        if isinstance(c, str) and c.strip():
            return c.strip()

    raise KeyError("VMP_env not found in config.yml (tried: VMP_env, env.VMP_env, environments.VMP_env).")


class EnvRunner:
    def __init__(self, conda_prefix: str | None):
        self.conda_prefix = conda_prefix

    def run(self, cmd: str) -> int:
        if self.conda_prefix:
            wrapped = f'conda run --no-capture-output --prefix {self.conda_prefix} bash -lc "{cmd}"'
        else:
            wrapped = cmd
        return subprocess.call(wrapped, shell=True)


def Getfasta_list(input_dir):
    fasta_list = []
    exclude_dirs = {'contigs_stat'}
    for item in os.listdir(input_dir):
        item_path = os.path.join(input_dir, item)
        if os.path.isfile(item_path) and item not in exclude_dirs:
            match = re.search(r'(.*?)\.\w+\.fasta', item)
            if match:
                fasta_list.append(match.group(1))
    print('Your sample list:', fasta_list)
    return fasta_list


def do_cdhit_ANI(fasta_list, input_dir, out_dir, thread, drep_ani, memory, global_identity, band_width, accurate_mode, sort_by_size, sort_fasta_by_size, parallel, runner: EnvRunner | None = None):
    os.makedirs(os.path.join(out_dir, "vOTU/cd-hit_drp"), exist_ok=True)
    print("Using CD-HIT to remove redundancy...")
    print(f"Average Nucleotide Identity = {drep_ani}")

    def process_cdhit(item):
        cdhit_command = f"cd-hit -i {os.path.join(input_dir, item + '.*.fasta')} -o {os.path.join(out_dir, 'vOTU/cd-hit_drp', item + '.drp.fasta')} -T {thread}"
        if drep_ani:
            cdhit_command += f" -c {drep_ani}"
        if memory:
            cdhit_command += f" -M {memory}"
        if global_identity:
            cdhit_command += f" -G {global_identity}"
        if band_width:
            cdhit_command += f" -b {band_width}"
        if accurate_mode:
            cdhit_command += f" -g {accurate_mode}"
        if sort_by_size:
            cdhit_command += f" -sc {sort_by_size}"
        if sort_fasta_by_size:
            cdhit_command += f" -sf {sort_fasta_by_size}"
        exit_code = runner.run(cdhit_command)
        if exit_code != 0:
            print(f"Error occurred during CD-HIT processing for {item}. Command: {cdhit_command}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_cdhit, fasta_list)

    print("Done removing redundancy!")


def do_mmseqs(fasta_list, out_dir, clustering, seed_sub_mat, sensitivity, kmer_length, target_search_mode, max_seqs, cov_mode, alignment_mode, cluster_mode, parallel, runner: EnvRunner | None = None):
    os.makedirs(os.path.join(out_dir, "vOTU/cluster"), exist_ok=True)
    print("Using MMseqs2 to cluster viral contigs...")
    print(f"Alignment fraction = {clustering}")

    def process_mmseqs(item):
        item_db_dir = os.path.join(out_dir, "vOTU/cluster/vOTU", item)
        os.makedirs(item_db_dir, exist_ok=True)

        command1 = f"mmseqs createdb {os.path.join(out_dir, 'vOTU/cd-hit_drp', item + '.drp.fasta')} {item_db_dir}/{item}-contigsDB"
        runner.run(command1)

        command = f"mmseqs cluster {item_db_dir}/{item}-contigsDB  {out_dir}/vOTU/cluster/vOTU/{item}.cluster  tmp  --min-seq-id {clustering} "
        if seed_sub_mat:
            command += f" --seed-sub-mat {seed_sub_mat}"
        if sensitivity:
            command += f" -s {sensitivity}"
        if kmer_length:
            command += f" -k {kmer_length}"
        if target_search_mode:
            command += f" --target-search-mode {target_search_mode}"
        if max_seqs:
            command += f" --max-seqs {max_seqs}"
        if cov_mode is not None:
            command += f" --cov-mode {cov_mode}"
        if alignment_mode:
            command += f" --alignment-mode {alignment_mode}"
        if cluster_mode is not None:
            command += f" --cluster-mode {cluster_mode}"
        runner.run(command)

        os.makedirs(os.path.join(out_dir, "vOTU/cluster/vOTU", item), exist_ok=True)

        runner.run(f"mmseqs createtsv {item_db_dir}/{item}-contigsDB {out_dir}/vOTU/cluster/vOTU/{item}.cluster {out_dir}/vOTU/cluster/{item}.cluster.tsv")
        runner.run(f"mmseqs createseqfiledb {item_db_dir}/{item}-contigsDB  {out_dir}/vOTU/cluster/vOTU/{item}.cluster  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-DB_clu_seq")
        runner.run(f"mmseqs result2flat {item_db_dir}/{item}-contigsDB  {item_db_dir}/{item}-contigsDB  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-DB_clu_seq  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-clu_seq.fasta")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_mmseqs, fasta_list)

    print("Done clustering!")


def _edges_from_skani_triangle(af_path: str, fasta_path: str, edges_path: str, af_threshold_frac: float):
    ids = []
    with open(fasta_path, 'r') as f:
        for ln in f:
            if ln.startswith('>'):
                ids.append(ln[1:].strip())
    n = len(ids)
    if n == 0:
        open(edges_path, 'w').close()
        return

    tri = [[None] * n for _ in range(n)]
    with open(af_path, 'r') as f:
        first = f.readline().strip()
        try:
            _ = int(first)
        except Exception:
            f.seek(0)

        row = 0
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            _name = parts[0] if parts else ""
            values = [x for x in parts[1:] if x != ""]
            for col in range(len(values)):
                try:
                    val = float(values[col])
                except Exception:
                    val = 0.0
                tri[row][col] = val
                tri[col][row] = val
            row += 1
            if row >= n:
                break

    for i in range(n):
        tri[i][i] = 100.0

    thr = af_threshold_frac * 100.0  
    with open(edges_path, 'w') as w:
        for i in range(n):
            for j in range(i, n):
                af = tri[i][j] if tri[i][j] is not None else 0.0
                if af >= thr:
                    wt = af / 100.0
                    w.write(f"{ids[i]}\t{ids[j]}\t{wt:.4f}\n")
        for id_ in ids:
            w.write(f"{id_}\t{id_}\t1.0000\n")


def clean_and_format_fasta(fasta_list, out_dir, parallel, runner: EnvRunner | None = None):
    def process_clean(item):
        sequences = {}
        with open(os.path.join(out_dir, f"vOTU/cluster/vOTU/{item}", f"{item}-clu_seq.fasta"), 'r') as file:
            identifier = None
            sequence = []
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if identifier and ''.join(sequence) not in sequences:
                        sequences[identifier] = ''.join(sequence)
                    identifier = line[1:]
                    sequence = []
                else:
                    sequence.append(line)
            if identifier and ''.join(sequence) not in sequences:
                sequences[identifier] = ''.join(sequence)

        with open(os.path.join(out_dir, f"{item}.v.fasta"), 'w') as outfile:
            for identifier, seq in sequences.items():
                outfile.write(f">{identifier}\n{seq}\n")

        runner.run(f"seqkit seq -m 1 --only-id --upper-case -w 0 {os.path.join(out_dir, f'{item}.v.fasta')} -o {os.path.join(out_dir, f'{item}.vOTU.fasta')}")
        runner.run(f"rm {os.path.join(out_dir, f'{item}.v.fasta')}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_clean, fasta_list)


def pyleiden(fasta_list, input_dir, out_dir, thread, drep_ani, clustering, sparse, ci, detailed, diagonal, distance, min_af, full_matrix,
             marker_c, compression_factor, similarity_threshold, resolution, objective_function, beta, number_iterations, seed, quiet, parallel, runner: EnvRunner | None = None):

    os.makedirs(os.path.join(out_dir, "vOTU/pyleiden"), exist_ok=True)
    print("Using aksni and pyleiden to remove redundancy and cluster viral contigs ...")
    print(f"Average Nucleotide Identity = {drep_ani}")
    print(f"Alignment fraction = {clustering}")

    def process_pyleiden(item):
        sample_dir = os.path.join(out_dir, "vOTU/pyleiden", item)
        os.makedirs(sample_dir, exist_ok=True)

        viral_fa = os.path.join(sample_dir, f"{item}.viral_contigs.fasta")
        runner.run(
            f'seqkit seq --only-id --upper-case -w 0 {os.path.join(input_dir, item + ".*.fasta")} -o "{viral_fa}"'
        )

        skani_cmd = f'cd "{sample_dir}" && skani triangle -i -t {thread} "{viral_fa}"'
        if sparse:                skani_cmd += " --sparse"
        if ci:                    skani_cmd += " --ci"
        if detailed:              skani_cmd += " --detailed"
        if diagonal:              skani_cmd += " --diagonal"
        if distance:              skani_cmd += " --distance"
        if min_af is not None:    skani_cmd += f" --min-af {min_af}"
        if full_matrix:           skani_cmd += " --full-matrix"
        if marker_c is not None:  skani_cmd += f" -m {marker_c}"
        if compression_factor is not None: skani_cmd += f" -c {compression_factor}"
        if similarity_threshold is not None: skani_cmd += f" -s {similarity_threshold}"
        exit_code = runner.run(skani_cmd)
        if exit_code != 0:
            print(f"Error occurred during SKANI processing for {item}. Command: {skani_cmd}")

        af_path   = os.path.join(sample_dir, "skani_matrix.af")
        edges_out = os.path.join(sample_dir, f"{item}.mVC.skani_edges.tsv")
        if not os.path.exists(af_path):
            alt = os.path.join(sample_dir, f"{item}.mVC.skani_output.tsv")
            if os.path.exists(alt):
                af_path = alt
        _edges_from_skani_triangle(af_path, viral_fa, edges_out, clustering)

        clusters_out = os.path.join(sample_dir, f"{item}.mVC.clusters.txt")
        pyleiden_command = f'pyleiden "{edges_out}" "{clusters_out}"'
        if resolution:             pyleiden_command += f" -r {resolution}"
        if objective_function:     pyleiden_command += f" -o {objective_function}"
        if beta:                   pyleiden_command += f" -b {beta}"
        if number_iterations is not None: pyleiden_command += f" -n {number_iterations}"
        if seed:                   pyleiden_command += f" -s {seed}"
        if quiet:                  pyleiden_command += " -q"
        exit_code = runner.run(pyleiden_command)
        if exit_code != 0:
            print(f"Error occurred during pyleiden processing for {item}. Command: {pyleiden_command}")

        labels_out = os.path.join(sample_dir, f"{item}.mVC.labels.txt")
        runner.run(f'cut -f1 "{clusters_out}" > "{labels_out}"')
        runner.run(
            f'seqkit grep -w 0 -f "{labels_out}" "{viral_fa}" > {os.path.join(out_dir, f"{item}.vOTU.fasta")}'
        )

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_pyleiden, fasta_list)


def custom_help():
    console = Console()

    intro = (
        "This script performs de-replication (by ANI) and clustering (by AF) to generate high-quality vOTUs from viral contigs. "
        "You can choose between two pipelines:\n"
        "1) cd-hit + MMseqs2, or\n"
        "2) skani + pyleiden."
    )
    examples_md = Markdown(
        "\nExample:\n"
        "```\n"
        "vOTU -i /data/input -o /data/output --tool pyleiden --drep_ani 0.95 --clustering 0.85 -t 40 \\\n"
        "     -cf /path/to/config.yml\n"
        "```\n"
    )
    io_md = Markdown(
        "Input (-i): directory containing input FASTA files, e.g.\n"
        "```\n"
        "/path/to/input/sample1.fasta\n"
        "/path/to/input/sample2.fasta\n"
        "```\n"
        "Output (-o): directory for results.\n\n"
        "--drep_ani: de-replication threshold by ANI (default: 0.95)\n\n"
        "--clustering: alignment fraction (AF) threshold for clustering (default: 0.85)\n\n"
        "-t / --thread: number of threads (default: 40)\n"
        "-cf / --config: Path to config.yml for env/tools/dbs.\n"
        "--tool : 'cd-Hit + MMseqs2' or 'skani + pyLeiden'.\n"
    )

    console.print(Panel(intro, border_style="cyan", title="vOTU (De-replication & Clustering)", title_align="left"))
    console.print(examples_md)
    console.print(io_md)
    console.print()

    console.print("[bold]Detailed parameters[/bold]\n")
    console.print(Panel("Usage: vOTU [OPTIONS] -i IN_DIR -o OUT_DIR --tool {cd-hit_mmseq,skani_pyleiden}",
                        border_style="cyan", title="Usage & Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-i, --input_dir", "Directory containing input FASTA files.")
    g_tbl.add_row("-o, --out_dir", "Output directory for results.")
    g_tbl.add_row("-ani, --drep_ani", "De-replication ANI threshold (default: 0.95).")
    g_tbl.add_row("-cls, --clustering", "Alignment fraction (AF) for clustering (default: 0.85).")
    g_tbl.add_row("-t, --thread", "Number of threads (default: 40).")
    g_tbl.add_row("--tool", "Choose pipeline: 'cd-hit_mmseq' or 'skani_pyleiden'.")
    g_tbl.add_row("--parallel", "Number of samples to process in parallel (default: 1).")
    g_tbl.add_row("-cf, --config", "Path to config.yml.")
    console.print(g_tbl)
    console.print()

    console.print(Panel("Choose 1: cd-hit & MMseqs2", border_style="magenta"))
    t_cdhit = Table(show_header=True, header_style="bold magenta")
    t_cdhit.add_column("CD-HIT", style="cyan", no_wrap=True)
    t_cdhit.add_column("Description", style="white")
    t_cdhit.add_row("-m, --memory", "Memory limit (MB).")
    t_cdhit.add_row("-g, --global_identity", "Use global identity (0 local, 1 global; default: 1).")
    t_cdhit.add_row("-b, --band_width", "Band width of alignment (default: 20).")
    t_cdhit.add_row("-g_mode, --accurate_mode", "Accurate mode (0 fast, 1 accurate; default: 0).")
    t_cdhit.add_row("-sc, --sort_by_size", "Sort clusters by size (0/1; default: 0).")
    t_cdhit.add_row("-sf, --sort_fasta_by_size", "Sort FASTA by cluster size (0/1; default: 0).")

    t_mmseq = Table(show_header=True, header_style="bold magenta")
    t_mmseq.add_column("MMseqs2", style="cyan", no_wrap=True)
    t_mmseq.add_column("Description", style="white")
    t_mmseq.add_row("--seed-sub-mat", "Substitution matrix file for k-mer generation.")
    t_mmseq.add_row("-s, --sensitivity", "Sensitivity (1.0 faster; 4.0 fast; 7.5 sensitive; default: 4.0).")
    t_mmseq.add_row("-k, --kmer_length", "k-mer length (0 auto-optimal).")
    t_mmseq.add_row("--target-search-mode", "0 regular k-mer; 1 similar k-mer (default: 0).")
    t_mmseq.add_row("--max-seqs", "Max results per query passing prefilter (default: 20).")
    t_mmseq.add_row("--cov-mode", "Coverage mode (0: coverage of query & target; default: 0).")
    t_mmseq.add_row("--alignment-mode", "How to compute the alignment (default: 3).")
    t_mmseq.add_row("--cluster-mode", "0 Set-Cover; 1 Connected component (default: 0).")

    console.print(t_cdhit)
    console.print(t_mmseq)
    console.print()

    console.print(Panel("Choose 2: skani & pyleiden", border_style="green"))
    t_skani = Table(show_header=True, header_style="bold green")
    t_skani.add_column("skani", style="cyan", no_wrap=True)
    t_skani.add_column("Description", style="white")
    t_skani.add_row("--sparse", "Output row-by-row (sparse) comparisons.")
    t_skani.add_row("--ci", "Output [5%,95%] ANI confidence intervals (bootstrap).")
    t_skani.add_row("--detailed", "Print more stats (e.g., contig N50).")
    t_skani.add_row("--diagonal", "Output diagonal (self-self).")
    t_skani.add_row("--distance", "Output 100 - ANI (distance matrix).")
    t_skani.add_row("--min-af", "Keep pairs with aligned fraction > value (default: 15).")
    t_skani.add_row("--full-matrix", "Output full (not lower-triangular) matrix.")
    t_skani.add_row("-marker_c, --marker_c", "Marker k-mer compression factor (lower for small genomes).")
    t_skani.add_row("-comf, --compression_factor", "k-mer subsampling rate (default: 125).")
    t_skani.add_row("-sim, --similarity_threshold", "Screen out pairs below similarity (default: 80).")

    t_leiden = Table(show_header=True, header_style="bold green")
    t_leiden.add_column("pyleiden", style="cyan", no_wrap=True)
    t_leiden.add_column("Description", style="white")
    t_leiden.add_row("-of, --objective_function", "Objective: CPM or modularity (default: CPM).")
    t_leiden.add_row("-r, --resolution", "Higher → more/smaller clusters (default: 1.0).")
    t_leiden.add_row("-beta, --beta", "Affects randomness during refinement (default: 0.01).")
    t_leiden.add_row("-n, --number_iterations", "Iterations (negative → until stable; default: 2).")
    t_leiden.add_row("-seed, --seed", "Random seed (default: 2024).")
    t_leiden.add_row("-q, --quiet", "Suppress stdout.")
    console.print(t_skani)
    console.print(t_leiden)
    console.print(Rule(style="dim"))


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def main():
    parser = CustomArgumentParser(description="vOTU generation via de-replication (ANI) and clustering (AF).")
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    parser.add_argument("-t", "--thread", type=int, default=40)
    parser.add_argument("--tool", choices=['cd-hit_mmseq', 'skani_pyleiden'], required=True)
    parser.add_argument("-ani", "--drep_ani", type=float, default=0.95)
    parser.add_argument("-cls", "--clustering", type=float, default=0.85)
    parser.add_argument("--parallel", type=int, default=1)
    parser.add_argument("-cf", "--config", type=str, help="Path to config.yml for env/tools/dbs")

    parser.add_argument("-m", "--memory", default=8000, type=int)
    parser.add_argument("-g", "--global_identity", type=int, choices=[0, 1], default=1)
    parser.add_argument("-b", "--band_width", type=int, default=20)
    parser.add_argument("-g_mode", "--accurate_mode", type=int, choices=[0, 1], default=0)
    parser.add_argument("-sc", "--sort_by_size", type=int, choices=[0, 1], default=0)
    parser.add_argument("-sf", "--sort_fasta_by_size", type=int, choices=[0, 1], default=0)

    parser.add_argument("--seed-sub-mat", type=str)
    parser.add_argument("-s", "--sensitivity", type=float, default=4.0)
    parser.add_argument("-k", "--kmer_length", type=int, default=0)
    parser.add_argument("--target-search-mode", type=int, default=0)
    parser.add_argument("--max-seqs", type=int, default=20)
    parser.add_argument("--cov-mode", type=int, default=0)
    parser.add_argument("--alignment-mode", type=int, default=3)
    parser.add_argument("--cluster-mode", type=int, default=0)

    parser.add_argument("--sparse", action='store_true')
    parser.add_argument("--ci", action='store_true')
    parser.add_argument("--detailed", action='store_true')
    parser.add_argument("--diagonal", action='store_true')
    parser.add_argument("--distance", action='store_true')
    parser.add_argument("--min-af", type=float, default=15)
    parser.add_argument("--full-matrix", action='store_true')
    parser.add_argument("-marker_c", "--marker_c", type=int)
    parser.add_argument("-comf", "--compression_factor", type=int, default=125)
    parser.add_argument("-sim", "--similarity_threshold", type=float, default=80)

    parser.add_argument("-of", "--objective_function", choices=['CPM', 'modularity'], default='CPM')
    parser.add_argument("-r", "--resolution", type=float, default=1.0)
    parser.add_argument("-beta", "--beta", type=float, default=0.01)
    parser.add_argument("-n", "--number_iterations", type=int, default=2)
    parser.add_argument("-seed", "--seed", type=int, default=2024)
    parser.add_argument("-q", "--quiet", action='store_true')

    args = parser.parse_args()

    if not args.config:
        sys.exit("ERROR: Please provide -cf/--config pointing to config.yml that contains VMP_env.")
    try:
        vmp_env_prefix = _load_vmp_env_prefix(args.config)
    except Exception as e:
        sys.exit(f"ERROR: Failed to read VMP_env from config: {e}")
    runner = EnvRunner(vmp_env_prefix)

    fasta_list = Getfasta_list(args.input_dir)

    if args.tool == 'cd-hit_mmseq':
        do_cdhit_ANI(
            fasta_list, args.input_dir, args.out_dir, args.thread, args.drep_ani,
            args.memory, args.global_identity, args.band_width, args.accurate_mode,
            args.sort_by_size, args.sort_fasta_by_size, args.parallel, runner
        )
        do_mmseqs(
            fasta_list, args.out_dir, args.clustering, args.seed-sub-mat if hasattr(args, "seed-sub-mat") else None, args.sensitivity,
            args.kmer_length, args.target_search_mode, args.max_seqs, args.cov_mode,
            args.alignment_mode, args.cluster_mode, args.parallel, runner
        )
        clean_and_format_fasta(fasta_list, args.out_dir, args.parallel, runner)

    elif args.tool == 'skani_pyleiden':
        pyleiden(
            fasta_list, args.input_dir, args.out_dir, args.thread, args.drep_ani, args.clustering,
            args.sparse, args.ci, args.detailed, args.diagonal, args.distance, args.min_af,
            args.full_matrix, args.marker_c, args.compression_factor, args.similarity_threshold,
            args.resolution, args.objective_function, args.beta, args.number_iterations,
            args.seed, args.quiet, args.parallel, runner
        )


if __name__ == "__main__":
    main()
