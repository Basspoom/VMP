import os
import glob
import datetime
import types
import csv
import concurrent.futures
import subprocess
import argparse
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import re
import shutil
import argparse
import subprocess


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


def do_cdhit_ANI(fasta_list, input_dir, out_dir, thread, drep_ani, memory, global_identity, band_width, accurate_mode, sort_by_size, sort_fasta_by_size, parallel):
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

        exit_code = os.system(cdhit_command)
        if exit_code != 0:
            print(f"Error occurred during CD-HIT processing for {item}. Command: {cdhit_command}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_cdhit, fasta_list)

    print("Done removing redundancy!")


def do_mmseqs(fasta_list, out_dir, clustering, seed_sub_mat, sensitivity, kmer_length, target_search_mode, max_seqs, cov_mode, alignment_mode, cluster_mode, parallel):
    os.makedirs(os.path.join(out_dir, "vOTU/cluster"), exist_ok=True)
    print("Using MMseqs2 to cluster viral contigs...")
    print(f"Alignment fraction = {clustering}")

    def process_mmseqs(item):
        item_db_dir = os.path.join(out_dir, "vOTU/cluster/vOTU", item)
        os.makedirs(item_db_dir, exist_ok=True)

        command1 = f"mmseqs createdb {os.path.join(out_dir, 'vOTU/cd-hit_drp', item + '.drp.fasta')} {item_db_dir}/{item}-contigsDB"
        os.system(command1)

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

        os.system(command)

        os.makedirs(os.path.join(out_dir, "vOTU/cluster/vOTU", item), exist_ok=True)

        os.system(f"mmseqs createtsv {item_db_dir}/{item}-contigsDB {out_dir}/vOTU/cluster/vOTU/{item}.cluster {out_dir}/vOTU/cluster/{item}.cluster.tsv")
        os.system(f"mmseqs createseqfiledb {item_db_dir}/{item}-contigsDB  {out_dir}/vOTU/cluster/vOTU/{item}.cluster  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-DB_clu_seq")
        os.system(f"mmseqs result2flat {item_db_dir}/{item}-contigsDB  {item_db_dir}/{item}-contigsDB  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-DB_clu_seq  {out_dir}/vOTU/cluster/vOTU/{item}/{item}-clu_seq.fasta")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_mmseqs, fasta_list)

    print("Done clustering!")


def clean_and_format_fasta(fasta_list, out_dir, parallel):
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

        os.system(f"seqkit seq -m 1 --only-id --upper-case -w 0 {os.path.join(out_dir, f'{item}.v.fasta')} -o {os.path.join(out_dir, f'{item}.vOTU.fasta')}")
        os.system(f"rm {os.path.join(out_dir, f'{item}.v.fasta')}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_clean, fasta_list)





def pyleiden(fasta_list, input_dir, out_dir, thread, drep_ani, clustering, sparse, ci, detailed, diagonal, distance, min_af, full_matrix, 
             marker_c, compression_factor, similarity_threshold, resolution, objective_function, beta, number_iterations, seed, quiet, parallel):

    os.makedirs(os.path.join(out_dir, "vOTU/pyleiden"), exist_ok=True)
    print("Using aksni and pyleiden to remove redundancy and cluster viral contigs ...")
    print(f"Average Nucleotide Identity = {drep_ani}")
    print(f"Alignment fraction = {clustering}")

    def process_pyleiden(item):
        os.makedirs(os.path.join(out_dir, "vOTU/pyleiden", item), exist_ok=True)
        
        os.system(f"seqkit seq --only-id --upper-case -w 0 {os.path.join(input_dir, item + '.*.fasta')} -o {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.viral_contigs.fasta')}")
        
        skani_command = f"skani triangle -i -t {thread} {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.viral_contigs.fasta')} > {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.skani_output.tsv')}"

        if sparse:
            skani_command += " --sparse"
        if ci:
            skani_command += " --ci"
        if detailed:
            skani_command += " --detailed"
        if diagonal:
            skani_command += " --diagonal"
        if distance:
            skani_command += " --distance"
        if min_af is not None:
            skani_command += f" --min-af {min_af}"
        if full_matrix:
            skani_command += " --full-matrix"
        if marker_c is not None:
            skani_command += f" -m {marker_c}"
        if compression_factor is not None:
            skani_command += f" -c {compression_factor}"
        if similarity_threshold is not None:
            skani_command += f" -s {similarity_threshold}"

        exit_code = os.system(skani_command)
        if exit_code != 0:
            print(f"Error occurred during SKANI processing for {item}. Command: {skani_command}")

        awk_command = f"awk -v drep_ani={drep_ani} -v clustering={clustering} 'NR>1 && $3 >= drep_ani && ($4 >= clustering || $5 >= clustering) {{ printf(\"%s\\t%s\\t%.4f\\n\", $6, $7, $3 * ($4 > $5 ? $4 : $5) / 10000) }}' {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.skani_output.tsv')} > {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.skani_edges.tsv')}"
        os.system(awk_command)

        seqkit_command = (
            f"seqkit fx2tab -ni {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.viral_contigs.fasta')} | "
            "awk '{ print $1 \"\\t\" $1 \"\\t1.0000\" }' >> "
            f"{os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.skani_edges.tsv')}"
        )
        os.system(seqkit_command)

        pyleiden_command = f"pyleiden {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.skani_edges.tsv')} {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.clusters.txt')}"

        if resolution:
            pyleiden_command += f" -r {resolution}"
        if objective_function:
            pyleiden_command += f" -o {objective_function}"
        if beta:
            pyleiden_command += f" -b {beta}"
        if number_iterations is not None:
            pyleiden_command += f" -n {number_iterations}"
        if seed:
            pyleiden_command += f" -s {seed}"
        if quiet:
            pyleiden_command += " -q"

        exit_code = os.system(pyleiden_command)
        if exit_code != 0:
            print(f"Error occurred during pyleiden processing for {item}. Command: {pyleiden_command}")


        os.system(f"cut -f1 {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.clusters.txt')} > {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.labels.txt')}")
        os.system(f"seqkit grep -w 0 -f  {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.mVC.labels.txt')} {os.path.join(out_dir, 'vOTU/pyleiden', item, f'{item}.viral_contigs.fasta')} > {os.path.join(out_dir, f'{item}.vOTU.fasta')}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_pyleiden, fasta_list)


def custom_help():
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "This script is for de-duplication (by ANI) and clustering (by AF) to generate high-quality vOTUs from viral contigs. Two methods you can choose: 'cd-hit + MMseqs2' or 'skani + pyleiden'." + RESET)

    print("\n" + PURPLE + "Examples:" + RESET)
    print("  vOTU -i /data/input -o /data/output --tool pyleiden  --drep_ani 0.95 --clustering 0.85 -t 40")

    print("\n" + GREEN + "The -i parameter specifies the directory containing input FASTA files." + RESET)
    print("  " + BLUE + "/path/to/input/sample1.fasta" + RESET)
    print("  " + BLUE + "/path/to/input/sample2.fasta" + RESET)

    print("\n" + GREEN + "The -o parameter specifies the output directory for the results." + RESET)

    print("\n" + GREEN + "The --drep_ani parameter sets the threshold for de-duplication by Average Nucleotide Identity (ANI) (default: 0.95)." + RESET)

    print("\n" + GREEN + "The --clustering parameter defines the Alignment Fraction (AF) for clustering (default: 0.85)." + RESET)

    print("\n" + GREEN + "The -t parameter specifies the number of threads to use for processing (default: 40)." + RESET)

    print("\n" + PURPLE + "Detailed parameters")
    print("" + GREEN + "=" * 50 + " " + RESET + YELLOW + "Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)

    print("\n" + YELLOW + "Global parameters:" + RESET)
    print(f"  {'-i, --input_dir':<40} Directory containing the input fasta files.")
    print(f"  {'-o, --out_dir':<40} Output directory for results.")
    print(f"  {'-ani, --drep_ani':<40} Threshold for de-duplication (default: 0.95).")
    print(f"  {'-cls, --clustering':<40} Alignment fraction for clustering (default: 0.85).")
    print(f"  {'-t, --thread':<40} Number of threads for processing (default: 40).")
    print(f"  {'--tool':<40} Assembly tool to use: 'cd-hit_mmseq' or 'pyleiden'.")


    print("\n" + GREEN + "=" * 46 + " " + RESET + YELLOW + "Choose1: cd-hit & MMseqs2" + RESET + " " + GREEN + "=" * 46 + RESET)

    print("\n" + YELLOW + "CD-HIT Parameters:" + RESET)
    print(f"  {'-m, --memory':<40} Memory limit for CD-HIT (in MB).")
    print(f"  {'-g, --global_identity':<40} Use global sequence identity (0 for local, 1 for global; default: 1).")
    print(f"  {'-b, --band_width':<40} Band width of alignment (default: 20).")
    print(f"  {'-g_mode, --accurate_mode':<40} Use accurate clustering mode (0 for fast, 1 for accurate; default: 0).")
    print(f"  {'-sc, --sort_by_size':<40} Sort clusters by size (0 for no sorting, 1 for sorting; default: 0).")
    print(f"  {'-sf, --sort_fasta_by_size':<40} Sort fasta/fastq by cluster size (0 for no sorting, 1 for sorting; default: 0).")

    print("\n" + YELLOW + "MMseqs2 Parameters:" + RESET)
    print(f"  {'--seed-sub-mat':<40} Substitution matrix file for k-mer generation.")
    print(f"  {'-s, --sensitivity':<40} Sensitivity (1.0 faster; 4.0 fast; 7.5 sensitive; default: 4.0).")
    print(f"  {'-k, --kmer_length':<40} k-mer length (0: automatically set to optimum).")
    print(f"  {'--target-search-mode':<40} Target search mode (0: regular k-mer, 1: similar k-mer; default: 0).")
    print(f"  {'--max-seqs':<40} Maximum results per query sequence allowed to pass the prefilter (default: 20).")
    print(f"  {'--cov-mode':<40} Coverage mode (0: coverage of query and target; default: 0).")
    print(f"  {'--alignment-mode':<40} How to compute the alignment (default: 3).")
    print(f"  {'--cluster-mode':<40} Clustering mode (0: Set-Cover, 1: Connected component; default: 0).")


    print("\n" + GREEN + "=" * 46 + " " + RESET + YELLOW + "Choose2: skani & pyleiden" + RESET + " " + GREEN + "=" * 46 + RESET)

    print("\n" + YELLOW + "skani Parameters:" + RESET)
    print(f"  {'--sparse':<40} Output comparisons in a row-by-row form (sparse matrix).")
    print(f"  {'--ci':<40} Output [5%,95%] ANI confidence intervals using percentile bootstrap.")
    print(f"  {'--detailed':<40} Print additional info including contig N50s and more.")
    print(f"  {'--diagonal':<40} Output the diagonal of the ANI matrix (self-self comparisons).")
    print(f"  {'--distance':<40} Output 100 - ANI instead of ANI, creating a distance matrix.")
    print(f"  {'--min-af':<40} Only output ANI values where one genome has aligned fraction > this value (default: 15).")
    print(f"  {'--full-matrix':<40} Output full matrix instead of lower-triangular matrix.")
    print(f"  {'-marker_c':<40} Marker k-mer compression factor; consider decreasing for small genomes.")
    print(f"  {'-c, --compression_factor':<40} Compression factor (k-mer subsampling rate; default: 125).")
    print(f"  {'-sim, --similarity_threshold':<40} Screen out pairs with approximately < identity using k-mer sketching (default: 80).")

    print("\n" + YELLOW + "Pyleiden Parameters:" + RESET)
    print(f"  {'-of, --objective_function':<40} Use the Constant Potts Model (CPM) or modularity (default: CPM).")
    print(f"  {'-r, --resolution':<40} The resolution parameter to use (default: 1.0).")
    print(f"  {'-beta, --beta':<40} Parameter affecting the randomness in the Leiden algorithm (default: 0.01).")
    print(f"  {'-n, --number_iterations':<40} The number of iterations for the Leiden algorithm (default: 2).")
    print(f"  {'-seed, --seed':<40} Seed for the random number generator (default: 2024).")
    print(f"  {'-q, --quiet':<40} Run quietly (will not print output to stdout; default: False).")


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()



def main():
    parser = CustomArgumentParser(description="Upstream of metagenomics (quality control and contigs assembly).")
    # globalabout:blank#blocked
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="Directory containing the input fasta files.")
    parser.add_argument("-o", "--out_dir", type=str, required=True, help="Output directory for results.")
    parser.add_argument("-t", "--thread", type=int, default=40, help="Number of threads for processing.")

    parser.add_argument("--tool", choices=['cd-hit_mmseq', 'pyleiden'], required=True, help="Assembly tool to use: 'cd-hit_mmseq' or 'pyleiden'.")

    parser.add_argument("-ani", "--drep_ani", type=float, default=0.95, help="Threshold for CD-HIT de-duplication.")
    parser.add_argument("-cls", "--clustering", type=float, default=0.85, help="Alignment fraction for clustering.")


    # cd-hit
    parser.add_argument("-m", "--memory", default=8000, type=int, help="Memory limit for CD-HIT (in MB).")
    parser.add_argument("-g", "--global_identity", type=int, choices=[0, 1], default=1, help="Use global sequence identity (0 for local, 1 for global).")
    parser.add_argument("-b", "--band_width", type=int, default=20, help="Band width of alignment.")
    parser.add_argument("-g_mode", "--accurate_mode", type=int, choices=[0, 1], default=0, help="Use accurate clustering mode (0 for fast, 1 for accurate).")
    parser.add_argument("-sc", "--sort_by_size", type=int, choices=[0, 1], default=0, help="Sort clusters by size (0 for no sorting, 1 for sorting).")
    parser.add_argument("-sf", "--sort_fasta_by_size", type=int, choices=[0, 1], default=0, help="Sort fasta/fastq by cluster size (0 for no sorting, 1 for sorting).")

    # mmseqs2
    parser.add_argument("--seed-sub-mat", type=str, help="Substitution matrix file for k-mer generation.")
    parser.add_argument("-s", "--sensitivity", type=float, default=4.0, help="Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive.")
    parser.add_argument("-k", "--kmer_length", type=int, default=0, help="k-mer length (0: automatically set to optimum).")
    parser.add_argument("--target-search-mode", type=int, default=0, help="Target search mode (0: regular k-mer, 1: similar k-mer).")
    parser.add_argument("--max-seqs", type=int, default=20, help="Maximum results per query sequence allowed to pass the prefilter.")
    parser.add_argument("--cov-mode", type=int, default=0, help="Coverage mode (0: coverage of query and target).")
    parser.add_argument("--alignment-mode", type=int, default=3, help="How to compute the alignment.")
    parser.add_argument("--cluster-mode", type=int, default=0, help="Clustering mode (0: Set-Cover, 1: Connected component).")


    # skani
    parser.add_argument("--sparse", action='store_true', help="Output comparisons in a row-by-row form (sparse matrix).")
    parser.add_argument("--ci", action='store_true', help="Output [5%,95%] ANI confidence intervals using percentile bootstrap.")
    parser.add_argument("--detailed", action='store_true', help="Print additional info including contig N50s and more.")
    parser.add_argument("--diagonal", action='store_true', help="Output the diagonal of the ANI matrix (self-self comparisons).")
    parser.add_argument("--distance", action='store_true', help="Output 100 - ANI instead of ANI, creating a distance matrix.")
    parser.add_argument("--min-af", type=float, default=15, help="Only output ANI values where one genome has aligned fraction > this value.")
    parser.add_argument("--full-matrix", action='store_true', help="Output full matrix instead of lower-triangular matrix.")
    parser.add_argument("-marker_c", "--marker_c", type=int, help="Marker k-mer compression factor; consider decreasing for small genomes.")
    parser.add_argument("-cf", "--compression_factor", type=int, default=125, help="Compression factor (k-mer subsampling rate).")
    parser.add_argument("-sim", "--similarity_threshold", type=float, default=80, help="Screen out pairs with approximately < identity using k-mer sketching.")

    # pyleiden
    parser.add_argument("-of", "--objective_function", choices=['CPM', 'modularity'], default='CPM', help="Whether to use the Constant Potts Model (CPM) or modularity. (default: CPM)")
    parser.add_argument("-r", "--resolution", type=float, default=1.0, help="The resolution parameter to use. Higher resolutions lead to more smaller clusters, while lower resolutions lead to fewer larger clusters. (default: 1.0)")
    parser.add_argument("-beta", "--beta", type=float, default=0.01, help="Parameter affecting the randomness in the Leiden algorithm. This affects only the refinement step of the algorithm. (default: 0.01)")
    parser.add_argument("-n", "--number_iterations", type=int, default=2, help="The number of iterations to iterate the Leiden algorithm. Using a negative number will run until a stable iteration is encountered. (default: 2)")
    parser.add_argument("-seed", "--seed", type=int, default=2024, help="Seed for the random number generator. (default: 2024)")
    parser.add_argument("-q", "--quiet", action='store_true', help="Run quietly (will not print output to stdout). (default: False)")

    parser.add_argument("--parallel", type=int, default=1, help="Number of samples to process in parallel (default: 1)")



    args = parser.parse_args()




    fasta_list = Getfasta_list(args.input_dir)

    if args.tool == 'cd-hit_mmseq':
        do_cdhit_ANI(fasta_list, args.input_dir, args.out_dir, args.thread, args.drep_ani, args.memory, args.global_identity, args.band_width, args.accurate_mode, args.sort_by_size, args.sort_fasta_by_size, args.parallel)
        do_mmseqs(fasta_list, args.out_dir, args.clustering, args.seed_sub_mat, args.sensitivity, args.kmer_length, args.target_search_mode, args.max_seqs, args.cov_mode, args.alignment_mode, args.cluster_mode, args.parallel)
        clean_and_format_fasta(fasta_list, args.out_dir, args.parallel)
    
    elif args.tool == 'pyleiden':
        pyleiden(fasta_list, args.input_dir, args.out_dir, args.thread, args.drep_ani, args.clustering, args.sparse, args.ci, args.detailed, args.diagonal, args.distance, args.min_af, args.full_matrix, args.marker_c, args.compression_factor, args.similarity_threshold, args.resolution, args.objective_function, args.beta, args.number_iterations, args.seed, args.quiet, args.parallel)


if __name__ == "__main__":
    main()
