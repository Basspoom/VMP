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
    path = input_dir
    fasta_list = []
    for i in os.listdir(path):
        if os.path.isdir(i):
            continue
        else:
            match = re.search(r'(.*?)\.\w+\.fasta', i)
            if match:
                fasta_list.append(match.group(1))
            else:
                continue
    print('Your sample list:', fasta_list)
    return fasta_list



def binning(input_dir, out_dir, fastq_dir, fasta_list, thread, completion, contamination, metabat2, maxbin2, concoct, min_contig_length, run_checkm, single_end, interleaved, 
            skip_refinement, skip_checkm, skip_consolidation, keep_ambiguous, remove_ambiguous, quick, strict_cutoff, permissive_cutoff, skip_checkm_reassembly, parallel, mdmcleaner):
    
    os.makedirs(os.path.join(out_dir, "bin"), exist_ok=True)

    print("Process bin of non-viral contigs...(binning, bin_refinement, reassemble_bins, and classify_bins)")

    for item in fasta_list:
        command_binning = f"metawrap binning {fastq_dir}/{item}/*_1.fastq {fastq_dir}/{item}/*_2.fastq -o {out_dir}/bin/{item}_bin -t {thread} -a {input_dir}/{item}.*.fasta"

        if metabat2:
            command_binning += " --metabat2"
        if maxbin2:
            command_binning += " --maxbin2"
        if concoct:
            command_binning += " --concoct"
        if min_contig_length:
            command_binning += f" -l {min_contig_length}"
        if run_checkm:
            command_binning += " --run-checkm"
        if single_end:
            command_binning += " --single-end"
        if interleaved:
            command_binning += " --interleaved"

        os.system(command_binning)

        command_refinement = f"metawrap bin_refinement -o {out_dir}/bin/{item}_refbin -t {thread} -A {out_dir}/bin/{item}_bin/metabat2_bins/ -B {out_dir}/bin/{item}_bin/maxbin2_bins/ -C {out_dir}/bin/{item}_bin/concoct_bins/ -c {completion} -x {contamination}"

        if skip_refinement:
            command_refinement += " --skip-refinement"
        if skip_checkm:
            command_refinement += " --skip-checkm"
        if skip_consolidation:
            command_refinement += " --skip-consolidation"
        if keep_ambiguous:
            command_refinement += " --keep-ambiguous"
        if remove_ambiguous:
            command_refinement += " --remove-ambiguous"
        if quick:
            command_refinement += " --quick"

        os.system(command_refinement)

        command_reassembly = f"metawrap reassemble_bins -o {out_dir}/bin/{item}_binReass -1 {fastq_dir}/{item}/*_1.fastq -2 {fastq_dir}/{item}/*_2.fastq -t {thread} -m 800 -c {completion} -x {contamination} -b {out_dir}/bin/{item}_refbin/metawrap_{completion}_{contamination}_bins"

        if strict_cutoff:
            command_reassembly += " --strict-cut-off"
        if permissive_cutoff:
            command_reassembly += " --permissive-cut-off"
        if skip_checkm_reassembly:
            command_reassembly += " --skip-checkm"
        if parallel:
            command_reassembly += " --parallel"
        if mdmcleaner:
            command_reassembly += " --mdmcleaner"

        os.system(command_reassembly)

        command_classification = f"metawrap classify_bins -b {out_dir}/bin/{item}_binReass/reassembled_bins -o {out_dir}/bin/{item}_binClass -t {thread}"
        os.system(command_classification)



def drep(fasta_list, out_dir, thread, ani, completion_drep, contamination_drep, cov_thresh, length, checkm_method, clusterAlg, skip_plots):
    os.makedirs(os.path.join(out_dir, "bin/dRep_output"), exist_ok=True)

    print("Using dRep to remove redundancy...")

    for item in fasta_list:
        command_drep = f"dRep dereplicate {out_dir}/bin/dRep_output/{item} -g {out_dir}/bin/{item}_refbin/metawrap_{completion_drep}_{contamination_drep}_bins/*.fa -p {thread} "

        if ani:
            command_drep += f" -sa {ani}"
        if completion_drep:
            command_drep += f" -comp {completion_drep}"       
        if contamination_drep:
            command_drep += f" -con {contamination_drep}"
        if cov_thresh:
            command_drep += f" -nc {cov_thresh}"       
        if length:
            command_drep += f" -l {length}"
        if checkm_method:
            command_drep += f" --checkM_method {checkm_method}"
        if clusterAlg:
            command_drep += f" --clusterAlg {clusterAlg}"
        if skip_plots:
            command_drep += " --skip_plots"

        os.system(command_drep)

    print("Done remove redundancy!")




def custom_help():
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "This script is designed to use MetaWRAP for binning, dRep for de-duplication, and GTDB for annotation to generate high-quality Metagenome-Assembled Genomes (MAGs) from non-viral contigs." + RESET)

    print("\n" + PURPLE + "Examples:" + RESET)
    print("  meta_wrapper -i /data/input -o /data/output -t 50 --metabat2 --ani 0.95")

    print("\n" + GREEN + "The -i parameter specifies the directory containing input metagenomic non-viral contigs (mNVC) files." + RESET)
    print("  " + BLUE + "/path/to/input/sample1.mNVC.fasta" + RESET)
    print("  " + BLUE + "/path/to/input/sample2.mNVC.fasta" + RESET)

    print("\n" + GREEN + "The -o parameter specifies the output directory for the results." + RESET)

    print("\n" + GREEN + "The -t parameter specifies the number of threads to use for processing (default: 50)." + RESET)

    print("\n" + GREEN + "The -c parameter sets the minimum bin completion percentage (default: 70)." + RESET)

    print("\n" + GREEN + "The -x parameter sets the maximum bin contamination percentage (default: 10)." + RESET)

    print("\n" + PURPLE + "Detailed parameters")
    print("" + GREEN + "=" * 50 + " " + RESET + YELLOW + "Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)

    print("\n" + YELLOW + "Global parameters:" + RESET)
    print(f"  {'-i, --input_dir':<40} Directory containing the metagenomic non-viral contigs (mNVC) files.")
    print(f"  {'-o, --out_dir':<40} Output directory for high-quality non-viral genomes (MAGs).")
    print(f"  {'-t, --thread':<40} Number of parallel threads for processing (default: 50).")
    print(f"  {'-c, --completion':<40} Minimum bin completion percentage (default: 70).")
    print(f"  {'-x, --contamination':<40} Maximum bin contamination percentage (default: 10).")

    print("\n" + GREEN + "=" * 48 + " " + RESET + YELLOW + "Part1: Binning Options" + RESET + " " + GREEN + "=" * 47 + RESET)

    print("\n" + YELLOW + "MetaWRAP Binning Options:" + RESET)
    print(f"  {'--metabat2':<40} Use metaBAT2 for binning.")
    print(f"  {'--maxbin2':<40} Use MaxBin2 for binning.")
    print(f"  {'--concoct':<40} Use CONCOCT for binning.")
    print(f"  {'--min_contig_length':<40} Minimum contig length to bin.")
    print(f"  {'--run_checkm':<40} Run CheckM immediately on the bin results.")
    print(f"  {'--single_end':<40} Use single-end reads.")
    print(f"  {'--interleaved':<40} Input read files contain interleaved paired-end reads.")

    print("\n" + YELLOW + "Refinement Options:" + RESET)
    print(f"  {'--skip_refinement':<40} Don't use binning_refiner.")
    print(f"  {'--skip_checkm':<40} Don't run CheckM to assess bins.")
    print(f"  {'--skip_consolidation':<40} Skip consolidation of bins.")
    print(f"  {'--keep_ambiguous':<40} Keep ambiguous contigs in all bins.")
    print(f"  {'--remove_ambiguous':<40} Remove ambiguous contigs from all bins.")
    print(f"  {'--quick':<40} Quick mode for CheckM.")

    print("\n" + YELLOW + "Reassembly & Classification Options:" + RESET)
    print(f"  {'--strict_cutoff':<40} Use strict cutoff for reassembly.")
    print(f"  {'--permissive_cutoff':<40} Use permissive cutoff for reassembly.")
    print(f"  {'--skip_checkm_reassembly':<40} Don't run CheckM during reassembly.")
    print(f"  {'--parallel':<40} Run SPAdes reassembly in parallel.")
    print(f"  {'--mdmcleaner':<40} Use MDMcleaner results.")


    print("\n" + GREEN + "=" * 49 + " " + RESET + YELLOW + "Part2: dRep Options" + RESET + " " + GREEN + "=" * 49 + RESET)

    print("\n" + YELLOW + "dRep Parameters:" + RESET)
    print(f"  {'-ani, --ani':<40} Average nucleotide identity threshold for secondary clustering (default: 0.95).")
    print(f"  {'-ct, --cov_thresh':<40} Minimum level of overlap between genomes (default: 0.1).")
    print(f"  {'-comp, --completion_drep':<40} Minimum genome completeness (default: None).")
    print(f"  {'-con, --contamination_drep':<40} Maximum genome contamination_drep (default: None).")
    print(f"  {'-l, --length':<40} Minimum genome length (default: None).")
    print(f"  {'--checkm_method':<40} CheckM method to use (default: 'lineage_wf').")
    print(f"  {'--clusterAlg':<40} Clustering algorithm to use (default: 'average').")
    print(f"  {'--skip_plots':<40} If set, skip generating plots (default: False).")


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def main():
    parser = CustomArgumentParser(description="Use MetaWRAP for binning, dRep for de-duplication, and GTDB for annotation to obtain high-quality MAGs from non-viral contigs.")
    
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="The metagenomic non-viral contigs (mNVC, fasta) file directory, such as /xx/mNVC (mNVC/sample1.mNVC.fasta, mNVC/sample2.mNVC.fasta ...)")
    parser.add_argument("-o", "--out_dir", type=str, required=True, help="The location for output MAG (high quality non-viral genomes)")
    parser.add_argument("-t", "--thread", type=int, default=50, help="Number of parallel threads (default: 50)")
    parser.add_argument("-c", "--completion", type=int, default=70, help="Minimum bin completion percentage (default: 70)")
    parser.add_argument("-x", "--contamination", type=int, default=10, help="Maximum bin contamination percentage (default: 10)")

    #### metawrap 
    # binning
    parser.add_argument("--metabat2", action='store_true', help="Use metaBAT2 for binning.")
    parser.add_argument("--maxbin2", action='store_true', help="Use MaxBin2 for binning.")
    parser.add_argument("--concoct", action='store_true', help="Use CONCOCT for binning.")
    parser.add_argument("--min_contig_length", type=int, help="Minimum contig length to bin.")
    parser.add_argument("--run_checkm", action='store_true', help="Run CheckM immediately on the bin results.")
    parser.add_argument("--single_end", action='store_true', help="Use single-end reads.")
    parser.add_argument("--interleaved", action='store_true', help="Input read files contain interleaved paired-end reads.")

    # refinement
    parser.add_argument("--skip_refinement", action='store_true', help="Don't use binning_refiner.")
    parser.add_argument("--skip_checkm", action='store_true', help="Don't run CheckM to assess bins.")
    parser.add_argument("--skip_consolidation", action='store_true', help="Skip consolidation of bins.")
    parser.add_argument("--keep_ambiguous", action='store_true', help="Keep ambiguous contigs in all bins.")
    parser.add_argument("--remove_ambiguous", action='store_true', help="Remove ambiguous contigs from all bins.")
    parser.add_argument("--quick", action='store_true', help="Quick mode for CheckM.")

    # reassembly & classification (Na)
    parser.add_argument("--strict_cutoff", action='store_true', help="Use strict cutoff for reassembly.")
    parser.add_argument("--permissive_cutoff", action='store_true', help="Use permissive cutoff for reassembly.")
    parser.add_argument("--skip_checkm_reassembly", action='store_true', help="Don't run CheckM during reassembly.")
    parser.add_argument("--parallel", action='store_true', help="Run SPAdes reassembly in parallel.")
    parser.add_argument("--mdmcleaner", action='store_true', help="Use MDMcleaner results.")


    #### dRep
    parser.add_argument("-ani", "--ani", type=float, default=0.95, help="Average nucleotide identity threshold for secondary clustering (default: 0.95).")
    parser.add_argument("-ct", "--cov_thresh", type=float, default=0.1, help="Minimum level of overlap between genomes when doing secondary comparisons (default: 0.1).")
    parser.add_argument("-comp", "--completion_drep", type=float, default=None, help="Minimum genome completeness (default: None).")
    parser.add_argument("-con", "--contamination_drep", type=float, default=None, help="Maximum genome contamination (default: None).")
    parser.add_argument("-l", "--length", type=int, default=None, help="Minimum genome length (default: None).")
    parser.add_argument("--checkm_method", type=str, default='lineage_wf', help="CheckM method to use ('lineage_wf' or 'taxonomy_wf', default: 'lineage_wf').")
    parser.add_argument("--clusterAlg", type=str, default='average', help="Clustering algorithm to use (default: 'average').")
    parser.add_argument("--skip_plots", action='store_true',  help="If set, skip generating plots (default: False).")
    
    args = parser.parse_args()

    

    fasta_list = Getfasta_list(args.input_dir)

    binning(args.input_dir, args.out_dir, args.fastq_dir, fasta_list, args.thread, args.completion, args.contamination, args.metabat2, args.maxbin2, args.concoct, args.min_contig_length, args.run_checkm, args.single_end, args.interleaved, 
            args.skip_refinement, args.skip_checkm, args.skip_consolidation, args.keep_ambiguous, args.remove_ambiguous, args.quick, args.strict_cutoff, args.permissive_cutoff, args.skip_checkm_reassembly, args.parallel, args.mdmcleaner)

    drep(fasta_list, args.out_dir, args.thread, args.ani, args.completion, args.completion_drep, args.cov_thresh, args.length, args.checkm_method, args.clusterAlg, args.skip_plots)





if __name__ == "__main__":
    main()





 


        