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
import subprocess
import re
import pandas as pd
from tqdm import tqdm

def Getfasta_list(input_dir):
    path = input_dir
    fasta_list=[]
    for i in os.listdir(path):
        if os.path.isdir(i):
            continue
        else:
            if 'contigs.fa' in i:
                fasta_list.append(i[0:i.index('.')])
            else:
                continue
    print('Your sample list:', fasta_list)
    return fasta_list


    
def run_seqkit(fasta_list, input_dir, out_dir, min_length, max_length, thread):
    print(f"Filering contig length during {min_length} ~ {max_length} kb ...")

    if not os.path.exists(f"{out_dir}/filtered_contigs"):
        os.makedirs(f"{out_dir}/filtered_contigs") 

    for item in fasta_list:
        os.system(f"seqkit seq -m {min_length} -M {max_length} -w 0 -j {thread} {input_dir}/{item}.contigs.fa > {out_dir}/filtered_contigs/{item}.contigs.fasta")


    print(f"Done filer!")


def run_deepvirfinder_single(fasta_list, out_dir, thread, parallel):
    print("Begin DeepVirFinder...")
    
    os.makedirs(f"{out_dir}/dvf", exist_ok=True)
    os.makedirs(f"{out_dir}/dvf/dvf_all", exist_ok=True)

    def process_item(item):
        os.system(f"dvf.py -i {out_dir}/filtered_contigs/{item}.contigs.fasta -o {out_dir}/dvf/{item} -c {thread}")
        os.system(f"cp {out_dir}/dvf/{item}/{item}.*.txt {out_dir}/dvf/dvf_all")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done DeepVirFinder!")


def run_virsorter2_single(fasta_list, out_dir, thread, parallel):
    print("Begin VirSorter2...")
    
    os.makedirs(f"{out_dir}/vs2", exist_ok=True)
    os.makedirs(f"{out_dir}/vs2/vs2_all", exist_ok=True)

    def process_item(item):
        os.system(f"virsorter run -w {out_dir}/vs2/{item}.out -i {out_dir}/filtered_contigs/{item}.contigs.fasta -j {thread} all --keep-original-seq")
        os.system(f"cp {out_dir}/vs2/{item}.out/final-viral-score.tsv {out_dir}/vs2/vs2_all/{item}-final-viral-score.tsv")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done VirSorter2!")


def run_vibrant_single(fasta_list, out_dir, vibrant_db, thread, parallel):
    print("Begin VIBRANT...")
    
    os.makedirs(f"{out_dir}/VIBRANT", exist_ok=True)

    def process_item(item):
        os.system(f"VIBRANT_run.py -t {thread} -d {vibrant_db} -f nucl -virome -i {out_dir}/filtered_contigs/{item}.contigs.fasta -folder {out_dir}/VIBRANT/{item}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done VIBRANT!")


#export CHECKVDB=/home/zly/databases/checkv-db-v1.5
def do_checkV_single(fasta_list, out_dir, checkv_db, thread, parallel):
    print("Begin CheckV...")
    
    os.makedirs(f"{out_dir}/CheckV", exist_ok=True)

    def process_item(item):
        os.system(f"checkv end_to_end -d {checkv_db} -t {thread} {out_dir}/filtered_contigs/{item}.contigs.fasta {out_dir}/CheckV/{item}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done CheckV!")  


def do_kaiju_single(fasta_list, out_dir, kaiju_nodes, kaiju_fmi, kaiju_names, thread, parallel):
    print("Begin kaiju...")
    
    os.makedirs(f"{out_dir}/kaiju", exist_ok=True)

    def process_item(item):
        os.system(f"kaiju -t {kaiju_nodes} -f {kaiju_fmi} -i {out_dir}/filtered_contigs/{item}.contigs.fasta -o {out_dir}/kaiju/{item}.kaiju.tax.out -z {thread} -v")
        os.system(f"kaiju-addTaxonNames -i {out_dir}/kaiju/{item}.kaiju.tax.out -t {kaiju_nodes} -n {kaiju_names} -o {out_dir}/kaiju/{item}.kaiju.names.out -r superkingdom")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done kaiju!")


def run_genomad(fasta_list, genomad_env, genomad_db, out_dir, thread, parallel):
    print("Begin geNomad...")
    
    os.makedirs(f"{out_dir}/genomad", exist_ok=True)

    def process_item(item):
        os.system(f"conda run --prefix {genomad_env} genomad end-to-end {out_dir}/filtered_contigs/{item}.contigs.fasta  {out_dir}/genomad  {genomad_db}  -t {thread}")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        executor.map(process_item, fasta_list)
    
    print("Done geNomad!")



def deal_results(fasta_list, out_dir):

    dirs = [f"{out_dir}/results", f"{out_dir}/results/dvf", f"{out_dir}/results/vs2",
            f"{out_dir}/results/vibr", f"{out_dir}/results/checkV", f"{out_dir}/results/kaiju", f"{out_dir}/results/genomad"]

    for dir_path in dirs:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    print("Integrating key output files...")

    for item in fasta_list:
        os.system(f"cp {out_dir}/dvf/dvf_all/*  {out_dir}/results/dvf")
        os.system(f"cp {out_dir}/vs2/vs2_all/*  {out_dir}/results/vs2")
        os.system(f"cp {out_dir}/VIBRANT/{item}/VIBRANT_{item}.contigs/VIBRANT_results_{item}.contigs/VIBRANT_genome_quality_{item}.contigs.tsv  {out_dir}/results/vibr")
        os.system(f"cp {out_dir}/CheckV/{item}/quality_summary.tsv  {out_dir}/results/checkV/{item}.quality_summary.tsv")
        os.system(f"cp {out_dir}/kaiju/{item}.kaiju.names.out  {out_dir}/results/kaiju")
        
        # with open(f"{out_dir}/genomad/{item}.contigs_marker_classification/{item}.contigs_features.tsv", "r") as file:
        #     content = file.read().replace(" \t", "\t")
        # with open(f"{out_dir}/results/genomad/{item}.contigs_features.tsv", "w") as file:
        #     file.write(content)

        os.system(f"cp {out_dir}/genomad/{item}.contigs_marker_classification/{item}.contigs_features.tsv  {out_dir}/results/genomad")
        os.system(f"cp {out_dir}/genomad/{item}.contigs_aggregated_classification/{item}.contigs_aggregated_classification.tsv  {out_dir}/results/genomad")

        # os.system(f"python3 deal_kaiju_out.py")
        

def get_scores(out_dir):

    print("Perform three rounds of scoring: single tool, tunning addition and tunning remove.")

    print("Step1: Single tool")
    os.system(f"VPAC_part2_dvf.py  {out_dir}/results/dvf") 
    os.system(f"VPAC_part2_vs2.py  {out_dir}/results/vs2")  
    os.system(f"VPAC_part2_vibr.py  {out_dir}/results/vibr")  
    os.system(f"VPAC_part2_genomad.py  {out_dir}/results/genomad") 

    print("Step2: Tunning addition")
    os.system(f"VPAC_part2_rmv.py  -c {out_dir}/results/checkV  -v {out_dir}/results/vs2 -g {out_dir}/results/genomad")  

    print("Step3: Tunning remove")
    os.system(f"VPAC_part2_add.py -k {out_dir}/results/kaiju  -c {out_dir}/results/checkV  -v {out_dir}/results/vs2 -g {out_dir}/results/genomad")  



def score_summary(out_dir, input_dir, fasta_list):
    dirs = [f"{out_dir}/score", f"{out_dir}/summary",
            f"{out_dir}/score/part1", f"{out_dir}/score/part2", f"{out_dir}/score/part3"]

    for dir_path in dirs:
        os.makedirs(dir_path, exist_ok=True)

    print("Integrating the results of three rounds of scoring...")
    os.system(f"cp {out_dir}/results/dvf/*.tsv  {out_dir}/score/part1")
    os.system(f"cp {out_dir}/results/vs2/*_vs2.tsv  {out_dir}/score/part1")
    os.system(f"cp {out_dir}/results/vibr/*_vibrant.tsv  {out_dir}/score/part1")
    os.system(f"cp {out_dir}/results/genomad/*_genomad.tsv  {out_dir}/score/part1")

    os.system(f"cp {out_dir}/results/kaiju/*scores.tsv  {out_dir}/score/part2")
    os.system(f"cp {out_dir}/results/checkV/*_checkV_add.tsv  {out_dir}/score/part2")
    os.system(f"cp {out_dir}/results/vs2/*_vs2_add.tsv  {out_dir}/score/part2")
    os.system(f"cp {out_dir}/results/genomad/*_genomad_add.tsv  {out_dir}/score/part2")

    os.system(f"cp {out_dir}/results/checkV/*_checkV_rmv.tsv  {out_dir}/score/part3")
    os.system(f"cp {out_dir}/results/vs2/*_vs2_rmv.tsv  {out_dir}/score/part3")
    os.system(f"cp {out_dir}/results/genomad/*_genomad_rmv.tsv  {out_dir}/score/part3")

    contigs_folder = input_dir
    score_folders = [f"{out_dir}/score/part1", f"{out_dir}/score/part2", f"{out_dir}/score/part3"]

    score_files_info = {
        "part1": ["*.contigs.fasta_gt1bp_dvfpred_dvf.tsv", "*-final-viral-score_vs2.tsv", "VIBRANT_genome_quality_*.contigs_vibrant.tsv", "*.contigs_aggregated_classification_genomad.tsv"],
        "part2": ["*.kaiju.names_kaiju_scores.tsv", "*.quality_summary_checkV_add.tsv", "*-final-viral-score_vs2_add.tsv", "*.contigs_features_genomad_add.tsv"],
        "part3": ["*.quality_summary_checkV_rmv.tsv", "*-final-viral-score_vs2_rmv.tsv", "*.contigs_features_genomad_rmv.tsv"]
    }

    def custom_sort(seq_id):
        match = re.match(r"k\d+_(\d+)", seq_id)
        return int(match.group(1)) if match else 0

    for sample_name in tqdm(fasta_list, desc="Processing samples"):
        contigs_file = os.path.join(contigs_folder, f"{sample_name}.contigs.fa")
        sequence_ids = {line.strip().split()[0][1:] for line in open(contigs_file) if line.startswith('>')}
        sequence_ids = sorted(sequence_ids, key=custom_sort)

        summary_df = pd.DataFrame({"Sequence ID": sequence_ids})

        for step, score_folder in enumerate(tqdm(score_folders, desc=f"Processing steps for {sample_name}", leave=False), start=1):
            score_files_pattern = score_files_info.get(f"part{step}")
            if score_files_pattern is None:
                # print(f"未找到与步骤 part{step} 相关的得分文件信息")
                print("quit")
                continue

            for pattern in score_files_pattern:
                pattern_with_sample = pattern.replace("*", sample_name)
                score_files = [file_name for file_name in os.listdir(score_folder) if file_name.endswith(pattern_with_sample)]

                for score_file in score_files:
                    score_file_path = os.path.join(score_folder, score_file)
                    df = pd.read_csv(score_file_path, sep='\t')
                    df[df.columns[1:]] = df[df.columns[1:]].fillna(0).astype(int)
                    df.columns = [f"{step}_{score_file.split('.')[0]}_score{col}" if idx > 0 else col for idx, col in enumerate(df.columns)]
                    summary_df = summary_df.merge(df, how='left', left_on='Sequence ID', right_on=df.columns[0]).fillna(0)

        summary_df.iloc[:, 1:] = summary_df.iloc[:, 1:].astype(int)

        output_file = f"{out_dir}/summary/{sample_name}_summary.tsv"
        summary_df.to_csv(output_file, index=False, sep='\t')
        # print(f"已生成汇总表: {output_file}")



def get_viralseq(input_dir, out_dir, fasta_list, model):

    dirs = [f"{out_dir}/viral_score", f"{out_dir}/labels", f"{out_dir}/mVC", f"{out_dir}/mNVC"]

    for dir_path in dirs:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)


    if model == "best_f1":
        # cmd = "which model_best_f1.joblib"
        cmd = "/opt/user/basspoom/GUI/model/VPAC/random_forest_model.joblib"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        path = result.stdout.decode().strip()
        predict_model = path

    elif model == "best_precision":
        cmd = "which model_best_precision.joblib"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        path = result.stdout.decode().strip()
        predict_model = path


    elif model == "best_recall":
        cmd = "which model_best_recall.joblib"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        path = result.stdout.decode().strip()
        predict_model = path

    print(f"Predict viral and non-viral contigs by {model} model...")



    for item in fasta_list:
        # os.system(f"predict.py  {out_dir}/summary/{item}_summary.tsv  {predict_model}  {out_dir}/viral_score/{item}_score.txt")
        os.system(f"predict.py  {out_dir}/summary/{item}_summary.tsv  /opt/user/basspoom/GUI/model/VPAC/random_forest_model.joblib  {out_dir}/viral_score/{item}_score.txt")
        os.system(f"grep '>' {input_dir}/{item}.contigs.fa > {out_dir}/viral_score/temp1.txt")
        os.system(f"awk -F ' ' '{{print $1}}' {out_dir}/viral_score/temp1.txt > {out_dir}/viral_score/temp2.txt")
        os.system(f"sed 's/>//' {out_dir}/viral_score/temp2.txt > {out_dir}/labels/{item}_all_labels.txt")
        os.system(f'awk \'NR>1 && $3=="Virus"{{print $1}}\' {out_dir}/viral_score/{item}_score.txt > {out_dir}/labels/{item}_viral_labels.txt')
        os.system(f'grep -Fv -f {out_dir}/labels/{item}_viral_labels.txt {out_dir}/labels/{item}_all_labels.txt > {out_dir}/labels/{item}_non_viral_labels.txt')
        os.system(f'seqkit grep -w 0 -f {out_dir}/labels/{item}_viral_labels.txt  {input_dir}/{item}.contigs.fa > {out_dir}/mVC/{item}.mVC.fasta')
        os.system(f'seqkit grep -w 0 -f {out_dir}/labels/{item}_non_viral_labels.txt  {input_dir}/{item}.contigs.fa > {out_dir}/mNVC/{item}.mNVC.fasta')

    # print("All processes have been done!")
    print(f"Your viral contigs are in path '{out_dir}/mVC/' ")
    print(f"Your non-viral contigs are in path '{out_dir}/mNVC/'")




def check_viralseq(out_dir):

    import os
    import glob
    import subprocess
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from io import StringIO

    for mode in ['mVC', 'mNVC']:
        contigs_dir = os.path.join(out_dir, mode)
        stat_dir = os.path.join(contigs_dir, 'contigs_stat')
        
        os.makedirs(stat_dir, exist_ok=True)
        
        fa_files = glob.glob(os.path.join(contigs_dir, "*.fasta"))
        if not fa_files:
            print(f"No contig files found in {contigs_dir} for statistics.")
            continue
        
        try:
            stats_output = subprocess.check_output(
                f"seqkit stat -a {' '.join(fa_files)}",
                shell=True, text=True, stderr=subprocess.STDOUT
            )
            stats_df = pd.read_csv(StringIO(stats_output), sep='\t')
            stats_summary_path = os.path.join(stat_dir, 'contig_stats_summary.tsv')
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
            sample_name = os.path.basename(file_path).replace('.contigs.fa', '')
            print(f"[{mode}] Processing length distribution for {sample_name}...")
            try:
                cmd = f"seqkit fx2tab -n -l {file_path} | cut -f 2"
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
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
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        plot_pdf = os.path.join(stat_dir, 'contig_length_distribution.pdf')
        plot_png = os.path.join(stat_dir, 'contig_length_distribution.png')
        fig.savefig(plot_pdf, bbox_inches='tight')
        fig.savefig(plot_png, dpi=500, bbox_inches='tight')
        plt.close()
        print(f"[{mode}] Plots saved to:\n  - {plot_pdf}\n  - {plot_png}")

        print("All processes have been done!")




def custom_help():
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "This script integrates DeepVirFinder, VirSorter2, and VIBRANT to identify viral contigs, then runs kaiju and CheckV to check quality." + RESET)
    print("\n" + PURPLE + "Examples:" + RESET)
    print("  VPAC_predict -i /data/input -o /data/output --min_length 3000 --max_length 50000 --vibrant_db /opt/user/xxx/databases/vibrant/databases/ --checkv_db /opt/user/xxx/databases/checkv-db-v1.5/ --kaiju_nodes /opt/user/xxx/kaiju/nodes.dmp --kaiju_names /opt/user/xxx/kaiju/names.dmp --kaiju_fmi /opt/user/xxx/kaiju/kaiju_db_refseq.fmi --thread 40 --genomad_env /opt/user/xxx/tools/anaconda3/envs/genomad --genomad_db /opt/user/xxx/databases/genomad_db/ --model best_f1")
    print("\n" + GREEN + "The -i parameter specifies the directory containing input FASTA files." + RESET)
    print("  " + BLUE + "/path/to/input/sample1.fasta" + RESET)
    print("  " + BLUE + "/path/to/input/sample2.fasta" + RESET)
    print("\n" + GREEN + "The -o parameter specifies the output directory for the results." + RESET)
    print("\n" + GREEN + "The --min_length parameter specifies the minimum contig length for analysis (default: 3000)." + RESET)
    print("\n" + GREEN + "The --max_length parameter specifies the maximum contig length for analysis (default: 50000)." + RESET)
    print("\n" + GREEN + "The --vibrant_db parameter specifies the location of VIBRANT's database." + RESET)
    print("\n" + GREEN + "The --checkv_db parameter specifies the location of CheckV's database." + RESET)
    print("\n" + GREEN + "The --kaiju_nodes parameter specifies the location of kaiju's nodes.dmp file." + RESET)
    print("\n" + GREEN + "The --kaiju_names parameter specifies the location of kaiju's names.dmp file." + RESET)
    print("\n" + GREEN + "The --kaiju_fmi parameter specifies the location of kaiju's kaiju_db_refseq.fmi file." + RESET)
    print("\n" + GREEN + "The --thread parameter specifies the number of parallel threads (default: 40)." + RESET)
    print("\n" + GREEN + "The --genomad_env parameter specifies the location of geNomad's environment." + RESET)
    print("\n" + GREEN + "The --genomad_db parameter specifies the location of geNomad's database." + RESET)
    print("\n" + GREEN + "The --model parameter allows you to choose the prediction model (default: 'best_f1')." + RESET)
    print("\n" + PURPLE + "Detailed parameters")
    print("" + GREEN + "=" * 50 + " " + RESET + YELLOW + "Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)
    print("\n" + YELLOW + "Global parameters:" + RESET)
    print(f"  {'-i, --input_dir':<40} Directory containing the input fasta files.")
    print(f"  {'-o, --out_dir':<40} Output directory for results.")
    print(f"  {'-lmin, --min_length':<40} Minimum contig length for analysis (default: 3000).")
    print(f"  {'-lmax, --max_length':<40} Maximum contig length for analysis (default: 50000).")
    print(f"  {'--vibrant_db':<40} Location of VIBRANT's database.")
    print(f"  {'--checkv_db':<40} Location of CheckV's database.")
    print(f"  {'--kaiju_nodes':<40} Location of kaiju's nodes.dmp file.")
    print(f"  {'--kaiju_names':<40} Location of kaiju's names.dmp file.")
    print(f"  {'--kaiju_fmi':<40} Location of kaiju's kaiju_db_refseq.fmi file.")
    print(f"  {'-t, --thread':<40} Number of parallel threads (default: 40).")
    print(f"  {'--genomad_env':<40} Location of geNomad's environment.")
    print(f"  {'--genomad_db':<40} Location of geNomad's database.")
    print(f"  {'-m, --model':<40} Choose the prediction model (default: 'best_f1').")
    print(f"  {'--parallel':<40} Number of samples to process in parallel (default: 1).")


class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="This script integrates DeepVirFinder, VirSorter2, and VIBRANT to identify viral contigs, then runs kaiju and CheckV to check quality.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="The assembled contig (fasta) file directory, such as /xx/raw_contigs (raw_contigs/sample1.fa, raw_contigs/sample2.fa ...)")
    parser.add_argument("-o", "--out_dir", type=str, required=True, help="The location of viral prediction results and middle files")
    parser.add_argument("-lmin", "--min_length", type=int, default=3000, help="Minimum contig length for analysis (default: 3000)")
    parser.add_argument("-lmax", "--max_length", type=int, default=50000, help="Maximum contig length for analysis (default: 50000)")
    parser.add_argument("-vibrant_db", "--vibrant_db", type=str, required=True, help="The location of VIBRANT's database, for example '/opt/user/xxx/databases/vibrant/databases/' (path to original 'databases' directory that contains .HMM files)")
    parser.add_argument("-checkv_db", "--checkv_db", type=str, required=True, help="The location of CheckV's database, for example '/opt/user/xxx/databases/checkv-db-v1.5/' ")
    parser.add_argument("-kaiju_nodes", "--kaiju_nodes", type=str, required=True, help="The location of kaiju's nodes.dmp file")
    parser.add_argument("-kaiju_names", "--kaiju_names", type=str, required=True, help="The location of kaiju's names.dmp file")
    parser.add_argument("-kaiju_fmi", "--kaiju_fmi", type=str, required=True, help="The location of kaiju's kaiju_db_refseq.fmi file")
    parser.add_argument("-t", "--thread", type=int, default=40, help="Number of parallel threads (default: 40)")
    parser.add_argument("-genomad_env", "--genomad_env", type=str, required=True, help="The location of geNomad's env, for example '/opt/user/xxx/tools/anaconda3/envs/genomad' ")
    parser.add_argument("-genomad_db", "--genomad_db", type=str, required=True, help="The location of geNomad's database, for example '/opt/user/xxx/databases/genomad_db/' ")
    parser.add_argument("-m", "--model", type=str, choices=["best_f1", "best_precision", "best_recall"], default="best_f1", help='Choose your prediction model: "best_precision" to minimize false positives; "best_recall" to recall as many true virus sequences as possible; "best_f1" a model with a balanced precision and recall')
    
    parser.add_argument("--parallel", type=int, default=1, help="Number of samples to process in parallel (default: 1)")
    
    args = parser.parse_args()


    fasta_list = Getfasta_list(args.input_dir)
    run_seqkit(fasta_list, args.input_dir, args.out_dir, args.min_length, args.max_length, args.thread)

    run_deepvirfinder_single(fasta_list, args.out_dir, args.thread, args.parallel)
    run_vibrant_single(fasta_list, args.out_dir, args.vibrant_db, args.thread, args.parallel)
    do_checkV_single(fasta_list, args.out_dir, args.checkv_db, args.thread, args.parallel)
    do_kaiju_single(fasta_list, args.out_dir, args.kaiju_nodes, args.kaiju_fmi, args.kaiju_names, args.thread, args.parallel)
    run_genomad(fasta_list, args.genomad_env, args.genomad_db, args.out_dir, args.thread, args.parallel)
    run_virsorter2_single(fasta_list, args.out_dir, args.thread, args.parallel)

    deal_results(fasta_list, args.out_dir)
    get_scores(args.out_dir)
    score_summary(args.out_dir, args.input_dir, fasta_list)
    get_viralseq(args.input_dir, args.out_dir, fasta_list, args.model)
    check_viralseq(args.out_dir)




if __name__ == "__main__":
    main()