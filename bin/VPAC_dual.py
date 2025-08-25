import sys
import os
import torch
import pandas as pd
from tqdm import tqdm
import yaml
from Bio import SeqIO
import torch
import argparse
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import recall_score, precision_score, accuracy_score, f1_score
import torch.nn as nn
import numpy as np

import random

from Bio.SeqRecord import SeqRecord


evo_emb_path = "/opt/user/basspoom/GUI/model"
if evo_emb_path not in sys.path:
    sys.path.append(evo_emb_path)
        
from stripedhyena.tokenizer import CharLevelTokenizer
from stripedhyena.utils import dotdict
from model.evo_emb import StripedHyena_embedding


class DualPathModel(nn.Module):
    def __init__(self, embedding_dim):
        super(DualPathModel, self).__init__()

        self.scores_path = nn.Sequential(
            nn.Linear(24, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(),
            nn.Dropout(0.2)
        )

        self.embeddings_path = nn.Sequential(
            nn.Linear(embedding_dim, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(),
            nn.Dropout(0.2)
        )

        self.shared_representation = nn.Sequential(
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(),
            nn.Dropout(0.2)
        )
        
        self.classifier = nn.Sequential(
            nn.Linear(32, 1),
            nn.Sigmoid()
        )

    def forward(self, scores, embeddings):
        scores_embedding = self.scores_path(scores)
        embeddings_embedding = self.embeddings_path(embeddings)
        combined_embedding = torch.cat((scores_embedding, embeddings_embedding), dim=1)

        shared_embedding = self.shared_representation(combined_embedding)
        output = self.classifier(shared_embedding)
        return output


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


def read_fasta_biopython(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append([record.id, str(record.seq)])
    return pd.DataFrame(sequences, columns=['ID', 'Sequence'])


def run_filter(fasta_list, out_dir):
    for item in fasta_list:   
        file_path = f'{out_dir}/filtered_contigs/{item}.contigs.fasta'
        output_file = f'{out_dir}/filtered_contigs/{item}.contigs_8k.fasta'

        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(file_path, 'fasta'):
                if len(record.seq) > 8000:
                    start = random.randint(0, len(record.seq) - 8000)
                    new_seq = record.seq[start:start+8000]
                    new_record = SeqRecord(new_seq, id=record.id, description=record.description)
                    SeqIO.write(new_record, out_handle, 'fasta')
                else:
                    SeqIO.write(record, out_handle, 'fasta')




def load_evo_model_and_prepare(Evo_config, Evo_model, fasta_list, out_dir, device):
    
    if not os.path.exists(f"{out_dir}/evo"):
        os.makedirs(f"{out_dir}/evo") 

    config_path = f'{Evo_config}'
    state_dict_path = f'{Evo_model}'
    state_dict = torch.load(state_dict_path)

    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    global_config = dotdict(config, Loader=yaml.FullLoader)
    tokenizer = CharLevelTokenizer(512)
    model = StripedHyena_embedding(global_config)
    model.load_state_dict(state_dict, strict=True)

    model.to(device)
    model.to_bfloat16_except_poles_residues()
    model.eval()

    for item in fasta_list:   
            file_path = f'{out_dir}/filtered_contigs/{item}.contigs_8k.fasta'
            df = read_fasta_biopython(file_path)
            input_ids_list = []
            for sequence in tqdm(df['Sequence']):
                input_ids = torch.tensor(
                    tokenizer.tokenize(sequence),
                    dtype=torch.int,
                ).unsqueeze(0).tolist()
                input_ids_list.append(input_ids)
            df['input_ids'] = input_ids_list
            df.to_pickle(f'{out_dir}/filtered_contigs/{item}.contigs.pkl')
    
    print("The Evo model has been successfully loaded, and the input data is now fully prepared!")

    for item in fasta_list:   
        output_path = f'{out_dir}/evo/{item}.contigs.pkl'
        if os.path.exists(output_path):
            embedding_df = pd.read_pickle(output_path)
            completed_ids = set(embedding_df['ID'].tolist())
        else:
            embedding_df = pd.DataFrame(columns=['ID', 'embedding'])
            completed_ids = set()

        df = pd.read_pickle(f'{out_dir}/filtered_contigs/{item}.contigs.pkl')
        ids = df['ID'].tolist()  
        embeddings = []

        batch_size = 200
        for input_id, original_id in tqdm(zip(df['input_ids'], ids), total=len(ids)):
            if original_id in completed_ids:
                continue  
            
            with torch.no_grad():
                logits, _ = model(torch.tensor(input_id).to(device))  # (batch, length, vocab)
                embedding = logits.mean(dim=1)[0].cpu().tolist()  
            embeddings.append({'ID': original_id, 'embedding': embedding})

            if len(embeddings) >= batch_size:
                new_embeddings_df = pd.DataFrame(embeddings)
                embedding_df = pd.concat([embedding_df, new_embeddings_df], ignore_index=True)
                embedding_df.to_pickle(output_path)
                embeddings.clear()  
                torch.cuda.empty_cache()

        if embeddings:
            new_embeddings_df = pd.DataFrame(embeddings)
            embedding_df = pd.concat([embedding_df, new_embeddings_df], ignore_index=True)
            embedding_df.to_pickle(output_path)

        print(f"Embeddings saved successfully to '{out_dir}/evo/{item}.contigs.pkl'.")


def combine_and_predict(wpath_model,fasta_list, out_dir, device):
    
    if not os.path.exists(f"{out_dir}/viral_score_and_embedding"):
        os.makedirs(f"{out_dir}/viral_score_and_embedding") 
    
    for item in fasta_list:  
        df = pd.read_csv(f'{out_dir}/summary/{item}_summary.tsv', sep='\t')
        score_columns = [col for col in df.columns if 'score' in col.lower() and 'label' not in col.lower()]
        ids = df['Sequence ID']
        scores = df[score_columns]
        id_to_score_label = dict(zip(ids, scores.values.tolist()))
        df2 = pd.read_pickle(f'{out_dir}/evo/{item}.contigs.pkl')

        data = []
        for idx, row in tqdm(df2.iterrows(), total=df2.shape[0], desc="Processing rows"):
            id = row['ID']
            embedding = row['embedding']
            
            if id in id_to_score_label:
                score = id_to_score_label[id]
                data.append({'ID': id, 'embedding': embedding, 'scores': score})
                
        df_merged = pd.DataFrame(data)
        df_merged.to_pickle(f'{out_dir}/evo/{item}.contigs_merged.pkl')
        
    for item in fasta_list:
        embedding_dim = 4096
        model = DualPathModel(embedding_dim).to(device)
        model.load_state_dict(torch.load(f'{wpath_model}'))
        model.eval()
        
        test_dff = pd.read_pickle(f'{out_dir}/evo/{item}.contigs_merged.pkl')
        X_test_scores = torch.tensor(np.array(test_dff['scores'].tolist()), dtype=torch.float32).to(device)
        X_test_embeddings = torch.tensor(np.array(test_dff['embedding'].tolist()), dtype=torch.float32).to(device)
        test_ids = test_dff['ID'].tolist()
        
        dataset = TensorDataset(X_test_scores, X_test_embeddings)
        test_loader = DataLoader(dataset, batch_size=128, shuffle=False)
        
        all_scores = []
        all_predictions = []
        with torch.no_grad():
            for scores, embeddings in test_loader:
                outputs = model(scores, embeddings)
                scores = outputs.squeeze().cpu().numpy()
                all_scores.extend(scores)

                predictions = ["Virus" if score > 0.5 else "Not Virus" for score in scores]
                all_predictions.extend(predictions)
        
        with open(f"{out_dir}/viral_score_and_embedding/{item}_score.txt", "w") as f:
            f.write("ID\tScore\tPrediction\n")
            for id_, score, prediction in zip(test_ids, all_scores, all_predictions):
                f.write(f"{id_}\t{score}\t{prediction}\n")
                
                

def get_viralseq(input_dir, out_dir, fasta_list):

    dirs = [f"{out_dir}/viral_score_and_embedding", f"{out_dir}/labels", f"{out_dir}/mVC_evo", f"{out_dir}/mNVC_evo"]

    for dir_path in dirs:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            
    for item in fasta_list:
        os.system(f"grep '>' {input_dir}/{item}.contigs.fa > {out_dir}/viral_score_and_embedding/temp1.txt")
        os.system(f"awk -F ' ' '{{print $1}}' {out_dir}/viral_score_and_embedding/temp1.txt > {out_dir}/viral_score_and_embedding/temp2.txt")
        os.system(f"sed 's/>//' {out_dir}/viral_score_and_embedding/temp2.txt > {out_dir}/labels/{item}_all_labels.txt")
        os.system(f'awk \'NR>1 && $3=="Virus"{{print $1}}\' {out_dir}/viral_score_and_embedding/{item}_score.txt > {out_dir}/labels/{item}_viral_labels.txt')
        os.system(f'grep -Fv -f {out_dir}/labels/{item}_viral_labels.txt {out_dir}/labels/{item}_all_labels.txt > {out_dir}/labels/{item}_non_viral_labels.txt')
        os.system(f'seqkit grep -w 0 -f {out_dir}/labels/{item}_viral_labels.txt  {input_dir}/{item}.contigs.fa > {out_dir}/mVC_evo/{item}.mVC.fasta')
        os.system(f'seqkit grep -w 0 -f {out_dir}/labels/{item}_non_viral_labels.txt  {input_dir}/{item}.contigs.fa > {out_dir}/mNVC_evo/{item}.mNVC.fasta')

    # print("All processes have been done!")
    print(f"Your viral contigs are in path '{out_dir}/mVC_evo/' ")
    print(f"Your non-viral contigs are in path '{out_dir}/mNVC_evo/'")




def check_viralseq(out_dir):

    import os
    import glob
    import subprocess
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from io import StringIO

    for mode in ['mVC_evo', 'mNVC_evo']:
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
                
                bins = list(range(0, 40001, 100))
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
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
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
            ax.set_xlim(0, 40000)
            ax.set_xticks(np.arange(0, 40001, 5000))  
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
    # ANSI escape codes for colors
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "This script is for the viral and non-viral contigs prediction by VPAC and Evo" + RESET)
    print("\n" + PURPLE + "Examples:" + RESET)
    print(" VPAC-evo -i /data/input -o /data/output --device cuda:0 --Evo_path /opt/user/xxx/evo-model --Evo_model /opt/user/xxx/model/evo-1-8k.pt --Evo_config /opt/user/xxx/model/evo-1-8k-base_inference.yml --wpath_model /opt/user/xxx/w-model/epoch-991.pt")
    print("\n" + GREEN + "The -i parameter specifies the directory containing input FASTA files." + RESET)
    print("  " + BLUE + "/path/to/input/sample1.fasta" + RESET)
    print("  " + BLUE + "/path/to/input/sample2.fasta" + RESET)
    print("\n" + GREEN + "The -o parameter specifies the output directory for the results." + RESET)
    print("\n" + GREEN + "The --device parameter specifies the GPU device to use (default: 'cuda:0')." + RESET)
    print("\n" + GREEN + "The --Evo_path parameter specifies the location of Evo's model.py file." + RESET)
    print("\n" + GREEN + "The --Evo_model parameter specifies the location of Evo's model file." + RESET)
    print("\n" + GREEN + "The --Evo_config parameter specifies the location of Evo's config file." + RESET)
    print("\n" + GREEN + "The --wpath_model parameter specifies the location of the 2-path model file." + RESET)
    print("\n" + PURPLE + "Detailed parameters")
    print("" + GREEN + "=" * 50 + " " + RESET + YELLOW + "Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)
    print("\n" + YELLOW + "Global parameters:" + RESET)
    print(f"  {'-i, --input_dir':<40} Directory containing the input fasta files.")
    print(f"  {'-o, --out_dir':<40} Output directory for results.")
    print(f"  {'-d, --device':<40} GPU device to use (default: 'cuda:0').")
    print(f"  {'--Evo_path':<40} Location of Evo's model.py file.")
    print(f"  {'--Evo_model':<40} Location of Evo's model file.")
    print(f"  {'--Evo_config':<40} Location of Evo's config file.")
    print(f"  {'--wpath_model':<40} Location of the 2-path model file.")

class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()


def main():
    # parser = argparse.ArgumentParser(description="This script integrates DeepVirFinder, VirSorter2, and VIBRANT to identify viral contigs, then run kaiju and CheckV to check quality.", formatter_class=argparse.RawTextHelpFormatter)
    parser = CustomArgumentParser(description="Predict viral and non-viral contigs by VPAC + Evo.")
  
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="The assembled contig (fasta) file directory, such as /xx/raw_contigs (raw_contigs/sample1.fa, raw_contigs/sample2.fa ...)")
    parser.add_argument("-o", "--out_dir", type=str,  required=True, help="The location of viral prediction results and middle files")
    parser.add_argument("-d", "--device", type=str,  required=True, help="Your GPU name, default 'cuda:0'. ")

    # Evo
    # parser.add_argument("-Evo_path", "--Evo_path", type=str, required=True, help="The location of Evo's model.py, for example '/opt/user/xxx/evo-model' ")
    parser.add_argument("-Evo_model", "--Evo_model", type=str, required=True, help="The location of Evo's model file, for example '/opt/user/xxx/model/evo-1-8k.pt' ")
    parser.add_argument("-Evo_config", "--Evo_config", type=str, required=True, help="The location of Evo's config file, for example '/opt/user/xxx/model/evo-1-8k-base_inference.yml' ")


    parser.add_argument("-wpath_model", "--wpath_model", type=str, required=True, help="The location of 2 papth model file, for example '/opt/user/xxx/w-model/epoch-991.pt' ")


    args = parser.parse_args()


    fasta_list = Getfasta_list(args.input_dir)

    run_filter(fasta_list, args.out_dir)

    load_evo_model_and_prepare(args.Evo_config, args.Evo_model, fasta_list, args.out_dir, args.device)
    combine_and_predict(args.wpath_model,fasta_list, args.out_dir, args.device)
    get_viralseq(args.input_dir, args.out_dir, fasta_list)
    check_viralseq(args.out_dir)


if __name__ == "__main__":
    main()