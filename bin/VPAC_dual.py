import os, sys, argparse, yaml, subprocess, glob, random
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import recall_score, precision_score, accuracy_score, f1_score  # 保持兼容，不改
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.markdown import Markdown
import matplotlib.pyplot as plt
from io import StringIO
from typing import Union

CharLevelTokenizer = None
dotdict = None
StripedHyena_embedding = None


def read_yaml(p: str) -> dict:
    with open(p, "r") as f:
        return yaml.safe_load(f)

def get_from_cfg(cfg: dict, *keys, default=None):
    """支持顶层键或分组键，如 ('models','Evo_models')."""
    for k in keys:
        if isinstance(k, (list, tuple)):
            node = cfg
            ok = True
            for kk in k:
                if isinstance(node, dict) and kk in node:
                    node = node[kk]
                else:
                    ok = False
                    break
            if ok:
                return node
        else:
            if k in cfg:
                return cfg[k]
    return default

def build_paths_from_config(cfg: dict):
    vmp_dual = get_from_cfg(
        cfg,
        "VPAC_dual", "VMP_dual",
        ("environments", "VPAC_dual"),
        ("environments", "VMP_dual"),
    )
    if not vmp_dual:
        raise ValueError()
    vmp_dual = str(vmp_dual)

    evo_dir = get_from_cfg(cfg, "Evo_models", ("models", "Evo_models"))
    vpac_dir = get_from_cfg(cfg, "VPAC_models", ("models", "VPAC_models"))
    if not evo_dir or not vpac_dir:
        raise ValueError("Evo_models or VPAC_models not found in config.yml.")
    evo_dir, vpac_dir = str(evo_dir), str(vpac_dir)

    Evo_model = os.path.join(evo_dir, "evo-1-8k.pt")
    Evo_config = os.path.join(evo_dir, "evo-1-8k-base_inference.yml")
    wpath_model = os.path.join(vpac_dir, "VPAC-dual.pth")

    for p in (Evo_model, Evo_config, wpath_model):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required model/config file: {p}")

    return vmp_dual, Evo_model, Evo_config, wpath_model, evo_dir


def quote_bash(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"

def run_in_env(env_prefix: str, cmd: str, capture=False) -> Union[subprocess.CompletedProcess, int]:
    wrapper = f'conda run --no-capture-output --prefix "{env_prefix}" bash -lc {quote_bash(cmd)}'
    if capture:
        return subprocess.run(wrapper, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    else:
        return os.system(wrapper)

def import_evo_from_config(evo_dir: str):
    evo_emb_py = os.path.join(evo_dir, "evo_emb.py")
    if not os.path.exists(evo_emb_py):
        raise FileNotFoundError(f"evo_emb.py not found under Evo_models: {evo_dir}")

    if evo_dir not in sys.path:
        sys.path.insert(0, evo_dir)

    import evo_emb  

    try:
        import stripedhyena.model as sh_model
        from flash_attn.modules.mha import MHA as FlashMHA
        sh_model.MHA = FlashMHA
    except Exception as e:
        print(f"[WARN] Failed to patch stripedhyena.model.MHA: {e}")
    from stripedhyena.tokenizer import CharLevelTokenizer as _Tok
    from stripedhyena.utils import dotdict as _Dot
    from evo_emb import StripedHyena_embedding as _Model

    return _Tok, _Dot, _Model



class DualPathModel(nn.Module):
    def __init__(self, embedding_dim):
        super(DualPathModel, self).__init__()
        self.scores_path = nn.Sequential(
            nn.Linear(24, 32), nn.BatchNorm1d(32), nn.LeakyReLU(), nn.Dropout(0.2)
        )
        self.embeddings_path = nn.Sequential(
            nn.Linear(embedding_dim, 32), nn.BatchNorm1d(32), nn.LeakyReLU(), nn.Dropout(0.2)
        )
        self.shared_representation = nn.Sequential(
            nn.Linear(64, 32), nn.BatchNorm1d(32), nn.LeakyReLU(), nn.Dropout(0.2)
        )
        self.classifier = nn.Sequential(nn.Linear(32, 1), nn.Sigmoid())

    def forward(self, scores, embeddings):
        s = self.scores_path(scores)
        e = self.embeddings_path(embeddings)
        h = torch.cat((s, e), dim=1)
        h = self.shared_representation(h)
        return self.classifier(h)

def Getfasta_list(input_dir):
    fasta_list = []
    for i in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, i)):
            continue
        if 'contigs.fa' in i:
            fasta_list.append(i[: i.index('.')])
    print('Your sample list:', fasta_list)
    return fasta_list

def read_fasta_biopython(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append([record.id, str(record.seq)])
    return pd.DataFrame(sequences, columns=['ID', 'Sequence'])

def run_filter(fasta_list, out_dir):
    Path(f"{out_dir}/filtered_contigs").mkdir(parents=True, exist_ok=True)
    for item in fasta_list:
        fp = f'{out_dir}/filtered_contigs/{item}.contigs.fasta'
        in_guess = fp if os.path.exists(fp) else fp.replace(".fasta", "")
        if not os.path.exists(in_guess):
            continue
        file_path = in_guess
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
    global CharLevelTokenizer, dotdict, StripedHyena_embedding
    os.makedirs(f"{out_dir}/evo", exist_ok=True)

    state_dict = torch.load(Evo_model)
    with open(Evo_config, 'r') as f:
        config = yaml.safe_load(f)
    global_config = dotdict(config)

    tokenizer = CharLevelTokenizer(512)
    model = StripedHyena_embedding(global_config)
    model.load_state_dict(state_dict, strict=True)
    model.to(device)
    model.to_bfloat16_except_poles_residues()
    model.eval()

    for item in fasta_list:
        file_path = f'{out_dir}/filtered_contigs/{item}.contigs_8k.fasta'
        if not os.path.exists(file_path):
            continue
        df = read_fasta_biopython(file_path)
        input_ids_list = []
        for sequence in tqdm(df['Sequence']):
            input_ids = torch.tensor(tokenizer.tokenize(sequence), dtype=torch.int).unsqueeze(0).tolist()
            input_ids_list.append(input_ids)
        df['input_ids'] = input_ids_list
        df.to_pickle(f'{out_dir}/filtered_contigs/{item}.contigs.pkl')

    print("Evo model loaded; inputs prepared.")

    for item in fasta_list:
        output_path = f'{out_dir}/evo/{item}.contigs.pkl'
        if os.path.exists(output_path):
            embedding_df = pd.read_pickle(output_path)
            completed_ids = set(embedding_df['ID'].tolist())
        else:
            embedding_df = pd.DataFrame(columns=['ID', 'embedding'])
            completed_ids = set()

        src = f'{out_dir}/filtered_contigs/{item}.contigs.pkl'
        if not os.path.exists(src):
            continue
        df = pd.read_pickle(src)
        ids = df['ID'].tolist()
        embeddings = []
        batch_size = 200

        for input_id, original_id in tqdm(zip(df['input_ids'], ids), total=len(ids)):
            if original_id in completed_ids:
                continue
            with torch.no_grad():
                logits, _ = model(torch.tensor(input_id).to(device))
                embedding = logits.mean(dim=1)[0].cpu().tolist()
            embeddings.append({'ID': original_id, 'embedding': embedding})
            if len(embeddings) >= batch_size:
                embedding_df = pd.concat([embedding_df, pd.DataFrame(embeddings)], ignore_index=True)
                embedding_df.to_pickle(output_path)
                embeddings.clear()
                torch.cuda.empty_cache()

        if embeddings:
            embedding_df = pd.concat([embedding_df, pd.DataFrame(embeddings)], ignore_index=True)
            embedding_df.to_pickle(output_path)

        print(f"Embeddings saved: {output_path}")

def combine_and_predict(wpath_model, fasta_list, out_dir, device):
    os.makedirs(f"{out_dir}/viral_score_and_embedding", exist_ok=True)
    for item in fasta_list:
        summ = f'{out_dir}/summary/{item}_summary.tsv'
        evo_pkl = f'{out_dir}/evo/{item}.contigs.pkl'
        if not (os.path.exists(summ) and os.path.exists(evo_pkl)):
            continue
        df = pd.read_csv(summ, sep='\t')
        score_columns = [c for c in df.columns if 'score' in c.lower() and 'label' not in c.lower()]
        ids = df['Sequence ID']
        scores = df[score_columns]
        id_to_scores = dict(zip(ids, scores.values.tolist()))
        df2 = pd.read_pickle(evo_pkl)
        data = []
        for _, row in tqdm(df2.iterrows(), total=df2.shape[0], desc="Merging scores+embeddings"):
            _id = row['ID']
            emb = row['embedding']
            if _id in id_to_scores:
                data.append({'ID': _id, 'embedding': emb, 'scores': id_to_scores[_id]})
        df_merged = pd.DataFrame(data)
        df_merged.to_pickle(f'{out_dir}/evo/{item}.contigs_merged.pkl')

    for item in fasta_list:
        merged = f'{out_dir}/evo/{item}.contigs_merged.pkl'
        if not os.path.exists(merged):
            continue
        embedding_dim = 4096
        model = DualPathModel(embedding_dim).to(device)
        model.load_state_dict(torch.load(wpath_model))
        model.eval()

        test_dff = pd.read_pickle(merged)
        Xs = torch.tensor(np.array(test_dff['scores'].tolist()), dtype=torch.float32).to(device)
        Xe = torch.tensor(np.array(test_dff['embedding'].tolist()), dtype=torch.float32).to(device)
        ids = test_dff['ID'].tolist()
        loader = DataLoader(TensorDataset(Xs, Xe), batch_size=128, shuffle=False)

        all_scores, all_preds = [], []
        with torch.no_grad():
            for s, e in loader:
                out = model(s, e).squeeze().cpu().numpy()
                all_scores.extend(out)
                all_preds.extend(["Virus" if v > 0.5 else "Not Virus" for v in out])
        with open(f"{out_dir}/viral_score_and_embedding/{item}_score.txt", "w") as f:
            f.write("ID\tScore\tPrediction\n")
            for _id, sc, pr in zip(ids, all_scores, all_preds):
                f.write(f"{_id}\t{sc}\t{pr}\n")

def get_viralseq(input_dir, out_dir, fasta_list, env_prefix):
    for d in [f"{out_dir}/viral_score_and_embedding", f"{out_dir}/labels", f"{out_dir}/mVC_dual", f"{out_dir}/mNVC_dual"]:
        os.makedirs(d, exist_ok=True)

    for item in fasta_list:
        run_in_env(env_prefix, f"grep '>' {quote_bash(os.path.join(input_dir, item+'.contigs.fa'))} > {quote_bash(os.path.join(out_dir,'viral_score_and_embedding','temp1.txt'))}")
        run_in_env(env_prefix, f"awk -F ' ' '{{print $1}}' {quote_bash(os.path.join(out_dir,'viral_score_and_embedding','temp1.txt'))} > {quote_bash(os.path.join(out_dir,'viral_score_and_embedding','temp2.txt'))}")
        run_in_env(env_prefix, f"sed 's/>//' {quote_bash(os.path.join(out_dir,'viral_score_and_embedding','temp2.txt'))} > {quote_bash(os.path.join(out_dir,'labels',f'{item}_all_labels.txt'))}")
        run_in_env(env_prefix, f"awk 'NR>1 && $3==\"Virus\"{{print $1}}' {quote_bash(os.path.join(out_dir,'viral_score_and_embedding',f'{item}_score.txt'))} > {quote_bash(os.path.join(out_dir,'labels',f'{item}_viral_labels.txt'))}")
        run_in_env(env_prefix, f"grep -Fv -f {quote_bash(os.path.join(out_dir,'labels',f'{item}_viral_labels.txt'))} {quote_bash(os.path.join(out_dir,'labels',f'{item}_all_labels.txt'))} > {quote_bash(os.path.join(out_dir,'labels',f'{item}_non_viral_labels.txt'))}")
        run_in_env(env_prefix, f"seqkit grep -w 0 -f {quote_bash(os.path.join(out_dir,'labels',f'{item}_viral_labels.txt'))} {quote_bash(os.path.join(input_dir,item+'.contigs.fa'))} > {quote_bash(os.path.join(out_dir,'mVC_dual',f'{item}.mVC.fasta'))}")
        run_in_env(env_prefix, f"seqkit grep -w 0 -f {quote_bash(os.path.join(out_dir,'labels',f'{item}_non_viral_labels.txt'))} {quote_bash(os.path.join(input_dir,item+'.contigs.fa'))} > {quote_bash(os.path.join(out_dir,'mNVC_dual',f'{item}.mNVC.fasta'))}")

    print(f"Your viral contigs are in '{out_dir}/mVC_dual/'")
    print(f"Your non-viral contigs are in '{out_dir}/mNVC_dual/'")

def check_viralseq(out_dir, env_prefix):
    for mode in ['mVC_dual', 'mNVC_dual']:
        contigs_dir = os.path.join(out_dir, mode)
        stat_dir = os.path.join(contigs_dir, 'contigs_stat')
        os.makedirs(stat_dir, exist_ok=True)

        fa_files = glob.glob(os.path.join(contigs_dir, "*.fasta"))
        if not fa_files:
            print(f"No contig files found in {contigs_dir} for statistics.")
            continue

        cmd = f"seqkit stat -a {' '.join(map(quote_bash, fa_files))}"
        proc = run_in_env(env_prefix, cmd, capture=True)
        if proc.stdout is None:
            print(f"[{mode}] Failed to run seqkit stat.")
            continue
        stats_output = proc.stdout
        try:
            stats_df = pd.read_csv(StringIO(stats_output), sep='\t')
        except Exception as e:
            print(f"[{mode}] Could not parse seqkit output: {e}")
            continue
        stats_summary_path = os.path.join(stat_dir, 'contig_stats_summary.tsv')
        stats_df.to_csv(stats_summary_path, sep='\t', index=False)
        print(f"[{mode}] Contig statistics saved to {stats_summary_path}")

        all_samples_data = []
        for fp in fa_files:
            sample_name = os.path.basename(fp).replace('.contigs.fa', '').replace('.fasta','')
            cmd_len = f"seqkit fx2tab -n -l {quote_bash(fp)} | cut -f 2"
            p2 = run_in_env(env_prefix, cmd_len, capture=True)
            if p2.stdout is None:
                continue
            try:
                lengths = [int(x) for x in p2.stdout.strip().splitlines() if x.strip()]
                bins = list(range(0, 40001, 100))
                counts, bin_edges = np.histogram(lengths, bins=bins)
                all_samples_data.append((sample_name, counts, bin_edges))
            except Exception:
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
            row, col = divmod(idx, n_cols)
            ax = axes[row, col] if n_rows > 1 else (axes[col] if n_cols > 1 else axes)
            x = bin_edges[:-1]
            ax.bar(x, counts, width=100, align='edge', alpha=0.7, edgecolor='k', linewidth=0.3)
            ax.set_title(sample_name, fontsize=10)
            ax.set_xlabel('Length (bp)', fontsize=8)
            ax.set_ylabel('Count', fontsize=8)
            ax.set_xlim(0, 40000)
            ax.set_xticks(np.arange(0, 40001, 5000))
            ax.tick_params(axis='both', labelsize=8)

        for idx in range(n_samples, n_rows * n_cols):
            row, col = divmod(idx, n_cols)
            if n_rows > 1:
                axes[row, col].axis('off')
            else:
                (axes[col] if n_cols > 1 else axes).axis('off')

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(stat_dir, 'contig_length_distribution.pdf'), bbox_inches='tight')
        fig.savefig(os.path.join(stat_dir, 'contig_length_distribution.png'), dpi=500, bbox_inches='tight')
        plt.close()
        print(f"[{mode}] Plots saved.")


def custom_help():
    console = Console()
    intro = "The 'VPAC-evo' script predicts viral vs non-viral contigs using Evo embeddings plus a 2-path VPAC classifier."
    examples_md = Markdown(
        "**Example:**\n```bash\n"
        "VPAC-evo -i /data/contigs -o /data/vpac_evo_results -d cuda:0 "
        "-cf /path/to/config.yml\n```"
    )
    fasta_md = Markdown("**Input (-i / --input_dir):** directory with per-sample FASTA files.")
    out_md = Markdown("**Output (-o / --out_dir):** root folder for results.")

    console.print(Panel(intro, border_style="cyan", title="VPAC-evo", title_align="left"))
    console.print(examples_md)
    console.print(fasta_md)
    console.print(out_md)
    console.print()

    out_tbl = Table(show_header=True, header_style="bold blue")
    out_tbl.add_column("Files/Dirs", style="cyan", no_wrap=True)
    out_tbl.add_column("Description", style="white")
    out_tbl.add_row("~/outdir/mNVC_dual/", "VPAC dual-path classifier predicted non-viral contigs.")
    out_tbl.add_row("~/outdir/mVC_dual/",  "VPAC dual-path classifier predicted viral contigs.")

    console.print(Panel(out_tbl, border_style="blue", title="Outputs", title_align="left"))
    console.print(Panel("Models are auto-resolved from the config file; you mainly specify input/output and device.", border_style="cyan"))
    console.print()

    console.print("[bold]Detailed parameters[/bold]\n")
    console.print(Panel("Usage: VPAC-evo [OPTIONS]\n"
                        "1) Provide contig FASTA directory\n"
                        "2) Provide config (-cf) which sets env & model paths\n"
                        "3) Export viral-only contigs and prediction table.",
                        border_style="cyan", title="Usage & Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-i, --input_dir", "Directory containing input contig FASTA files (per sample).")
    g_tbl.add_row("-o, --out_dir",
    "Output directory for all results and intermediates.\n"
    "[bold yellow]Important:[/bold yellow] This must be the SAME output directory used by the VPAC single-path classifier.\n"
    "The script loads per-contig scoring files from the VPAC single-path classifier's output (e.g., the 'summary/' folder).")

    g_tbl.add_row("-d, --device","Compute device, e.g. 'cuda:0' or 'cpu' (default: 'cuda:0').")
    g_tbl.add_row("-cf, --config","Path to config.yml (provides VPAC_dual environment prefix and Evo_models / VPAC_models paths).")
    
    console.print(g_tbl)
    console.print()
    console.print(Rule(style="dim"))

class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="Predict viral and non-viral contigs using Evo embeddings and a VPAC dual-path classifier.")
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    parser.add_argument("-d", "--device", type=str, default="cuda:0")
    parser.add_argument("-cf", "--config", type=str, required=True)
    args = parser.parse_args()

    cfg = read_yaml(args.config)
    VMP_DUAL, Evo_model, Evo_config, wpath_model, evo_dir = build_paths_from_config(cfg)

    global CharLevelTokenizer, dotdict, StripedHyena_embedding
    CharLevelTokenizer, dotdict, StripedHyena_embedding = import_evo_from_config(evo_dir)

    fasta_list = Getfasta_list(args.input_dir)
    run_filter(fasta_list, args.out_dir)
    load_evo_model_and_prepare(Evo_config, Evo_model, fasta_list, args.out_dir, args.device)
    combine_and_predict(wpath_model, fasta_list, args.out_dir, args.device)

    get_viralseq(args.input_dir, args.out_dir, fasta_list, VMP_DUAL)
    check_viralseq(args.out_dir, VMP_DUAL)

if __name__ == "__main__":
    main()
