import os
import sys
import yaml
import json
import argparse
import random
import numpy as np
import pandas as pd
import torch

from tqdm import tqdm
from pathlib import Path
from evo import Evo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out-pkl", required=True)
    parser.add_argument("--config-yml", required=True)
    return parser.parse_args()


def load_config(config_yml):
    with open(config_yml, "r") as f:
        return yaml.safe_load(f)


def read_fasta(fasta_path):
    records = []
    seq_id = None
    seq_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(seq_parts).upper()))
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if seq_id is not None:
            records.append((seq_id, "".join(seq_parts).upper()))

    return records


def clean_sequence(seq):
    seq = seq.upper().replace("U", "T")
    allowed = set("ACGTN")
    return "".join([base for base in seq if base in allowed])


def split_sequence(seq, max_len):
    if len(seq) <= max_len:
        return [seq]

    chunks = []
    for start in range(0, len(seq), max_len):
        chunk = seq[start:start + max_len]
        if chunk:
            chunks.append(chunk)
    return chunks


def resolve_evo_model_name_or_path(config):
    evo_models_dir = Path(config["Evo_models"]).expanduser().resolve()

    candidates = [
        evo_models_dir / "evo-1-131k-base",
        evo_models_dir / "evo-1-8k-base",
        evo_models_dir / "evo-1.5-8k-base"
    ]

    for path in candidates:
        if path.exists():
            return str(path)

    return "evo-1-131k-base"


def find_target_module(model, preferred_patterns):
    named_modules = dict(model.named_modules())

    for pattern in preferred_patterns:
        if pattern in named_modules:
            return named_modules[pattern], pattern

    for name, module in named_modules.items():
        lname = name.lower()
        if "blocks" in lname and ("mlp" in lname or "mixer" in lname):
            return module, name

    return None, None


def extract_last_hidden_by_hook(model, tokenizer, seq, device, layer_name=None):
    captured = {}

    if layer_name is not None:
        target_module = dict(model.named_modules()).get(layer_name)
        target_name = layer_name
    else:
        preferred_patterns = [
            "backbone.blocks.28.mlp.l3",
            "blocks.28.mlp.l3",
            "backbone.blocks.28",
            "blocks.28",
            "backbone.blocks.31",
            "blocks.31"
        ]
        target_module, target_name = find_target_module(model, preferred_patterns)

    if target_module is None:
        raise RuntimeError("Cannot find a valid Evo1 internal layer for embedding extraction")

    def hook_fn(module, inputs, output):
        if isinstance(output, tuple):
            output = output[0]
        captured["embedding"] = output.detach()

    handle = target_module.register_forward_hook(hook_fn)

    try:
        input_ids = torch.tensor(
            tokenizer.tokenize(seq),
            dtype=torch.int
        ).to(device).unsqueeze(0)

        with torch.no_grad():
            model(input_ids)

        if "embedding" not in captured:
            raise RuntimeError(f"No embedding captured from layer: {target_name}")

        emb = captured["embedding"]

        if emb.ndim == 3:
            emb = emb[0]
        elif emb.ndim == 2:
            pass
        else:
            raise RuntimeError(f"Unexpected embedding tensor shape: {tuple(emb.shape)}")

        emb = emb.float().cpu()
        seq_embedding = emb.mean(dim=0).numpy().astype(np.float32)

    finally:
        handle.remove()

    return seq_embedding, target_name


def get_sequence_embedding(model, tokenizer, seq, device, max_tokens, layer_name=None):
    chunks = split_sequence(seq, max_tokens)
    chunk_embeddings = []
    chunk_lengths = []
    used_layer_name = None

    for chunk in chunks:
        emb, used_layer_name = extract_last_hidden_by_hook(
            model=model,
            tokenizer=tokenizer,
            seq=chunk,
            device=device,
            layer_name=layer_name
        )
        chunk_embeddings.append(emb)
        chunk_lengths.append(len(chunk))

    chunk_embeddings = np.vstack(chunk_embeddings).astype(np.float32)
    chunk_lengths = np.array(chunk_lengths, dtype=np.float32)
    weights = chunk_lengths / chunk_lengths.sum()

    sequence_embedding = np.sum(chunk_embeddings * weights[:, None], axis=0).astype(np.float32)
    return sequence_embedding, len(chunks), used_layer_name


def main():
    args = parse_args()

    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    config = load_config(args.config_yml)

    evo_models_dir = Path(config["Evo_models"]).expanduser().resolve()
    os.environ["HF_HOME"] = str(evo_models_dir)
    os.environ["TRANSFORMERS_CACHE"] = str(evo_models_dir)
    os.environ["HF_HUB_CACHE"] = str(evo_models_dir)

    model_name_or_path = resolve_evo_model_name_or_path(config)

    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    max_tokens = 131072 if "131k" in str(model_name_or_path) else 8192

    records = read_fasta(args.fasta)
    if len(records) == 0:
        raise ValueError("No FASTA records found")

    evo_model = Evo(model_name_or_path)
    model = evo_model.model
    tokenizer = evo_model.tokenizer

    model.to(device)
    model.eval()

    rows = []

    for seq_id, seq in tqdm(records, desc="Generating Evo1 embeddings"):
        seq = clean_sequence(seq)

        if len(seq) == 0:
            continue

        embedding, n_chunks, used_layer_name = get_sequence_embedding(
            model=model,
            tokenizer=tokenizer,
            seq=seq,
            device=device,
            max_tokens=max_tokens,
            layer_name=None
        )

        rows.append({
            "ID": seq_id,
            "embedding": embedding,
            "sequence_length": len(seq),
            "n_chunks": n_chunks,
            "embedding_dim": int(embedding.shape[0]),
            "model_name_or_path": str(model_name_or_path),
            "layer_name": used_layer_name
        })

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    embedding_df = pd.DataFrame(rows)

    if len(embedding_df) == 0:
        raise ValueError("No valid embeddings generated")

    out_pkl = Path(args.out_pkl)
    out_pkl.parent.mkdir(parents=True, exist_ok=True)
    embedding_df.to_pickle(out_pkl)

    summary = {
        "input_fasta": args.fasta,
        "output_pkl": str(out_pkl),
        "config_yml": args.config_yml,
        "n_records": int(len(embedding_df)),
        "embedding_dim": int(embedding_df.iloc[0]["embedding_dim"]),
        "model_name_or_path": str(model_name_or_path),
        "device": device,
        "max_tokens": int(max_tokens),
        "layer_name": str(embedding_df.iloc[0]["layer_name"])
    }

    summary_path = out_pkl.with_suffix(".summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Saved pkl: {out_pkl}")
    print(f"Saved summary: {summary_path}")
    print(f"Records: {len(embedding_df)}")
    print(f"Embedding dim: {embedding_df.iloc[0]['embedding_dim']}")


if __name__ == "__main__":
    main()