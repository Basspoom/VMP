# VPAC Training Scripts README

This directory contains the scripts and input files used to train and validate VPAC classifiers.

## Input data

| File | Description |
|---|---|
| `train_shff.tsv` | Training table. The first column is the sequence or contig ID, columns 2–25 are the 24 VPAC scoring features, and column 26 is the binary label. |
| `vali.tsv` | Independent validation table with the same format as `train_shff.tsv`. |

The label column should contain:

```text
Virus
Non-Virus
```

All training scripts read `train_shff.tsv` and `vali.tsv` from the current working directory.  
All model outputs, metrics, prediction tables, ROC curves, and confusion matrices are saved to:

```text
./vali/
```

## Single-path classifier scripts

The following scripts train VPAC single-path classifiers using only the 24-dimensional VPAC scoring feature vector.

| Script | Model | Function |
|---|---|---|
| `VPAC-CNN.py` | 1D Convolutional Neural Network | Treats the 24 features as a one-dimensional sequence and learns local feature patterns with convolutional layers. |
| `VPAC-FNN.py` | Feedforward Neural Network | Trains a multilayer fully connected classifier using the 24 VPAC features. |
| `VPAC-AE.py` | Autoencoder-based classifier | Learns a denoised representation of the 24 features before classification. |
| `VPAC-VAE.py` | Variational Autoencoder-based classifier | Uses one-hot encoded feature scores and latent probabilistic representation learning for classification. |
| `VPAC-KAN.py` | Kolmogorov-Arnold Network | Uses spline-based KAN layers to model nonlinear relationships among VPAC features. |
| `VPAC-GB.py` | XGBoost / Gradient Boosting | Trains an XGBoost classifier using randomized hyperparameter search and the selected best parameter set. |
| `VPAC-RF.py` | Random Forest | Trains a random forest classifier using randomized hyperparameter search and the selected best parameter set. |
| `VPAC-SVM.py` | Support Vector Machine | Trains an SVM classifier with feature standardization and randomized hyperparameter search. |

Run examples:

```bash
python VPAC-CNN.py
python VPAC-FNN.py
python VPAC-AE.py
python VPAC-VAE.py
python VPAC-KAN.py
python VPAC-GB.py
python VPAC-RF.py
python VPAC-SVM.py
```

## Dual-path classifier script

| Script | Model | Function |
|---|---|---|
| `VPAC-dual.py` | VPAC dual-path classifier | Integrates the 24 VPAC scoring features and Evo sequence embeddings. The VPAC feature branch and Evo embedding branch are separately projected to 64 dimensions, concatenated, compressed to a 32-dimensional shared representation, and used for binary classification. |

Additional required input for `VPAC-dual.py`:

```text
evo_embeddings.pkl
```

The pkl file should contain at least:

| Column | Description |
|---|---|
| `ID` | Sequence or contig ID matching the first column of `train_shff.tsv` and `vali.tsv`. |
| `embedding` | Evo embedding vector for the corresponding sequence. |

Run example:

```bash
python VPAC-dual.py
```

## Evo embedding generation script

| Script | Function |
|---|---|
| `Embedding_and_feature.py` | Generates Evo embeddings from a FASTA file and saves them as a pkl file for dual-path VPAC training. |

Required inputs:

```text
input.fasta
config.yml
```

Example command:

```bash
python Embedding_and_feature.py \
  --fasta input.fasta \
  --out-pkl evo_embeddings.pkl \
  --config-yml config.yml
```

Example `config.yml` model section:

```yaml
VPAC_models: '/path/to/VMP/dependency/models/VPAC_models'
Evo_models: '/path/to/VMP/dependency/models/Evo'
```

Expected output:

```text
evo_embeddings.pkl
evo_embeddings.summary.json
```

## Recommended workflow

1. Prepare `train_shff.tsv` and `vali.tsv`.
2. Train single-path classifiers as needed.
3. Generate Evo embeddings from FASTA using `Embedding_and_feature.py`.
4. Train the dual-path classifier with `VPAC-dual.py`.

## Main output files

Most scripts generate the following types of outputs in `./vali/`:

| Output type | Description |
|---|---|
| Model files | Trained model checkpoints, such as `.pth`, `.pt`, `.pkl`, or `.keras`. |
| Metrics tables | Training and validation metrics in `.csv` or `.json` format. |
| Prediction tables | Per-sequence validation predictions, including true labels, predicted labels, and virus probabilities. |
| ROC curves | Model ROC curve PDFs. |
| Confusion matrices | Confusion matrix PDFs and related summary tables. |
