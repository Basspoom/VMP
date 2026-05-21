import os
import json
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import seaborn as sns

from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef,
    roc_curve,
    auc,
    confusion_matrix
)
from tqdm import tqdm

seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
embedding_pkl_path = "evo_embeddings.pkl"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

epochs = 1500
batch_size = 128
learning_rate = 0.01
min_learning_rate = 0.001
dropout_rate = 0.3

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")
embedding_df = pd.read_pickle(embedding_pkl_path)

if "ID" not in embedding_df.columns:
    raise ValueError("evo_embeddings.pkl must contain column: ID")

if "embedding" not in embedding_df.columns:
    raise ValueError("evo_embeddings.pkl must contain column: embedding")

embedding_df = embedding_df[["ID", "embedding"]].copy()
embedding_df["ID"] = embedding_df["ID"].astype(str)
embedding_df = embedding_df.drop_duplicates(subset=["ID"], keep="first")

train_id_col = train_df.columns[0]
vali_id_col = vali_df.columns[0]

train_df[train_id_col] = train_df[train_id_col].astype(str)
vali_df[vali_id_col] = vali_df[vali_id_col].astype(str)

train_merged = train_df.merge(
    embedding_df,
    left_on=train_id_col,
    right_on="ID",
    how="inner"
)

vali_merged = vali_df.merge(
    embedding_df,
    left_on=vali_id_col,
    right_on="ID",
    how="inner"
)

if len(train_merged) != len(train_df):
    missing_train = set(train_df[train_id_col]) - set(train_merged[train_id_col])
    raise ValueError(f"Missing embeddings for {len(missing_train)} training records")

if len(vali_merged) != len(vali_df):
    missing_vali = set(vali_df[vali_id_col]) - set(vali_merged[vali_id_col])
    raise ValueError(f"Missing embeddings for {len(missing_vali)} validation records")

X_train_scores = train_merged.iloc[:, 1:25].values.astype(np.float32)
X_val_scores = vali_merged.iloc[:, 1:25].values.astype(np.float32)

X_train_embeddings = np.vstack(train_merged["embedding"].values).astype(np.float32)
X_val_embeddings = np.vstack(vali_merged["embedding"].values).astype(np.float32)

if X_train_scores.shape[1] != 24:
    raise ValueError(f"Expected 24 score features, got {X_train_scores.shape[1]}")

if X_train_embeddings.shape[1] != 4096:
    raise ValueError(f"Expected 4096-dimensional Evo embeddings, got {X_train_embeddings.shape[1]}")

def encode_labels(series):
    labels = series.astype(str).str.strip().str.lower()
    return (labels == "virus").astype(np.float32).values

y_train = encode_labels(train_merged.iloc[:, 25])
y_val = encode_labels(vali_merged.iloc[:, 25])

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

X_train_scores_tensor = torch.tensor(X_train_scores, dtype=torch.float32)
X_train_embeddings_tensor = torch.tensor(X_train_embeddings, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train, dtype=torch.float32)

X_val_scores_tensor = torch.tensor(X_val_scores, dtype=torch.float32)
X_val_embeddings_tensor = torch.tensor(X_val_embeddings, dtype=torch.float32)
y_val_tensor = torch.tensor(y_val, dtype=torch.float32)

train_dataset = TensorDataset(
    X_train_scores_tensor,
    X_train_embeddings_tensor,
    y_train_tensor
)

val_dataset = TensorDataset(
    X_val_scores_tensor,
    X_val_embeddings_tensor,
    y_val_tensor
)

drop_last_train = len(train_dataset) % batch_size == 1 and len(train_dataset) > batch_size

train_loader = DataLoader(
    train_dataset,
    batch_size=batch_size,
    shuffle=True,
    drop_last=drop_last_train
)

val_loader = DataLoader(
    val_dataset,
    batch_size=batch_size,
    shuffle=False,
    drop_last=False
)

class DualPathModel(nn.Module):
    def __init__(self, embedding_dim=4096, dropout_rate=0.3):
        super().__init__()

        self.scores_path = nn.Sequential(
            nn.Linear(24, 64),
            nn.BatchNorm1d(64),
            nn.LeakyReLU(),
            nn.Dropout(dropout_rate)
        )

        self.embeddings_path = nn.Sequential(
            nn.Linear(embedding_dim, 64),
            nn.BatchNorm1d(64),
            nn.LeakyReLU(),
            nn.Dropout(dropout_rate)
        )

        self.integration = nn.Sequential(
            nn.Linear(128, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(),
            nn.Dropout(dropout_rate)
        )

        self.classifier = nn.Sequential(
            nn.Linear(32, 1),
            nn.Sigmoid()
        )

    def forward(self, scores, embeddings):
        score_repr = self.scores_path(scores)
        embedding_repr = self.embeddings_path(embeddings)
        fused = torch.cat([score_repr, embedding_repr], dim=1)
        shared = self.integration(fused)
        output = self.classifier(shared).view(-1)
        return output

model = DualPathModel(
    embedding_dim=X_train_embeddings.shape[1],
    dropout_rate=dropout_rate
).to(device)

criterion = nn.BCELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
    optimizer,
    T_max=epochs,
    eta_min=min_learning_rate
)

metrics_history = {
    "epoch": [],
    "learning_rate": [],
    "train_loss": [],
    "val_loss": [],
    "train_accuracy": [],
    "val_accuracy": [],
    "train_precision": [],
    "val_precision": [],
    "train_recall": [],
    "val_recall": [],
    "train_f1": [],
    "val_f1": [],
    "train_mcc": [],
    "val_mcc": []
}

best_val_loss = float("inf")
best_val_f1 = -1.0

for epoch in tqdm(range(epochs), desc="Training"):
    model.train()

    train_loss_sum = 0.0
    train_probs = []
    train_labels = []

    for scores, embeddings, labels in train_loader:
        scores = scores.to(device)
        embeddings = embeddings.to(device)
        labels = labels.to(device)

        optimizer.zero_grad()
        outputs = model(scores, embeddings)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        train_loss_sum += loss.item() * scores.size(0)
        train_probs.extend(outputs.detach().cpu().numpy())
        train_labels.extend(labels.detach().cpu().numpy())

    train_loss = train_loss_sum / len(train_loader.dataset)
    train_probs = np.array(train_probs)
    train_labels = np.array(train_labels).astype(int)
    train_preds = (train_probs >= 0.5).astype(int)

    train_accuracy = accuracy_score(train_labels, train_preds)
    train_precision = precision_score(train_labels, train_preds, zero_division=0)
    train_recall = recall_score(train_labels, train_preds, zero_division=0)
    train_f1 = f1_score(train_labels, train_preds, zero_division=0)
    train_mcc = matthews_corrcoef(train_labels, train_preds)

    model.eval()

    val_loss_sum = 0.0
    val_probs = []
    val_labels = []

    with torch.no_grad():
        for scores, embeddings, labels in val_loader:
            scores = scores.to(device)
            embeddings = embeddings.to(device)
            labels = labels.to(device)

            outputs = model(scores, embeddings)
            loss = criterion(outputs, labels)

            val_loss_sum += loss.item() * scores.size(0)
            val_probs.extend(outputs.detach().cpu().numpy())
            val_labels.extend(labels.detach().cpu().numpy())

    val_loss = val_loss_sum / len(val_loader.dataset)
    val_probs = np.array(val_probs)
    val_labels = np.array(val_labels).astype(int)
    val_preds = (val_probs >= 0.5).astype(int)

    val_accuracy = accuracy_score(val_labels, val_preds)
    val_precision = precision_score(val_labels, val_preds, zero_division=0)
    val_recall = recall_score(val_labels, val_preds, zero_division=0)
    val_f1 = f1_score(val_labels, val_preds, zero_division=0)
    val_mcc = matthews_corrcoef(val_labels, val_preds)

    current_lr = optimizer.param_groups[0]["lr"]

    metrics_history["epoch"].append(epoch + 1)
    metrics_history["learning_rate"].append(current_lr)
    metrics_history["train_loss"].append(train_loss)
    metrics_history["val_loss"].append(val_loss)
    metrics_history["train_accuracy"].append(train_accuracy)
    metrics_history["val_accuracy"].append(val_accuracy)
    metrics_history["train_precision"].append(train_precision)
    metrics_history["val_precision"].append(val_precision)
    metrics_history["train_recall"].append(train_recall)
    metrics_history["val_recall"].append(val_recall)
    metrics_history["train_f1"].append(train_f1)
    metrics_history["val_f1"].append(val_f1)
    metrics_history["train_mcc"].append(train_mcc)
    metrics_history["val_mcc"].append(val_mcc)

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), os.path.join(outdir, "vpac_dual_path_best_loss_model.pth"))

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        torch.save(model.state_dict(), os.path.join(outdir, "vpac_dual_path_best_f1_model.pth"))

    if (epoch + 1) % 50 == 0 or epoch == 0:
        torch.save(model.state_dict(), os.path.join(outdir, f"vpac_dual_path_epoch_{epoch + 1}.pth"))

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"LR: {current_lr:.6f}, "
        f"Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}, "
        f"Train Acc: {train_accuracy:.6f}, Val Acc: {val_accuracy:.6f}, "
        f"Train Precision: {train_precision:.6f}, Val Precision: {val_precision:.6f}, "
        f"Train Recall: {train_recall:.6f}, Val Recall: {val_recall:.6f}, "
        f"Train F1: {train_f1:.6f}, Val F1: {val_f1:.6f}, "
        f"Train MCC: {train_mcc:.6f}, Val MCC: {val_mcc:.6f}"
    )

    scheduler.step()

torch.save(model.state_dict(), os.path.join(outdir, "vpac_dual_path_final_model.pth"))

model_config = {
    "score_input_dim": 24,
    "embedding_input_dim": int(X_train_embeddings.shape[1]),
    "score_branch_dim": 64,
    "embedding_branch_dim": 64,
    "fused_dim": 128,
    "shared_dim": 32,
    "dropout_rate": dropout_rate,
    "learning_rate": learning_rate,
    "min_learning_rate": min_learning_rate,
    "epochs": epochs,
    "batch_size": batch_size,
    "optimizer": "Adam",
    "scheduler": "CosineAnnealingLR",
    "loss": "BCELoss"
}

with open(os.path.join(outdir, "vpac_dual_path_model_config.json"), "w") as f:
    json.dump(model_config, f, indent=2)

metrics_df = pd.DataFrame(metrics_history)
metrics_df.to_csv(os.path.join(outdir, "vpac_dual_path_training_metrics.csv"), index=False)

model.eval()

train_probs = []
train_labels = []

with torch.no_grad():
    for scores, embeddings, labels in train_loader:
        scores = scores.to(device)
        embeddings = embeddings.to(device)
        outputs = model(scores, embeddings)
        train_probs.extend(outputs.detach().cpu().numpy())
        train_labels.extend(labels.numpy())

train_probs = np.array(train_probs)
train_labels = np.array(train_labels).astype(int)
train_preds = (train_probs >= 0.5).astype(int)

val_probs = []
val_labels = []

with torch.no_grad():
    for scores, embeddings, labels in val_loader:
        scores = scores.to(device)
        embeddings = embeddings.to(device)
        outputs = model(scores, embeddings)
        val_probs.extend(outputs.detach().cpu().numpy())
        val_labels.extend(labels.numpy())

val_probs = np.array(val_probs)
val_labels = np.array(val_labels).astype(int)
val_preds = (val_probs >= 0.5).astype(int)

final_metrics = {
    "train_accuracy": float(accuracy_score(train_labels, train_preds)),
    "val_accuracy": float(accuracy_score(val_labels, val_preds)),
    "train_precision": float(precision_score(train_labels, train_preds, zero_division=0)),
    "val_precision": float(precision_score(val_labels, val_preds, zero_division=0)),
    "train_recall": float(recall_score(train_labels, train_preds, zero_division=0)),
    "val_recall": float(recall_score(val_labels, val_preds, zero_division=0)),
    "train_f1": float(f1_score(train_labels, train_preds, zero_division=0)),
    "val_f1": float(f1_score(val_labels, val_preds, zero_division=0)),
    "train_mcc": float(matthews_corrcoef(train_labels, train_preds)),
    "val_mcc": float(matthews_corrcoef(val_labels, val_preds))
}

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "vpac_dual_path_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "vpac_dual_path_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_merged.copy()
vali_pred_df["true_label"] = val_labels
vali_pred_df["virus_probability"] = val_probs
vali_pred_df["pred_label"] = val_preds
vali_pred_df["pred_class"] = np.where(val_preds == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "vpac_dual_path_validation_predictions.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(metrics_history["epoch"], metrics_history["train_loss"], label="Train Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_loss"], label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("VPAC Dual-Path Loss")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vpac_dual_path_loss_curve.pdf"), format="pdf")
plt.close()

plt.figure(figsize=(8, 6))
plt.plot(metrics_history["epoch"], metrics_history["learning_rate"], label="Learning Rate")
plt.xlabel("Epoch")
plt.ylabel("Learning Rate")
plt.title("Cosine Annealing Learning Rate")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vpac_dual_path_learning_rate_curve.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(val_labels, val_probs)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "vpac_dual_path_roc_curve_values.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vpac_dual_path_roc_curve.pdf"), format="pdf")
plt.close()

threshold_records = []
best_threshold = 0.5
best_f1 = -1.0

for threshold in np.arange(0.1, 0.9, 0.01):
    pred = (val_probs >= threshold).astype(int)

    accuracy = accuracy_score(val_labels, pred)
    precision = precision_score(val_labels, pred, zero_division=0)
    recall = recall_score(val_labels, pred, zero_division=0)
    f1 = f1_score(val_labels, pred, zero_division=0)
    mcc = matthews_corrcoef(val_labels, pred)

    threshold_records.append({
        "threshold": round(float(threshold), 2),
        "accuracy": float(accuracy),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "mcc": float(mcc)
    })

    if f1 > best_f1:
        best_f1 = f1
        best_threshold = threshold

threshold_df = pd.DataFrame(threshold_records)
threshold_df.to_csv(os.path.join(outdir, "vpac_dual_path_threshold_search.csv"), index=False)

val_preds_best = (val_probs >= best_threshold).astype(int)

best_threshold_metrics = {
    "best_threshold": float(best_threshold),
    "accuracy": float(accuracy_score(val_labels, val_preds_best)),
    "precision": float(precision_score(val_labels, val_preds_best, zero_division=0)),
    "recall": float(recall_score(val_labels, val_preds_best, zero_division=0)),
    "f1": float(f1_score(val_labels, val_preds_best, zero_division=0)),
    "mcc": float(matthews_corrcoef(val_labels, val_preds_best))
}

pd.DataFrame([best_threshold_metrics]).to_csv(os.path.join(outdir, "vpac_dual_path_best_threshold_metrics.csv"), index=False)

with open(os.path.join(outdir, "vpac_dual_path_best_threshold_metrics.json"), "w") as f:
    json.dump(best_threshold_metrics, f, indent=2)

cm = confusion_matrix(val_labels, val_preds, labels=[0, 1])
TN, FP, FN, TP = cm.ravel()
total = cm.sum()

confusion_metrics = {
    "TN": int(TN),
    "FP": int(FP),
    "FN": int(FN),
    "TP": int(TP),
    "TN_ratio": float(TN / total),
    "FP_ratio": float(FP / total),
    "FN_ratio": float(FN / total),
    "TP_ratio": float(TP / total)
}

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "vpac_dual_path_confusion_matrix_values.csv"), index=False)

fig, ax = plt.subplots(figsize=(6, 5.5))
sns.heatmap(
    cm,
    annot=False,
    fmt="d",
    cmap="Blues",
    xticklabels=["Negative", "Positive"],
    yticklabels=["Negative", "Positive"],
    cbar=False
)

for i in range(2):
    for j in range(2):
        value = cm[i, j]
        ratio = value / total
        font_color = "white" if value > cm.max() / 2 else "black"

        if i == 0 and j == 0:
            label = f"TN: {value}\n({ratio:.2%})"
        elif i == 0 and j == 1:
            label = f"FP: {value}\n({ratio:.2%})"
        elif i == 1 and j == 0:
            label = f"FN: {value}\n({ratio:.2%})"
        else:
            label = f"TP: {value}\n({ratio:.2%})"

        ax.text(
            j + 0.5,
            i + 0.5,
            label,
            ha="center",
            va="center",
            color=font_color,
            fontsize=12
        )

plt.title("Confusion Matrix")
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vpac_dual_path_confusion_matrix.pdf"), format="pdf")
plt.close()

print("Final validation metrics:")
print(final_metrics)
print("Best threshold metrics:")
print(best_threshold_metrics)
print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")