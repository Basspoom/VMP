import os
import json
import random
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import seaborn as sns

from kan import *
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

seed =777
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

epochs = 20
batch_size = 6000
learning_rate = 0.1

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")

X_train = train_df.iloc[:, 1:25].values.astype(np.float32)
y_train = (train_df.iloc[:, 25].values == "Virus").astype(np.float32)

X_val = vali_df.iloc[:, 1:25].values.astype(np.float32)
y_val = (vali_df.iloc[:, 25].values == "Virus").astype(np.float32)

device_str = "cuda" if torch.cuda.is_available() else "cpu"
device = torch.device(device_str)

X_train = torch.tensor(X_train, dtype=torch.float32).to(device)
y_train = torch.tensor(y_train, dtype=torch.float32).view(-1, 1).to(device)

X_val = torch.tensor(X_val, dtype=torch.float32).to(device)
y_val = torch.tensor(y_val, dtype=torch.float32).view(-1, 1).to(device)

dataset = {
    "train_input": X_train,
    "train_label": y_train,
    "test_input": X_val,
    "test_label": y_val
}

model = KAN(
    width=[24, 49, 1],
    grid=5,
    k=3,
    seed=seed,
    device=device_str
)

loss_fn = torch.nn.BCEWithLogitsLoss()

metrics_history = {
    "epoch": [],
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

best_val_f1 = -1.0
best_val_precision = -1.0
best_val_recall = -1.0

for epoch in tqdm(range(epochs), desc="Training Epochs"):
    results = model.fit(
        dataset,
        opt="LBFGS",
        steps=1,
        lr=learning_rate,
        batch=batch_size,
        loss_fn=loss_fn
    )

    with torch.no_grad():
        train_logits = model(dataset["train_input"])
        val_logits = model(dataset["test_input"])

        train_loss = loss_fn(train_logits, dataset["train_label"]).item()
        val_loss = loss_fn(val_logits, dataset["test_label"]).item()

        train_prob = torch.sigmoid(train_logits).detach().cpu().numpy().reshape(-1)
        val_prob = torch.sigmoid(val_logits).detach().cpu().numpy().reshape(-1)

        train_pred = (train_prob >= 0.5).astype(int)
        val_pred = (val_prob >= 0.5).astype(int)

        train_true = dataset["train_label"].detach().cpu().numpy().reshape(-1).astype(int)
        val_true = dataset["test_label"].detach().cpu().numpy().reshape(-1).astype(int)

    train_accuracy = accuracy_score(train_true, train_pred)
    val_accuracy = accuracy_score(val_true, val_pred)

    train_precision = precision_score(train_true, train_pred, zero_division=0)
    val_precision = precision_score(val_true, val_pred, zero_division=0)

    train_recall = recall_score(train_true, train_pred, zero_division=0)
    val_recall = recall_score(val_true, val_pred, zero_division=0)

    train_f1 = f1_score(train_true, train_pred, zero_division=0)
    val_f1 = f1_score(val_true, val_pred, zero_division=0)

    train_mcc = matthews_corrcoef(train_true, train_pred)
    val_mcc = matthews_corrcoef(val_true, val_pred)

    metrics_history["epoch"].append(epoch + 1)
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

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}, "
        f"Train Acc: {train_accuracy:.6f}, Val Acc: {val_accuracy:.6f}, "
        f"Train Precision: {train_precision:.6f}, Val Precision: {val_precision:.6f}, "
        f"Train Recall: {train_recall:.6f}, Val Recall: {val_recall:.6f}, "
        f"Train F1: {train_f1:.6f}, Val F1: {val_f1:.6f}, "
        f"Train MCC: {train_mcc:.6f}, Val MCC: {val_mcc:.6f}"
    )

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        torch.save(model.state_dict(), os.path.join(outdir, "best_kan_f1_model.pt"))

    if val_precision > best_val_precision:
        best_val_precision = val_precision
        torch.save(model.state_dict(), os.path.join(outdir, "best_kan_precision_model.pt"))

    if val_recall > best_val_recall:
        best_val_recall = val_recall
        torch.save(model.state_dict(), os.path.join(outdir, "best_kan_recall_model.pt"))

torch.save(model.state_dict(), os.path.join(outdir, "final_kan_model.pt"))

metrics_df = pd.DataFrame(metrics_history)
metrics_df.to_csv(os.path.join(outdir, "kan_training_metrics.csv"), index=False)

with torch.no_grad():
    train_logits = model(dataset["train_input"])
    val_logits = model(dataset["test_input"])

    train_prob = torch.sigmoid(train_logits).detach().cpu().numpy().reshape(-1)
    val_prob = torch.sigmoid(val_logits).detach().cpu().numpy().reshape(-1)

train_true = dataset["train_label"].detach().cpu().numpy().reshape(-1).astype(int)
val_true = dataset["test_label"].detach().cpu().numpy().reshape(-1).astype(int)

train_pred = (train_prob >= 0.5).astype(int)
val_pred = (val_prob >= 0.5).astype(int)

final_metrics = {
    "train_accuracy": float(accuracy_score(train_true, train_pred)),
    "val_accuracy": float(accuracy_score(val_true, val_pred)),
    "train_precision": float(precision_score(train_true, train_pred, zero_division=0)),
    "val_precision": float(precision_score(val_true, val_pred, zero_division=0)),
    "train_recall": float(recall_score(train_true, train_pred, zero_division=0)),
    "val_recall": float(recall_score(val_true, val_pred, zero_division=0)),
    "train_f1": float(f1_score(train_true, train_pred, zero_division=0)),
    "val_f1": float(f1_score(val_true, val_pred, zero_division=0)),
    "train_mcc": float(matthews_corrcoef(train_true, train_pred)),
    "val_mcc": float(matthews_corrcoef(val_true, val_pred))
}

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "kan_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "kan_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_df.copy()
vali_pred_df["true_label"] = val_true
vali_pred_df["virus_probability"] = val_prob
vali_pred_df["pred_label"] = val_pred
vali_pred_df["pred_class"] = np.where(val_pred == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "kan_validation_predictions.csv"), index=False)

plt.figure()
plt.plot(metrics_history["epoch"], metrics_history["train_loss"], label="Train Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_loss"], label="Validation Loss")
plt.title("Loss vs Epoch")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.savefig(os.path.join(outdir, "kan_loss_vs_epoch.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(val_true, val_prob)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "kan_roc_curve_values.csv"), index=False)

plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.savefig(os.path.join(outdir, "kan_roc_curve.pdf"), format="pdf")
plt.close()

cm = confusion_matrix(val_true, val_pred)
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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "kan_confusion_matrix_values.csv"), index=False)

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
plt.savefig(os.path.join(outdir, "kan_confusion_matrix.pdf"), format="pdf")
plt.close()

print("Final validation metrics:")
print(final_metrics)
print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")