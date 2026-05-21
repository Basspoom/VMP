import os
import json
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
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

seed =777
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

epochs = 100
batch_size = 888
learning_rate = 0.001
weight_decay = 1e-5
dropout_rate = 0.3

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")

X_train = train_df.iloc[:, 1:25].values.astype(np.float32)
y_train = train_df.iloc[:, 25].map({"Virus": 1, "Non-Virus": 0}).values.astype(np.float32)

X_val = vali_df.iloc[:, 1:25].values.astype(np.float32)
y_val = vali_df.iloc[:, 25].map({"Virus": 1, "Non-Virus": 0}).values.astype(np.float32)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train, dtype=torch.float32).view(-1, 1)

X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
y_val_tensor = torch.tensor(y_val, dtype=torch.float32).view(-1, 1)

train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
val_dataset = TensorDataset(X_val_tensor, y_val_tensor)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=False)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, drop_last=False)

class FeatureExpansion(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc = nn.Linear(1, 128)

    def forward(self, x):
        return self.fc(x)

class Encoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc = nn.Linear(128, 512)

    def forward(self, x):
        return self.fc(x)

class Decoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc = nn.Linear(512, 128)

    def forward(self, x):
        return self.fc(x)

class MLP(nn.Module):
    def __init__(self, dropout_rate=0.3):
        super().__init__()
        self.fc1 = nn.Linear(3072, 64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout_rate)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.sigmoid(self.fc3(x))
        return x

class VirusPredictionModel(nn.Module):
    def __init__(self, dropout_rate=0.3):
        super().__init__()
        self.feature_expansion = FeatureExpansion()
        self.encoder = Encoder()
        self.decoder = Decoder()
        self.mlp = MLP(dropout_rate=dropout_rate)

    def forward(self, x):
        x = x.view(-1, 24, 1)
        x = self.feature_expansion(x)
        x = self.encoder(x)
        x = self.decoder(x)
        x = x.view(-1, 24 * 128)
        x = self.mlp(x)
        return x

model = VirusPredictionModel(dropout_rate=dropout_rate).to(device)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

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

best_val_loss = float("inf")
best_val_f1 = -1.0

for epoch in tqdm(range(epochs), desc="Training progress"):
    model.train()
    train_loss_sum = 0.0
    train_probs = []
    train_labels = []

    for batch_x, batch_y in train_loader:
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)

        optimizer.zero_grad()
        output = model(batch_x)
        loss = criterion(output, batch_y)
        loss.backward()
        optimizer.step()

        train_loss_sum += loss.item()
        train_probs.append(output.detach().cpu().numpy())
        train_labels.append(batch_y.detach().cpu().numpy())

    train_loss = train_loss_sum / len(train_loader)
    train_probs = np.concatenate(train_probs, axis=0).reshape(-1)
    train_labels = np.concatenate(train_labels, axis=0).reshape(-1).astype(int)
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
        for batch_x, batch_y in val_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)

            output = model(batch_x)
            loss = criterion(output, batch_y)

            val_loss_sum += loss.item()
            val_probs.append(output.detach().cpu().numpy())
            val_labels.append(batch_y.detach().cpu().numpy())

    val_loss = val_loss_sum / len(val_loader)
    val_probs = np.concatenate(val_probs, axis=0).reshape(-1)
    val_labels = np.concatenate(val_labels, axis=0).reshape(-1).astype(int)
    val_preds = (val_probs >= 0.5).astype(int)

    val_accuracy = accuracy_score(val_labels, val_preds)
    val_precision = precision_score(val_labels, val_preds, zero_division=0)
    val_recall = recall_score(val_labels, val_preds, zero_division=0)
    val_f1 = f1_score(val_labels, val_preds, zero_division=0)
    val_mcc = matthews_corrcoef(val_labels, val_preds)

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

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), os.path.join(outdir, "ae_best_loss_model.pth"))

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        torch.save(model.state_dict(), os.path.join(outdir, "ae_best_f1_model.pth"))

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}, "
        f"Train Acc: {train_accuracy:.6f}, Val Acc: {val_accuracy:.6f}, "
        f"Train Precision: {train_precision:.6f}, Val Precision: {val_precision:.6f}, "
        f"Train Recall: {train_recall:.6f}, Val Recall: {val_recall:.6f}, "
        f"Train F1: {train_f1:.6f}, Val F1: {val_f1:.6f}, "
        f"Train MCC: {train_mcc:.6f}, Val MCC: {val_mcc:.6f}"
    )

torch.save(model.state_dict(), os.path.join(outdir, "ae_final_model.pth"))

metrics_df = pd.DataFrame(metrics_history)
metrics_df.to_csv(os.path.join(outdir, "ae_training_metrics.csv"), index=False)

model.eval()
val_probs = []
val_labels = []

with torch.no_grad():
    for batch_x, batch_y in val_loader:
        batch_x = batch_x.to(device)
        output = model(batch_x)
        val_probs.append(output.detach().cpu().numpy())
        val_labels.append(batch_y.detach().cpu().numpy())

val_probs = np.concatenate(val_probs, axis=0).reshape(-1)
val_labels = np.concatenate(val_labels, axis=0).reshape(-1).astype(int)
val_preds = (val_probs >= 0.5).astype(int)

train_probs = []
train_labels = []

with torch.no_grad():
    for batch_x, batch_y in train_loader:
        batch_x = batch_x.to(device)
        output = model(batch_x)
        train_probs.append(output.detach().cpu().numpy())
        train_labels.append(batch_y.detach().cpu().numpy())

train_probs = np.concatenate(train_probs, axis=0).reshape(-1)
train_labels = np.concatenate(train_labels, axis=0).reshape(-1).astype(int)
train_preds = (train_probs >= 0.5).astype(int)

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

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "ae_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "ae_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_df.copy()
vali_pred_df["true_label"] = val_labels
vali_pred_df["virus_probability"] = val_probs
vali_pred_df["pred_label"] = val_preds
vali_pred_df["pred_class"] = np.where(val_preds == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "ae_validation_predictions.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(metrics_history["epoch"], metrics_history["train_loss"], label="Training Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_loss"], label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Training and Validation Loss")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "ae_loss_curve.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(val_labels, val_probs)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "ae_roc_curve_values.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "ae_roc_curve.pdf"), format="pdf")
plt.close()

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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "ae_confusion_matrix_values.csv"), index=False)

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
plt.savefig(os.path.join(outdir, "ae_confusion_matrix.pdf"), format="pdf")
plt.close()

print("Final validation metrics:")
print(final_metrics)
print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")