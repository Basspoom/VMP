import os
import json
import joblib
import random
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader, Dataset
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_curve,
    auc,
    confusion_matrix,
    matthews_corrcoef
)
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

seed =777
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

batch_size = 32
epochs = 200
learning_rate = 0.001
dropout_rate = 0.3

train_data = pd.read_csv(train_path, sep="\t")
vali_data = pd.read_csv(vali_path, sep="\t")

X_train = train_data.iloc[:, 1:25].values.astype(np.float32)
X_val = vali_data.iloc[:, 1:25].values.astype(np.float32)

def encode_labels(y):
    if y.dtype == object:
        return y.map({"Virus": 1, "Non-Virus": 0}).values.astype(np.int64)
    return y.values.astype(np.int64)

y_train = encode_labels(train_data.iloc[:, 25])
y_val = encode_labels(vali_data.iloc[:, 25])

scaler = StandardScaler()
X_train = scaler.fit_transform(X_train).astype(np.float32)
X_val = scaler.transform(X_val).astype(np.float32)

joblib.dump(scaler, os.path.join(outdir, "cnn_standard_scaler.pkl"))

X_train = X_train.reshape(-1, 1, 24)
X_val = X_val.reshape(-1, 1, 24)

class CustomDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class CNNClassifier(nn.Module):
    def __init__(self, dropout_rate=0.3):
        super().__init__()
        self.conv_layers = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(64, 128, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2)
        )
        self.fc_layers = nn.Sequential(
            nn.Linear(128 * 3, 256),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(128, 2)
        )

    def forward(self, x):
        x = self.conv_layers(x)
        x = torch.flatten(x, start_dim=1)
        x = self.fc_layers(x)
        return x

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

train_dataset = CustomDataset(X_train, y_train)
val_dataset = CustomDataset(X_val, y_val)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

model = CNNClassifier(dropout_rate=dropout_rate).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = Adam(model.parameters(), lr=learning_rate)

metrics_path = os.path.join(outdir, "training_metrics_cnn.csv")
with open(metrics_path, "w") as f:
    f.write(
        "epoch,train_loss,val_loss,"
        "train_acc,val_acc,"
        "train_precision,val_precision,"
        "train_recall,val_recall,"
        "train_f1,val_f1,"
        "train_mcc,val_mcc\n"
    )

train_losses = []
val_losses = []
best_val_loss = float("inf")
best_state_dict = None

for epoch in range(epochs):
    model.train()
    train_loss_sum = 0.0
    train_preds = []
    train_labels = []

    for features, labels in tqdm(train_loader, desc=f"Epoch {epoch + 1}/{epochs}", ncols=100, unit="batch"):
        features = features.to(device)
        labels = labels.to(device)

        optimizer.zero_grad()
        logits = model(features)
        loss = criterion(logits, labels)
        loss.backward()
        optimizer.step()

        train_loss_sum += loss.item()
        preds = torch.argmax(logits, dim=1)

        train_preds.extend(preds.detach().cpu().numpy())
        train_labels.extend(labels.detach().cpu().numpy())

    train_loss = train_loss_sum / len(train_loader)
    train_losses.append(train_loss)

    model.eval()
    val_loss_sum = 0.0
    val_preds = []
    val_labels = []

    with torch.no_grad():
        for features, labels in tqdm(val_loader, desc=f"Validation {epoch + 1}/{epochs}", ncols=100, unit="batch"):
            features = features.to(device)
            labels = labels.to(device)

            logits = model(features)
            loss = criterion(logits, labels)

            val_loss_sum += loss.item()
            preds = torch.argmax(logits, dim=1)

            val_preds.extend(preds.detach().cpu().numpy())
            val_labels.extend(labels.detach().cpu().numpy())

    val_loss = val_loss_sum / len(val_loader)
    val_losses.append(val_loss)

    train_acc = accuracy_score(train_labels, train_preds)
    val_acc = accuracy_score(val_labels, val_preds)
    train_precision = precision_score(train_labels, train_preds, zero_division=0)
    val_precision = precision_score(val_labels, val_preds, zero_division=0)
    train_recall = recall_score(train_labels, train_preds, zero_division=0)
    val_recall = recall_score(val_labels, val_preds, zero_division=0)
    train_f1 = f1_score(train_labels, train_preds, zero_division=0)
    val_f1 = f1_score(val_labels, val_preds, zero_division=0)
    train_mcc = matthews_corrcoef(train_labels, train_preds)
    val_mcc = matthews_corrcoef(val_labels, val_preds)

    with open(metrics_path, "a") as f:
        f.write(
            f"{epoch + 1},"
            f"{train_loss:.8f},{val_loss:.8f},"
            f"{train_acc:.8f},{val_acc:.8f},"
            f"{train_precision:.8f},{val_precision:.8f},"
            f"{train_recall:.8f},{val_recall:.8f},"
            f"{train_f1:.8f},{val_f1:.8f},"
            f"{train_mcc:.8f},{val_mcc:.8f}\n"
        )

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, "
        f"Train Acc: {train_acc:.4f}, Val Acc: {val_acc:.4f}, "
        f"Train Precision: {train_precision:.4f}, Val Precision: {val_precision:.4f}, "
        f"Train Recall: {train_recall:.4f}, Val Recall: {val_recall:.4f}, "
        f"Train F1: {train_f1:.4f}, Val F1: {val_f1:.4f}, "
        f"Train MCC: {train_mcc:.4f}, Val MCC: {val_mcc:.4f}"
    )

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        best_state_dict = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        torch.save(best_state_dict, os.path.join(outdir, "best_cnn_model.pth"))

model.load_state_dict(best_state_dict)
model.to(device)
model.eval()

plt.figure()
plt.plot(range(1, epochs + 1), train_losses, label="Train Loss")
plt.plot(range(1, epochs + 1), val_losses, label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Loss Per Epoch")
plt.legend(loc="upper right")
plt.savefig(os.path.join(outdir, "cnn_loss_per_epoch.pdf"), format="pdf")
plt.close()

all_probs = []
all_preds_default = []
all_labels_final = []

with torch.no_grad():
    for features, labels in val_loader:
        features = features.to(device)
        logits = model(features)
        probs = torch.softmax(logits, dim=1)[:, 1]
        preds = torch.argmax(logits, dim=1)

        all_probs.extend(probs.detach().cpu().numpy())
        all_preds_default.extend(preds.detach().cpu().numpy())
        all_labels_final.extend(labels.numpy())

all_probs = np.array(all_probs)
all_labels_final = np.array(all_labels_final)

threshold_records = []
best_threshold = 0.5
best_precision = -1.0

for threshold in np.arange(0.1, 0.9, 0.01):
    y_pred_threshold = (all_probs >= threshold).astype(int)
    precision = precision_score(all_labels_final, y_pred_threshold, zero_division=0)
    recall = recall_score(all_labels_final, y_pred_threshold, zero_division=0)
    f1 = f1_score(all_labels_final, y_pred_threshold, zero_division=0)
    mcc = matthews_corrcoef(all_labels_final, y_pred_threshold)

    threshold_records.append({
        "threshold": round(float(threshold), 2),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "mcc": mcc
    })

    if precision > best_precision:
        best_precision = precision
        best_threshold = threshold

threshold_df = pd.DataFrame(threshold_records)
threshold_df.to_csv(os.path.join(outdir, "threshold_search.csv"), index=False)

y_pred = (all_probs >= best_threshold).astype(int)

cm = confusion_matrix(all_labels_final, y_pred)
TN, FP, FN, TP = cm.ravel()

accuracy = accuracy_score(all_labels_final, y_pred)
precision = precision_score(all_labels_final, y_pred, zero_division=0)
recall = recall_score(all_labels_final, y_pred, zero_division=0)
f1 = f1_score(all_labels_final, y_pred, zero_division=0)
mcc = matthews_corrcoef(all_labels_final, y_pred)

final_metrics = {
    "best_threshold": float(best_threshold),
    "accuracy": float(accuracy),
    "precision": float(precision),
    "recall": float(recall),
    "f1": float(f1),
    "mcc": float(mcc),
    "TP": int(TP),
    "TN": int(TN),
    "FP": int(FP),
    "FN": int(FN)
}

with open(os.path.join(outdir, "final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "final_metrics.csv"), index=False)

vali_pred_df = vali_data.copy()
vali_pred_df["true_label"] = all_labels_final
vali_pred_df["virus_probability"] = all_probs
vali_pred_df["pred_label"] = y_pred
vali_pred_df["pred_class"] = np.where(y_pred == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "validation_predictions.csv"), index=False)

print(f"Best Threshold: {best_threshold:.2f}, Best Precision: {best_precision:.4f}")
print(f"Optimized Accuracy: {accuracy:.4f}")
print(f"Optimized Precision: {precision:.4f}")
print(f"Optimized Recall: {recall:.4f}")
print(f"Optimized F1: {f1:.4f}")
print(f"Optimized MCC: {mcc:.4f}")
print(f"Confusion Matrix:\n{cm}")

total = cm.sum()
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
        elif i == 1 and j == 1:
            label = f"TP: {value}\n({ratio:.2%})"
        elif i == 0 and j == 1:
            label = f"FP: {value}\n({ratio:.2%})"
        else:
            label = f"FN: {value}\n({ratio:.2%})"

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
plt.savefig(os.path.join(outdir, "cnn_confusion_matrix.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(all_labels_final, all_probs)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.savefig(os.path.join(outdir, "cnn_roc_curve.pdf"), format="pdf")
plt.close()

print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")