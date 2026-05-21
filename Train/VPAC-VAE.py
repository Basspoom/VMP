import os
import json
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
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

epochs = 200
batch_size = 1145
learning_rate = 0.001
latent_dim = 4
dropout_prob = 0.5
beta = 1.0
input_dim = 72

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")

X_train_raw = train_df.iloc[:, 1:25].values
X_val_raw = vali_df.iloc[:, 1:25].values

y_train = train_df.iloc[:, 25].map({"Virus": 1, "Non-Virus": 0}).values.astype(np.float32)
y_val = vali_df.iloc[:, 25].map({"Virus": 1, "Non-Virus": 0}).values.astype(np.float32)

def manual_one_hot_encode(x):
    x = np.asarray(x)
    encoded = np.zeros((x.shape[0], x.shape[1] * 3), dtype=np.float32)
    x_round = np.rint(x).astype(int)

    for i in range(x_round.shape[0]):
        for j in range(x_round.shape[1]):
            value = x_round[i, j]
            base = j * 3
            if value == -1:
                encoded[i, base] = 1.0
            elif value == 0:
                encoded[i, base + 1] = 1.0
            elif value == 1:
                encoded[i, base + 2] = 1.0
            else:
                raise ValueError(f"Unexpected score value {value} at row {i}, feature {j}")

    return encoded

X_train = manual_one_hot_encode(X_train_raw)
X_val = manual_one_hot_encode(X_val_raw)

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

class VAE(nn.Module):
    def __init__(self, input_dim=72, latent_dim=4, dropout_prob=0.5):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, 128)
        self.dropout1 = nn.Dropout(dropout_prob)
        self.fc2 = nn.Linear(128, 512)
        self.fc_mu = nn.Linear(512, latent_dim)
        self.fc_logvar = nn.Linear(512, latent_dim)

        self.fc3 = nn.Linear(latent_dim, 512)
        self.fc4 = nn.Linear(512, 128)
        self.fc5 = nn.Linear(128, input_dim)

        self.fc_class = nn.Linear(latent_dim, 1)

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        h1 = self.dropout1(h1)
        h2 = F.relu(self.fc2(h1))
        mu = self.fc_mu(h2)
        logvar = self.fc_logvar(h2)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        h3 = F.relu(self.fc3(z))
        h4 = F.relu(self.fc4(h3))
        recon_x = torch.sigmoid(self.fc5(h4))
        return recon_x

    def classify(self, z):
        return torch.sigmoid(self.fc_class(z))

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        class_pred = self.classify(z)
        return recon_x, class_pred, mu, logvar, z

def loss_function(recon_x, x, class_pred, y, mu, logvar, beta=1.0):
    recon_loss = F.binary_cross_entropy(recon_x, x, reduction="sum")
    kld = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    class_loss = F.binary_cross_entropy(class_pred, y, reduction="sum")
    total_loss = recon_loss + beta * kld + class_loss
    return total_loss, recon_loss, kld, class_loss

model = VAE(input_dim=input_dim, latent_dim=latent_dim, dropout_prob=dropout_prob).to(device)
optimizer = optim.Adam(model.parameters(), lr=learning_rate)

model_config = {
    "input_dim": input_dim,
    "latent_dim": latent_dim,
    "encoder_widths": [72, 128, 512],
    "decoder_widths": [4, 512, 128, 72],
    "dropout_prob": dropout_prob,
    "learning_rate": learning_rate,
    "batch_size": batch_size,
    "epochs": epochs,
    "beta": beta
}

with open(os.path.join(outdir, "vae_model_config.json"), "w") as f:
    json.dump(model_config, f, indent=2)

metrics_history = {
    "epoch": [],
    "train_loss": [],
    "val_loss": [],
    "train_recon_loss": [],
    "val_recon_loss": [],
    "train_kld": [],
    "val_kld": [],
    "train_class_loss": [],
    "val_class_loss": [],
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
    train_recon_sum = 0.0
    train_kld_sum = 0.0
    train_class_sum = 0.0
    train_probs = []
    train_labels = []

    for batch_x, batch_y in train_loader:
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)

        optimizer.zero_grad()
        recon_x, class_pred, mu, logvar, z = model(batch_x)
        loss, recon_loss, kld, class_loss = loss_function(
            recon_x,
            batch_x,
            class_pred,
            batch_y,
            mu,
            logvar,
            beta=beta
        )
        loss.backward()
        optimizer.step()

        train_loss_sum += loss.item()
        train_recon_sum += recon_loss.item()
        train_kld_sum += kld.item()
        train_class_sum += class_loss.item()

        train_probs.append(class_pred.detach().cpu().numpy())
        train_labels.append(batch_y.detach().cpu().numpy())

    train_probs = np.concatenate(train_probs, axis=0).reshape(-1)
    train_labels = np.concatenate(train_labels, axis=0).reshape(-1).astype(int)
    train_preds = (train_probs >= 0.5).astype(int)

    train_loss = train_loss_sum / len(train_dataset)
    train_recon_loss = train_recon_sum / len(train_dataset)
    train_kld = train_kld_sum / len(train_dataset)
    train_class_loss = train_class_sum / len(train_dataset)

    train_accuracy = accuracy_score(train_labels, train_preds)
    train_precision = precision_score(train_labels, train_preds, zero_division=0)
    train_recall = recall_score(train_labels, train_preds, zero_division=0)
    train_f1 = f1_score(train_labels, train_preds, zero_division=0)
    train_mcc = matthews_corrcoef(train_labels, train_preds)

    model.eval()

    val_loss_sum = 0.0
    val_recon_sum = 0.0
    val_kld_sum = 0.0
    val_class_sum = 0.0
    val_probs = []
    val_labels = []

    with torch.no_grad():
        for batch_x, batch_y in val_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)

            recon_x, class_pred, mu, logvar, z = model(batch_x)
            loss, recon_loss, kld, class_loss = loss_function(
                recon_x,
                batch_x,
                class_pred,
                batch_y,
                mu,
                logvar,
                beta=beta
            )

            val_loss_sum += loss.item()
            val_recon_sum += recon_loss.item()
            val_kld_sum += kld.item()
            val_class_sum += class_loss.item()

            val_probs.append(class_pred.detach().cpu().numpy())
            val_labels.append(batch_y.detach().cpu().numpy())

    val_probs = np.concatenate(val_probs, axis=0).reshape(-1)
    val_labels = np.concatenate(val_labels, axis=0).reshape(-1).astype(int)
    val_preds = (val_probs >= 0.5).astype(int)

    val_loss = val_loss_sum / len(val_dataset)
    val_recon_loss = val_recon_sum / len(val_dataset)
    val_kld = val_kld_sum / len(val_dataset)
    val_class_loss = val_class_sum / len(val_dataset)

    val_accuracy = accuracy_score(val_labels, val_preds)
    val_precision = precision_score(val_labels, val_preds, zero_division=0)
    val_recall = recall_score(val_labels, val_preds, zero_division=0)
    val_f1 = f1_score(val_labels, val_preds, zero_division=0)
    val_mcc = matthews_corrcoef(val_labels, val_preds)

    metrics_history["epoch"].append(epoch + 1)
    metrics_history["train_loss"].append(train_loss)
    metrics_history["val_loss"].append(val_loss)
    metrics_history["train_recon_loss"].append(train_recon_loss)
    metrics_history["val_recon_loss"].append(val_recon_loss)
    metrics_history["train_kld"].append(train_kld)
    metrics_history["val_kld"].append(val_kld)
    metrics_history["train_class_loss"].append(train_class_loss)
    metrics_history["val_class_loss"].append(val_class_loss)
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
        torch.save(model.state_dict(), os.path.join(outdir, "vae_best_loss_model.pth"))

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        torch.save(model.state_dict(), os.path.join(outdir, "vae_best_f1_model.pth"))

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}, "
        f"Train Recon: {train_recon_loss:.6f}, Val Recon: {val_recon_loss:.6f}, "
        f"Train KLD: {train_kld:.6f}, Val KLD: {val_kld:.6f}, "
        f"Train Class Loss: {train_class_loss:.6f}, Val Class Loss: {val_class_loss:.6f}, "
        f"Train Acc: {train_accuracy:.6f}, Val Acc: {val_accuracy:.6f}, "
        f"Train Precision: {train_precision:.6f}, Val Precision: {val_precision:.6f}, "
        f"Train Recall: {train_recall:.6f}, Val Recall: {val_recall:.6f}, "
        f"Train F1: {train_f1:.6f}, Val F1: {val_f1:.6f}, "
        f"Train MCC: {train_mcc:.6f}, Val MCC: {val_mcc:.6f}"
    )

torch.save(model.state_dict(), os.path.join(outdir, "vae_final_model.pth"))

metrics_df = pd.DataFrame(metrics_history)
metrics_df.to_csv(os.path.join(outdir, "vae_training_metrics.csv"), index=False)

model.eval()

train_probs = []
train_labels = []
train_mu_list = []
train_logvar_list = []
train_z_list = []

with torch.no_grad():
    for batch_x, batch_y in train_loader:
        batch_x = batch_x.to(device)
        recon_x, class_pred, mu, logvar, z = model(batch_x)
        train_probs.append(class_pred.detach().cpu().numpy())
        train_labels.append(batch_y.detach().cpu().numpy())
        train_mu_list.append(mu.detach().cpu().numpy())
        train_logvar_list.append(logvar.detach().cpu().numpy())
        train_z_list.append(z.detach().cpu().numpy())

train_probs = np.concatenate(train_probs, axis=0).reshape(-1)
train_labels = np.concatenate(train_labels, axis=0).reshape(-1).astype(int)
train_preds = (train_probs >= 0.5).astype(int)
train_mu = np.concatenate(train_mu_list, axis=0)
train_logvar = np.concatenate(train_logvar_list, axis=0)
train_z = np.concatenate(train_z_list, axis=0)

val_probs = []
val_labels = []
val_mu_list = []
val_logvar_list = []
val_z_list = []

with torch.no_grad():
    for batch_x, batch_y in val_loader:
        batch_x = batch_x.to(device)
        recon_x, class_pred, mu, logvar, z = model(batch_x)
        val_probs.append(class_pred.detach().cpu().numpy())
        val_labels.append(batch_y.detach().cpu().numpy())
        val_mu_list.append(mu.detach().cpu().numpy())
        val_logvar_list.append(logvar.detach().cpu().numpy())
        val_z_list.append(z.detach().cpu().numpy())

val_probs = np.concatenate(val_probs, axis=0).reshape(-1)
val_labels = np.concatenate(val_labels, axis=0).reshape(-1).astype(int)
val_preds = (val_probs >= 0.5).astype(int)
val_mu = np.concatenate(val_mu_list, axis=0)
val_logvar = np.concatenate(val_logvar_list, axis=0)
val_z = np.concatenate(val_z_list, axis=0)

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

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "vae_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "vae_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

latent_cols = [f"z_{i + 1}" for i in range(latent_dim)]
mu_cols = [f"mu_{i + 1}" for i in range(latent_dim)]
logvar_cols = [f"logvar_{i + 1}" for i in range(latent_dim)]

train_latent_df = pd.concat(
    [
        train_df.reset_index(drop=True),
        pd.DataFrame(train_z, columns=latent_cols),
        pd.DataFrame(train_mu, columns=mu_cols),
        pd.DataFrame(train_logvar, columns=logvar_cols)
    ],
    axis=1
)
train_latent_df.to_csv(os.path.join(outdir, "vae_train_latent.csv"), index=False)

vali_pred_df = pd.concat(
    [
        vali_df.reset_index(drop=True),
        pd.DataFrame(val_z, columns=latent_cols),
        pd.DataFrame(val_mu, columns=mu_cols),
        pd.DataFrame(val_logvar, columns=logvar_cols)
    ],
    axis=1
)
vali_pred_df["true_label"] = val_labels
vali_pred_df["virus_probability"] = val_probs
vali_pred_df["pred_label"] = val_preds
vali_pred_df["pred_class"] = np.where(val_preds == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "vae_validation_predictions.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(metrics_history["epoch"], metrics_history["train_loss"], label="Train Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_loss"], label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("VAE Total Loss")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vae_loss_plot.pdf"), format="pdf")
plt.close()

plt.figure(figsize=(8, 6))
plt.plot(metrics_history["epoch"], metrics_history["train_recon_loss"], label="Train Reconstruction Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_recon_loss"], label="Validation Reconstruction Loss")
plt.plot(metrics_history["epoch"], metrics_history["train_kld"], label="Train KL Divergence")
plt.plot(metrics_history["epoch"], metrics_history["val_kld"], label="Validation KL Divergence")
plt.plot(metrics_history["epoch"], metrics_history["train_class_loss"], label="Train Classification Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_class_loss"], label="Validation Classification Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("VAE Loss Components")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vae_loss_components.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(val_labels, val_probs)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "vae_roc_curve_values.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "vae_roc_curve.pdf"), format="pdf")
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
threshold_df.to_csv(os.path.join(outdir, "vae_threshold_search.csv"), index=False)

val_preds_best = (val_probs >= best_threshold).astype(int)

best_threshold_metrics = {
    "best_threshold": float(best_threshold),
    "accuracy": float(accuracy_score(val_labels, val_preds_best)),
    "precision": float(precision_score(val_labels, val_preds_best, zero_division=0)),
    "recall": float(recall_score(val_labels, val_preds_best, zero_division=0)),
    "f1": float(f1_score(val_labels, val_preds_best, zero_division=0)),
    "mcc": float(matthews_corrcoef(val_labels, val_preds_best))
}

pd.DataFrame([best_threshold_metrics]).to_csv(os.path.join(outdir, "vae_best_threshold_metrics.csv"), index=False)

with open(os.path.join(outdir, "vae_best_threshold_metrics.json"), "w") as f:
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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "vae_confusion_matrix_values.csv"), index=False)

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
plt.savefig(os.path.join(outdir, "vae_confusion_matrix.pdf"), format="pdf")
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