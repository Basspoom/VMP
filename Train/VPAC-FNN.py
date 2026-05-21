import os
import json
import random
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
import seaborn as sns

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Input, Dense, Dropout
from tensorflow.keras import regularizers
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef,
    roc_curve,
    auc,
    confusion_matrix,
    log_loss
)
from tqdm import tqdm

seed =777
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

epochs = 200
batch_size = 32
learning_rate = 0.1
dropout_rate = 0.3
l2_rate = 0.01
hidden_units = [48, 96, 48, 24]

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")

X_train = train_df.iloc[:, 1:25].values.astype(np.float32)
y_train = (train_df.iloc[:, 25].values == "Virus").astype(np.int32)

X_val = vali_df.iloc[:, 1:25].values.astype(np.float32)
y_val = (vali_df.iloc[:, 25].values == "Virus").astype(np.int32)

model = Sequential()
model.add(Input(shape=(24,)))

for units in hidden_units:
    model.add(Dense(units, activation="relu", kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(Dropout(dropout_rate))

model.add(Dense(1, activation="sigmoid"))

model.compile(
    optimizer=Adam(learning_rate=learning_rate),
    loss="binary_crossentropy",
    metrics=["accuracy"]
)

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
best_val_loss = float("inf")

for epoch in tqdm(range(epochs), desc="Training Epochs"):
    history = model.fit(
        X_train,
        y_train,
        validation_data=(X_val, y_val),
        epochs=1,
        batch_size=batch_size,
        verbose=0
    )

    train_prob = model.predict(X_train, verbose=0).reshape(-1)
    val_prob = model.predict(X_val, verbose=0).reshape(-1)

    train_pred = (train_prob >= 0.5).astype(int)
    val_pred = (val_prob >= 0.5).astype(int)

    train_loss = log_loss(y_train, train_prob)
    val_loss = log_loss(y_val, val_prob)

    train_accuracy = accuracy_score(y_train, train_pred)
    val_accuracy = accuracy_score(y_val, val_pred)

    train_precision = precision_score(y_train, train_pred, zero_division=0)
    val_precision = precision_score(y_val, val_pred, zero_division=0)

    train_recall = recall_score(y_train, train_pred, zero_division=0)
    val_recall = recall_score(y_val, val_pred, zero_division=0)

    train_f1 = f1_score(y_train, train_pred, zero_division=0)
    val_f1 = f1_score(y_val, val_pred, zero_division=0)

    train_mcc = matthews_corrcoef(y_train, train_pred)
    val_mcc = matthews_corrcoef(y_val, val_pred)

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

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        model.save(os.path.join(outdir, "fnn_best_f1_model.keras"))

    if val_precision > best_val_precision:
        best_val_precision = val_precision
        model.save(os.path.join(outdir, "fnn_best_precision_model.keras"))

    if val_recall > best_val_recall:
        best_val_recall = val_recall
        model.save(os.path.join(outdir, "fnn_best_recall_model.keras"))

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        model.save(os.path.join(outdir, "fnn_best_loss_model.keras"))

    print(
        f"Epoch [{epoch + 1}/{epochs}], "
        f"Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}, "
        f"Train Acc: {train_accuracy:.6f}, Val Acc: {val_accuracy:.6f}, "
        f"Train Precision: {train_precision:.6f}, Val Precision: {val_precision:.6f}, "
        f"Train Recall: {train_recall:.6f}, Val Recall: {val_recall:.6f}, "
        f"Train F1: {train_f1:.6f}, Val F1: {val_f1:.6f}, "
        f"Train MCC: {train_mcc:.6f}, Val MCC: {val_mcc:.6f}"
    )

model.save(os.path.join(outdir, "fnn_final_model.keras"))

metrics_df = pd.DataFrame(metrics_history)
metrics_df.to_csv(os.path.join(outdir, "fnn_training_metrics.csv"), index=False)

val_prob = model.predict(X_val, verbose=0).reshape(-1)
val_pred = (val_prob >= 0.5).astype(int)

train_prob = model.predict(X_train, verbose=0).reshape(-1)
train_pred = (train_prob >= 0.5).astype(int)

final_metrics = {
    "train_loss": float(log_loss(y_train, train_prob)),
    "val_loss": float(log_loss(y_val, val_prob)),
    "train_accuracy": float(accuracy_score(y_train, train_pred)),
    "val_accuracy": float(accuracy_score(y_val, val_pred)),
    "train_precision": float(precision_score(y_train, train_pred, zero_division=0)),
    "val_precision": float(precision_score(y_val, val_pred, zero_division=0)),
    "train_recall": float(recall_score(y_train, train_pred, zero_division=0)),
    "val_recall": float(recall_score(y_val, val_pred, zero_division=0)),
    "train_f1": float(f1_score(y_train, train_pred, zero_division=0)),
    "val_f1": float(f1_score(y_val, val_pred, zero_division=0)),
    "train_mcc": float(matthews_corrcoef(y_train, train_pred)),
    "val_mcc": float(matthews_corrcoef(y_val, val_pred))
}

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "fnn_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "fnn_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_df.copy()
vali_pred_df["true_label"] = y_val
vali_pred_df["virus_probability"] = val_prob
vali_pred_df["pred_label"] = val_pred
vali_pred_df["pred_class"] = np.where(val_pred == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "fnn_validation_predictions.csv"), index=False)

plt.figure()
plt.plot(metrics_history["epoch"], metrics_history["train_loss"], label="Train Loss")
plt.plot(metrics_history["epoch"], metrics_history["val_loss"], label="Validation Loss")
plt.title("Loss vs Epoch")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.savefig(os.path.join(outdir, "fnn_loss_vs_epoch.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(y_val, val_prob)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "fnn_roc_curve_values.csv"), index=False)

plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.savefig(os.path.join(outdir, "fnn_roc_curve.pdf"), format="pdf")
plt.close()

cm = confusion_matrix(y_val, val_pred)
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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "fnn_confusion_matrix_values.csv"), index=False)

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
plt.savefig(os.path.join(outdir, "fnn_confusion_matrix.pdf"), format="pdf")
plt.close()

print("Final validation metrics:")
print(final_metrics)
print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")