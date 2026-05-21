import os
import json
import random
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import reciprocal, expon
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RandomizedSearchCV
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

seed =777
random.seed(seed)
np.random.seed(seed)

train_path = "train_shff.tsv"
vali_path = "vali.tsv"
outdir = "./vali"
os.makedirs(outdir, exist_ok=True)

train_df = pd.read_csv(train_path, sep="\t")
vali_df = pd.read_csv(vali_path, sep="\t")

X_train = train_df.iloc[:, 1:25].values
y_train = (train_df.iloc[:, 25].values == "Virus").astype(int)

X_val = vali_df.iloc[:, 1:25].values
y_val = (vali_df.iloc[:, 25].values == "Virus").astype(int)

param_distributions = [
    {
        "svc__C": reciprocal(0.01, 10),
        "svc__gamma": ["scale", "auto"],
        "svc__degree": [2, 3, 4, 5, 6, 7, 8, 9],
        "svc__coef0": np.linspace(0, 10, 10),
        "svc__shrinking": [True, False]
    },
    {
        "svc__C": reciprocal(0.01, 10),
        "svc__gamma": expon(scale=0.1),
        "svc__degree": [2, 3, 4, 5, 6, 7, 8, 9],
        "svc__coef0": np.linspace(0, 10, 10),
        "svc__shrinking": [True, False]
    }
]

base_model = make_pipeline(
    StandardScaler(),
    SVC(
        kernel="rbf",
        probability=True,
        random_state=seed
    )
)

random_search = RandomizedSearchCV(
    estimator=base_model,
    param_distributions=param_distributions,
    n_iter=200,
    cv=5,
    scoring="f1",
    n_jobs=-1,
    verbose=2,
    random_state=seed,
    return_train_score=True
)

random_search.fit(X_train, y_train)

cv_results = pd.DataFrame(random_search.cv_results_)
cv_results.to_csv(os.path.join(outdir, "svm_random_search_cv_results.csv"), index=False)

best_params_from_search = random_search.best_params_

fixed_best_params = {
    "svc__C": 1.9604,
    "svc__gamma": "auto",
    "svc__degree": 2,
    "svc__coef0": 2.2222,
    "svc__shrinking": True
}

with open(os.path.join(outdir, "svm_best_params_from_random_search.json"), "w") as f:
    json.dump(best_params_from_search, f, indent=2, default=str)

with open(os.path.join(outdir, "svm_fixed_best_params_used_for_final_model.json"), "w") as f:
    json.dump(fixed_best_params, f, indent=2)

model = make_pipeline(
    StandardScaler(),
    SVC(
        kernel="rbf",
        probability=True,
        random_state=seed
    )
)

model.set_params(**fixed_best_params)
model.fit(X_train, y_train)

joblib.dump(model, os.path.join(outdir, "best_svm_model.pkl"))

y_train_prob = model.predict_proba(X_train)[:, 1]
y_val_prob = model.predict_proba(X_val)[:, 1]

y_train_pred = (y_train_prob >= 0.5).astype(int)
y_val_pred = (y_val_prob >= 0.5).astype(int)

train_loss = log_loss(y_train, y_train_prob)
val_loss = log_loss(y_val, y_val_prob)

final_metrics = {
    "train_loss": float(train_loss),
    "val_loss": float(val_loss),
    "train_accuracy": float(accuracy_score(y_train, y_train_pred)),
    "val_accuracy": float(accuracy_score(y_val, y_val_pred)),
    "train_precision": float(precision_score(y_train, y_train_pred, zero_division=0)),
    "val_precision": float(precision_score(y_val, y_val_pred, zero_division=0)),
    "train_recall": float(recall_score(y_train, y_train_pred, zero_division=0)),
    "val_recall": float(recall_score(y_val, y_val_pred, zero_division=0)),
    "train_f1": float(f1_score(y_train, y_train_pred, zero_division=0)),
    "val_f1": float(f1_score(y_val, y_val_pred, zero_division=0)),
    "train_mcc": float(matthews_corrcoef(y_train, y_train_pred)),
    "val_mcc": float(matthews_corrcoef(y_val, y_val_pred))
}

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "svm_final_metrics.csv"), index=False)

with open(os.path.join(outdir, "svm_final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_df.copy()
vali_pred_df["true_label"] = y_val
vali_pred_df["virus_probability"] = y_val_prob
vali_pred_df["pred_label"] = y_val_pred
vali_pred_df["pred_class"] = np.where(y_val_pred == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "svm_validation_predictions.csv"), index=False)

threshold_records = []
best_threshold = 0.5
best_f1 = -1.0

for threshold in np.arange(0.1, 0.9, 0.01):
    pred = (y_val_prob >= threshold).astype(int)

    accuracy = accuracy_score(y_val, pred)
    precision = precision_score(y_val, pred, zero_division=0)
    recall = recall_score(y_val, pred, zero_division=0)
    f1 = f1_score(y_val, pred, zero_division=0)
    mcc = matthews_corrcoef(y_val, pred)

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
threshold_df.to_csv(os.path.join(outdir, "svm_threshold_search.csv"), index=False)

y_val_pred_best = (y_val_prob >= best_threshold).astype(int)

best_threshold_metrics = {
    "best_threshold": float(best_threshold),
    "accuracy": float(accuracy_score(y_val, y_val_pred_best)),
    "precision": float(precision_score(y_val, y_val_pred_best, zero_division=0)),
    "recall": float(recall_score(y_val, y_val_pred_best, zero_division=0)),
    "f1": float(f1_score(y_val, y_val_pred_best, zero_division=0)),
    "mcc": float(matthews_corrcoef(y_val, y_val_pred_best))
}

pd.DataFrame([best_threshold_metrics]).to_csv(os.path.join(outdir, "svm_best_threshold_metrics.csv"), index=False)

with open(os.path.join(outdir, "svm_best_threshold_metrics.json"), "w") as f:
    json.dump(best_threshold_metrics, f, indent=2)

cm = confusion_matrix(y_val, y_val_pred, labels=[0, 1])
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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "svm_confusion_matrix_values.csv"), index=False)

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
plt.savefig(os.path.join(outdir, "svm_confusion_matrix.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(y_val, y_val_prob)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "svm_roc_curve_values.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "svm_roc_curve.pdf"), format="pdf")
plt.close()

print("Best parameters from randomized search:")
print(best_params_from_search)
print("Fixed best parameters used for final model:")
print(fixed_best_params)
print("Final validation metrics:")
print(final_metrics)
print("Best threshold metrics:")
print(best_threshold_metrics)
print(f"ROC AUC: {roc_auc:.4f}")
print(f"True Positive (TP): {TP} ({TP / total:.2%})")
print(f"True Negative (TN): {TN} ({TN / total:.2%})")
print(f"False Positive (FP): {FP} ({FP / total:.2%})")
print(f"False Negative (FN): {FN} ({FN / total:.2%})")