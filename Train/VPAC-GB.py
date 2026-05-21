import os
import json
import joblib
import random
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns

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
y_train = train_df.iloc[:, 25].values
y_train = (y_train == "Virus").astype(int)

X_val = vali_df.iloc[:, 1:25].values
y_val = vali_df.iloc[:, 25].values
y_val = (y_val == "Virus").astype(int)

param_dist = {
    "learning_rate": [0.01, 0.05, 0.1, 0.2, 1],
    "n_estimators": [500, 1000, 2000],
    "max_depth": [None, 3, 5, 7, 9],
    "subsample": [0.8, 1.0],
    "min_child_weight": [1, 3, 5],
    "colsample_bytree": [0.6, 0.8, 1.0],
    "gamma": [0, 0.1, 0.2],
    "reg_alpha": [0, 0.5, 1]
}

base_model = xgb.XGBClassifier(
    objective="binary:logistic",
    eval_metric="logloss",
    random_state=seed,
    n_jobs=-1,
    tree_method="hist"
)

random_search = RandomizedSearchCV(
    estimator=base_model,
    param_distributions=param_dist,
    n_iter=100,
    scoring="f1",
    cv=10,
    random_state=seed,
    n_jobs=-1,
    verbose=2,
    return_train_score=True
)

random_search.fit(X_train, y_train)

cv_results = pd.DataFrame(random_search.cv_results_)
cv_results.to_csv(os.path.join(outdir, "random_search_cv_results.csv"), index=False)

best_params_from_search = random_search.best_params_

fixed_best_params = {
    "learning_rate": 0.05,
    "n_estimators": 2000,
    "max_depth": 9,
    "subsample": 1.0,
    "min_child_weight": 1,
    "colsample_bytree": 0.8,
    "gamma": 0.2,
    "reg_alpha": 0
}

with open(os.path.join(outdir, "best_params_from_random_search.json"), "w") as f:
    json.dump(best_params_from_search, f, indent=2)

with open(os.path.join(outdir, "fixed_best_params_used_for_final_model.json"), "w") as f:
    json.dump(fixed_best_params, f, indent=2)

model = xgb.XGBClassifier(
    objective="binary:logistic",
    eval_metric="logloss",
    random_state=seed,
    n_jobs=-1,
    tree_method="hist",
    **fixed_best_params
)

model.fit(
    X_train,
    y_train,
    eval_set=[(X_train, y_train), (X_val, y_val)],
    verbose=False
)

joblib.dump(model, os.path.join(outdir, "best_xgb_model.pkl"))
model.save_model(os.path.join(outdir, "best_xgb_model.json"))

y_train_prob = model.predict_proba(X_train)[:, 1]
y_val_prob = model.predict_proba(X_val)[:, 1]

y_train_pred = (y_train_prob >= 0.5).astype(int)
y_val_pred = (y_val_prob >= 0.5).astype(int)

train_loss = log_loss(y_train, y_train_prob)
val_loss = log_loss(y_val, y_val_prob)

train_accuracy = accuracy_score(y_train, y_train_pred)
val_accuracy = accuracy_score(y_val, y_val_pred)

train_precision = precision_score(y_train, y_train_pred, zero_division=0)
val_precision = precision_score(y_val, y_val_pred, zero_division=0)

train_recall = recall_score(y_train, y_train_pred, zero_division=0)
val_recall = recall_score(y_val, y_val_pred, zero_division=0)

train_f1 = f1_score(y_train, y_train_pred, zero_division=0)
val_f1 = f1_score(y_val, y_val_pred, zero_division=0)

train_mcc = matthews_corrcoef(y_train, y_train_pred)
val_mcc = matthews_corrcoef(y_val, y_val_pred)

final_metrics = {
    "train_loss": float(train_loss),
    "val_loss": float(val_loss),
    "train_accuracy": float(train_accuracy),
    "val_accuracy": float(val_accuracy),
    "train_precision": float(train_precision),
    "val_precision": float(val_precision),
    "train_recall": float(train_recall),
    "val_recall": float(val_recall),
    "train_f1": float(train_f1),
    "val_f1": float(val_f1),
    "train_mcc": float(train_mcc),
    "val_mcc": float(val_mcc)
}

pd.DataFrame([final_metrics]).to_csv(os.path.join(outdir, "final_metrics.csv"), index=False)

with open(os.path.join(outdir, "final_metrics.json"), "w") as f:
    json.dump(final_metrics, f, indent=2)

vali_pred_df = vali_df.copy()
vali_pred_df["true_label"] = y_val
vali_pred_df["virus_probability"] = y_val_prob
vali_pred_df["pred_label"] = y_val_pred
vali_pred_df["pred_class"] = np.where(y_val_pred == 1, "Virus", "Non-Virus")
vali_pred_df.to_csv(os.path.join(outdir, "validation_predictions.csv"), index=False)

threshold_records = []
best_threshold = 0.5
best_f1 = -1.0

for threshold in np.arange(0.1, 0.9, 0.01):
    pred = (y_val_prob >= threshold).astype(int)
    precision = precision_score(y_val, pred, zero_division=0)
    recall = recall_score(y_val, pred, zero_division=0)
    f1 = f1_score(y_val, pred, zero_division=0)
    mcc = matthews_corrcoef(y_val, pred)
    accuracy = accuracy_score(y_val, pred)

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
threshold_df.to_csv(os.path.join(outdir, "threshold_search.csv"), index=False)

y_val_pred_best = (y_val_prob >= best_threshold).astype(int)

best_threshold_metrics = {
    "best_threshold": float(best_threshold),
    "accuracy": float(accuracy_score(y_val, y_val_pred_best)),
    "precision": float(precision_score(y_val, y_val_pred_best, zero_division=0)),
    "recall": float(recall_score(y_val, y_val_pred_best, zero_division=0)),
    "f1": float(f1_score(y_val, y_val_pred_best, zero_division=0)),
    "mcc": float(matthews_corrcoef(y_val, y_val_pred_best))
}

pd.DataFrame([best_threshold_metrics]).to_csv(os.path.join(outdir, "best_threshold_metrics.csv"), index=False)

with open(os.path.join(outdir, "best_threshold_metrics.json"), "w") as f:
    json.dump(best_threshold_metrics, f, indent=2)

cm = confusion_matrix(y_val, y_val_pred)
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

pd.DataFrame([confusion_metrics]).to_csv(os.path.join(outdir, "confusion_matrix_values.csv"), index=False)

plt.figure(figsize=(6, 5.5))
ax = sns.heatmap(
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
plt.savefig(os.path.join(outdir, "xgb_confusion_matrix.pdf"), format="pdf")
plt.close()

fpr, tpr, _ = roc_curve(y_val, y_val_prob)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "fpr": fpr,
    "tpr": tpr
})
roc_df.to_csv(os.path.join(outdir, "roc_curve_values.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, lw=2, label=f"ROC curve (AUC = {roc_auc:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "xgb_roc_curve.pdf"), format="pdf")
plt.close()

feature_importance = pd.DataFrame({
    "feature": train_df.columns[1:25],
    "importance": model.feature_importances_
}).sort_values("importance", ascending=False)

feature_importance.to_csv(os.path.join(outdir, "feature_importance.csv"), index=False)

plt.figure(figsize=(8, 6))
plt.barh(feature_importance["feature"][::-1], feature_importance["importance"][::-1])
plt.xlabel("Feature Importance")
plt.ylabel("Feature")
plt.title("XGBoost Feature Importance")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "xgb_feature_importance.pdf"), format="pdf")
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