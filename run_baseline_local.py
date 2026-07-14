#!/usr/bin/env python3
"""
Everesteer Futures Hackathon — baseline run against locally downloaded parquet.

Same model/eval/output logic as futures_starter.py, but reads train.parquet /
validation.parquet from disk (downloaded via the MCP eiq_download_dataset
tool) instead of instantiating the everestapi SDK client.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent

print("Loading dataset...")
train = pd.read_parquet(HERE / "train.parquet")
val = pd.read_parquet(HERE / "validation.parquet")

EXPED_COL = "exped"
feat_cols = sorted(c for c in train.columns if c.startswith("feature_"))
target_col = "target_everest_20"


def exped_num(e):
    return int(str(e).split("_")[-1])


print(f"  Train:      {len(train):>8,} rows  |  {train[EXPED_COL].nunique()} expeds")
print(f"  Validation: {len(val):>8,} rows  |  {val[EXPED_COL].nunique()} expeds")
print(f"  Features:   {len(feat_cols)}")

# =====================================================================
# Split train into fit + embargoed holdout
# =====================================================================
EMBARGO = 20
expeds = sorted(train[EXPED_COL].unique(), key=exped_num)
holdout_expeds = set(expeds[-100:])
fit_expeds = set(expeds[: -100 - EMBARGO])

fit_df = train[train[EXPED_COL].isin(fit_expeds)]
holdout_df = train[train[EXPED_COL].isin(holdout_expeds)]
print(f"\nFit: {len(fit_df):,} rows | Holdout: {len(holdout_df):,} rows | Embargo: {EMBARGO} expeds")


def features_matrix(df):
    x = df[feat_cols].astype("float32")
    return x.where(x >= 0)


# =====================================================================
# Train a LightGBM model
# =====================================================================
print("\nTraining LightGBM...")

fit_df = fit_df[fit_df[target_col].notna()]

model = lgb.LGBMRegressor(
    n_estimators=2000,
    learning_rate=0.01,
    max_depth=6,
    num_leaves=64,
    colsample_bytree=0.10,
    subsample=0.80,
    min_child_samples=500,
    reg_alpha=0.1,
    reg_lambda=1.0,
    random_state=42,
    verbose=-1,
)
model.fit(features_matrix(fit_df), fit_df[target_col])

# =====================================================================
# Evaluate on the embargoed holdout
# =====================================================================
print("\nEvaluating on the embargoed holdout...")

holdout_df = holdout_df[holdout_df[target_col].notna()].copy()
holdout_df["prediction"] = model.predict(features_matrix(holdout_df))

corrs = []
for _, grp in holdout_df.groupby(EXPED_COL):
    if len(grp) >= 5:
        rho, _ = spearmanr(grp["prediction"], grp[target_col])
        if np.isfinite(rho):
            corrs.append(rho)

corr_arr = np.array(corrs)
print(f"  Mean CORR:     {corr_arr.mean():.4f}")
print(f"  Std CORR:      {corr_arr.std():.4f}")
print(f"  % Positive:    {(corr_arr > 0).mean() * 100:.1f}%")
print(f"  Sharpe (CORR): {corr_arr.mean() / corr_arr.std():.2f}")

# =====================================================================
# Predict the validation split (the leaderboard set) and save
# =====================================================================
print("\nGenerating validation predictions...")

val_preds = model.predict(features_matrix(val))

val_ids = val["id"] if "id" in val.columns else val.index
submission = pd.DataFrame({"prediction": val_preds}, index=pd.Index(val_ids, name="id"))

model_path = HERE / "baseline_model.pkl"
with open(model_path, "wb") as f:
    pickle.dump({"model": model, "feature_cols": feat_cols}, f)
print(f"\n  Model saved: {model_path} ({model_path.stat().st_size / 1024:.0f} KB)")

pred_path = HERE / "baseline_predictions.parquet"
submission.to_parquet(pred_path)
submission.to_csv(pred_path.with_suffix(".csv"))
print(f"  Predictions saved: {pred_path} ({len(submission)} rows)")

print("\nDone.")
