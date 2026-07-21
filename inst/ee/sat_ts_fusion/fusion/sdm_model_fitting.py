# ════════════════════════════════════════════════════════════════════════════
# SDM RF MODEL FITTING PIPELINE
# All outputs named after MODEL_ID — multiple runs never overwrite each other.
#
# Usage:
#   python sdm_model_fitting.py --species puma --model_id puma_modis_m2 \
#       --working_dir C:/Users/pj276/Projects/glad_ard \
#       [--ee_project geersprocessing] [--mode binary] ...
# ════════════════════════════════════════════════════════════════════════════

import numpy as np
import pandas as pd
import pickle
import sys
import os
import ee
import argparse
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import (roc_auc_score, r2_score,
                             mean_absolute_error, balanced_accuracy_score)
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import LabelEncoder
from mrmr import mrmr_classif, mrmr_regression
import warnings
import time
warnings.filterwarnings('ignore')

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM RF Model Fitting Pipeline')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--species',      type=str, required=True)
    p.add_argument('--model_id',     type=str, required=True,
                   help='Unique model run ID, e.g. puma_modis_m2. '
                        'All local and GEE outputs are named after this.')

    # ── Mode ──────────────────────────────────────────────────────────────────
    p.add_argument('--mode',         type=str, default='binary',
                   choices=['binary', 'regression', 'multiclass'])
    p.add_argument('--target_col',   type=str, default='presence')

    # ── Paths ─────────────────────────────────────────────────────────────────
    p.add_argument('--gee_assets',   type=str,
                   default='projects/geersprocessing/assets',
                   help='Root GEE assets folder (no trailing slash)')
    p.add_argument('--working_dir',  type=str, required=True,
                   help='Local working directory for model outputs')
    p.add_argument('--local_csv',    type=str, default=None,
                   help='Path to local training CSV to skip EE asset loading. '
                        'Default: {working_dir}/{model_id}_training.csv')

    # ── Variable selection ────────────────────────────────────────────────────
    p.add_argument('--imp_thresh',          type=float, default=0.10,
                   help='Relative importance threshold for feature pruning')
    p.add_argument('--categorical_threshold', type=int, default=20,
                   help='Max unique values to treat a variable as categorical')

    return p.parse_args()

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# All output paths derived from MODEL_ID.
# =============================================================================

SPECIES      = args.species
MODEL_ID     = args.model_id
MODE         = args.mode
TARGET_COL   = args.target_col
GEE_ASSETS   = args.gee_assets
WORKING_DIR  = args.working_dir
IMP_THRESH   = args.imp_thresh

# Training data asset folder — matches EXPORT_FOLDER from extraction script
TRAINING_ASSET_FOLDER = f'{GEE_ASSETS}/{SPECIES}_SDM_modis_exports_{MODEL_ID}'

# Local outputs — all named after MODEL_ID
LOCAL_CSV       = args.local_csv or os.path.join(WORKING_DIR,
                                                   f'{MODEL_ID}_training.csv')
OUTPUT_MODEL    = os.path.join(WORKING_DIR, f'{MODEL_ID}_rf_final.pkl')
OUTPUT_FEATURES = os.path.join(WORKING_DIR, f'{MODEL_ID}_selected_features.csv')
OUTPUT_SIZING   = os.path.join(WORKING_DIR, f'{MODEL_ID}_rf_sizing.csv')
OUTPUT_SUMMARY  = os.path.join(WORKING_DIR, f'{MODEL_ID}_rf_summary.csv')

# GEE outputs — all named after MODEL_ID
GEE_CLASSIFIER_ID = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier'
GEE_STRINGS_ID    = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier_strings'
GEE_FEATURES_ID   = f'{GEE_ASSETS}/{MODEL_ID}_selected_features'

# =============================================================================
# 2. FIXED PARAMETERS (not user-facing)
# =============================================================================

SAVE_LOCAL_CSV       = True
CATEGORICAL_VARS     = []
CATEGORICAL_THRESHOLD = args.categorical_threshold
FORCE_CONTINUOUS = [
    'lcluc_short_veg', 'lcluc_trees', 'lcluc_wet_short', 'lcluc_wet_trees',
    'lcluc_water', 'lcluc_cropland', 'lcluc_builtup',
    'lcluc_terra_short_pct', 'lcluc_wet_short_pct',
    'lcluc_terra_ht', 'lcluc_wet_ht',
]

INITIAL_N_VARS       = 20
INCREMENT            = 10
ACCURACY_THRESH      = 0.95
CAT_IMP_THRESHOLD    = 0.001
SCALE_CORR_THRESHOLD = 0.80
SIZE_TREES           = [50, 100, 150, 200, 300]
SIZE_NODE_SIZES      = [5, 10, 20, 50]
MAX_MODEL_MB         = 100
USE_SAMPLE_WEIGHTS   = False

META_COLS = [
    'id', 'FID', 'year', 'batch', 'presence', 'NObs',
    'presCount', 'longitude', 'latitude', 'date',
    'system:index', 'system.index', '.geo', 'weight',
    'nObs_clear', 'nObs_priorFill', 'nObs_interp',
    'nObs_clear_250m', 'nObs_clear_1km', 'nObs_clear_5km', 'nObs_clear_10km',
    'nObs_priorFill_250m', 'nObs_priorFill_1km',
    'nObs_priorFill_5km', 'nObs_priorFill_10km',
    'nObs_interp_250m', 'nObs_interp_1km',
    'nObs_interp_5km', 'nObs_interp_10km',
    'lcluc_class',
]

# =============================================================================
# 3. MODE CONFIGURATION
# =============================================================================

MODE_CONFIG = {
    'binary': {
        'metric_name'  : 'ORSS',
        'scoring'      : 'roc_auc',
        'rf_class'     : RandomForestClassifier,
        'rf_extra'     : {},
        'is_classifier': True
    },
    'regression': {
        'metric_name'  : 'R²',
        'scoring'      : 'r2',
        'rf_class'     : RandomForestRegressor,
        'rf_extra'     : {},
        'is_classifier': False
    },
    'multiclass': {
        'metric_name'  : 'Balanced Accuracy',
        'scoring'      : 'balanced_accuracy',
        'rf_class'     : RandomForestClassifier,
        'rf_extra'     : {'class_weight': 'balanced'},
        'is_classifier': True
    }
}
cfg = MODE_CONFIG[MODE]

# =============================================================================
# 4. DATA LOADING
# =============================================================================

def list_batch_assets(asset_folder):
    try:
        assets = ee.data.listAssets({'parent': asset_folder})
        ids    = [a['name'] for a in assets.get('assets', [])]
        print(f'Found {len(ids)} batch assets in {asset_folder}')
        return ids
    except Exception as e:
        print(f'Error listing assets: {e}')
        return []


def fc_to_dataframe(asset_id, max_retries=3):
    for attempt in range(max_retries):
        try:
            fc   = ee.FeatureCollection(asset_id)
            info = fc.getInfo()
            if not info or 'features' not in info:
                return None
            return pd.DataFrame([f['properties'] for f in info['features']])
        except Exception as e:
            if attempt < max_retries - 1:
                print(f'  Retry {attempt+1}/{max_retries}: {e}')
                time.sleep(5)
            else:
                print(f'  Failed: {asset_id}: {e}')
                return None


def load_training_data():
    if os.path.exists(LOCAL_CSV):
        print(f'Loading from local CSV: {LOCAL_CSV}')
        df = pd.read_csv(LOCAL_CSV)
        print(f'Loaded {len(df)} rows x {df.shape[1]} columns')
        return df

    print(f'Loading from EE assets: {TRAINING_ASSET_FOLDER}')
    asset_ids = list_batch_assets(TRAINING_ASSET_FOLDER)
    if not asset_ids:
        raise ValueError(f'No assets found in {TRAINING_ASSET_FOLDER}')

    dfs, failed = [], []
    for i, asset_id in enumerate(asset_ids):
        short = asset_id.split('/')[-1]
        print(f'  [{i+1}/{len(asset_ids)}] {short}', end='... ')
        df_batch = fc_to_dataframe(asset_id)
        if df_batch is not None and len(df_batch) > 0:
            dfs.append(df_batch)
            print(f'{len(df_batch)} rows')
        else:
            failed.append(short)
            print('FAILED or empty')

    if not dfs:
        raise ValueError('No data loaded.')

    df = pd.concat(dfs, ignore_index=True)
    print(f'\nMerged: {len(df)} rows x {df.shape[1]} columns')

    if failed:
        print(f'\n⚠️  {len(failed)} assets failed: {failed}')

    if SAVE_LOCAL_CSV:
        df.to_csv(LOCAL_CSV, index=False)
        print(f'Saved: {LOCAL_CSV}')

    return df

# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================

def identify_variable_types(X):
    scale_suffixes = ('_250m', '_1km', '_5km', '_10km')
    categorical, continuous = [], []
    for col in X.columns:
        base = col
        for sfx in scale_suffixes:
            if col.endswith(sfx):
                base = col[:-len(sfx)]
                break
        if base in FORCE_CONTINUOUS or col in FORCE_CONTINUOUS:
            continuous.append(col)
        elif base in CATEGORICAL_VARS or col in CATEGORICAL_VARS:
            categorical.append(col)
        elif X[col].nunique() <= CATEGORICAL_THRESHOLD:
            categorical.append(col)
        else:
            continuous.append(col)
    return categorical, continuous


def get_rf_model(n_trees, node_size):
    return cfg['rf_class'](
        n_estimators=n_trees, min_samples_leaf=node_size,
        n_jobs=-1, random_state=42, oob_score=True, **cfg['rf_extra'])


def _youden_threshold(y_true, y_prob):
    from sklearn.metrics import roc_curve
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    idx = np.argmax(tpr - fpr)
    return thresholds[idx], tpr[idx], 1.0 - fpr[idx]


def compute_orss(y_true, y_prob):
    from sklearn.metrics import confusion_matrix
    thresh, _, _ = _youden_threshold(y_true, y_prob)
    tn, fp, fn, tp = confusion_matrix(
        y_true, (y_prob >= thresh).astype(int)).ravel()
    denom = tp * tn + fp * fn
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0


def compute_sedi(y_true, y_prob):
    _, H, spec = _youden_threshold(y_true, y_prob)
    F   = 1.0 - spec
    eps = 1e-10
    H, F = np.clip(H, eps, 1-eps), np.clip(F, eps, 1-eps)
    num   = np.log(F) - np.log(H) + np.log(1-H) - np.log(1-F)
    denom = np.log(F) + np.log(H) + np.log(1-H) + np.log(1-F)
    return -num / denom if denom != 0 else 0.0


def compute_oob_metric(rf, X, y):
    if MODE == 'binary':
        return compute_orss(y, rf.oob_decision_function_[:, 1])
    elif MODE == 'regression':
        return r2_score(y, rf.oob_prediction_)
    elif MODE == 'multiclass':
        return balanced_accuracy_score(
            y, np.argmax(rf.oob_decision_function_, axis=1))


def fit_rf(X, y, n_trees=300, node_size=5, sample_weight=None):
    rf = get_rf_model(n_trees, node_size)
    rf.fit(X, y, sample_weight=sample_weight)
    return rf, compute_oob_metric(rf, X, y)


def get_permutation_importance(rf, X, y, n_repeats=5):
    if MODE == 'binary':
        def scorer(estimator, X, y):
            return compute_orss(y, estimator.predict_proba(X)[:, 1])
    else:
        scorer = cfg['scoring']
    result = permutation_importance(rf, X, y, n_repeats=n_repeats,
                                    random_state=42, n_jobs=1, scoring=scorer)
    imp = pd.DataFrame({'feature': X.columns,
                        'importance': result.importances_mean}
                       ).sort_values('importance', ascending=False)
    total = imp['importance'].clip(lower=0).sum()
    imp['importance_pct'] = imp['importance'].clip(lower=0) / total if total > 0 else 0.0
    return imp


def remove_low_importance(features, imp_df):
    mean_imp      = imp_df[imp_df['importance_pct'] > 0]['importance_pct'].mean()
    abs_threshold = IMP_THRESH * mean_imp
    low_imp   = imp_df[imp_df['importance_pct'] < abs_threshold]['feature'].tolist()
    remaining = [f for f in features if f not in low_imp]
    removed   = [f for f in features if f in low_imp]
    print(f'  Importance threshold: {IMP_THRESH:.0%} of mean '
          f'({mean_imp:.3%}) = {abs_threshold:.3%}')
    return remaining, removed


def get_model_size_mb(rf):
    return sys.getsizeof(pickle.dumps(rf)) / 1e6


def get_sample_weights(y):
    if not USE_SAMPLE_WEIGHTS:
        return None
    classes, counts = np.unique(y, return_counts=True)
    class_weights   = dict(zip(classes, 1.0 / counts))
    weights = np.array([class_weights[c] for c in y])
    return weights / weights.mean()


def compute_additional_metrics(rf, X, y):
    if MODE == 'binary':
        from sklearn.metrics import confusion_matrix
        oob_probs = rf.oob_decision_function_[:, 1]
        opt_thresh, sensitivity, specificity = _youden_threshold(y, oob_probs)
        predicted  = (oob_probs >= opt_thresh).astype(int)
        tn, fp, fn, tp = confusion_matrix(y, predicted).ravel()
        return {
            'threshold'  : opt_thresh,
            'sensitivity': sensitivity,
            'specificity': specificity,
            'accuracy'   : (tp + tn) / len(y),
            'orss'       : compute_orss(y, oob_probs),
            'tss'        : sensitivity + specificity - 1,
            'auc'        : roc_auc_score(y, oob_probs),
            'sedi'       : compute_sedi(y, oob_probs),
        }
    elif MODE == 'regression':
        oob_preds  = rf.oob_prediction_
        return {
            'rmse': float(np.sqrt(np.mean((y - oob_preds)**2))),
            'mae' : float(mean_absolute_error(y, oob_preds)),
            'bias': float(np.mean(oob_preds - y)),
        }
    elif MODE == 'multiclass':
        from sklearn.metrics import classification_report, confusion_matrix
        oob_preds = np.argmax(rf.oob_decision_function_, axis=1)
        return {
            'balanced_accuracy'    : balanced_accuracy_score(y, oob_preds),
            'classification_report': classification_report(y, oob_preds),
            'confusion_matrix'     : confusion_matrix(y, oob_preds),
        }


def print_metrics(metrics):
    if MODE == 'binary':
        print(f'  Threshold:   {metrics["threshold"]:.3f}')
        print(f'  Sensitivity: {metrics["sensitivity"]:.3f}')
        print(f'  Specificity: {metrics["specificity"]:.3f}')
        print(f'  Accuracy:    {metrics["accuracy"]:.3f}')
        print(f'  ORSS:        {metrics["orss"]:.3f}  <- primary')
        print(f'  TSS:         {metrics["tss"]:.3f}')
        print(f'  AUC:         {metrics["auc"]:.3f}')
        print(f'  SEDI:        {metrics["sedi"]:.3f}')
    elif MODE == 'regression':
        print(f'  RMSE: {metrics["rmse"]:.3f}')
        print(f'  MAE:  {metrics["mae"]:.3f}')
        print(f'  Bias: {metrics["bias"]:.3f}')
    elif MODE == 'multiclass':
        print(f'  Balanced accuracy: {metrics["balanced_accuracy"]:.3f}')
        print(metrics['classification_report'])

# =============================================================================
# 6. MRMR + VARIABLE SELECTION
# =============================================================================

def get_base_name(col):
    for sfx in ('_250m', '_1km', '_5km', '_10km'):
        if col.endswith(sfx):
            return col[:-len(sfx)]
    return col


def run_mrmr_with_scale_filter(X_cont, y, n_vars):
    k_request = len(X_cont.columns)
    try:
        if MODE == 'regression':
            ranked = mrmr_regression(X=X_cont, y=pd.Series(y),
                                     K=k_request, relevance='f', redundancy='c')
        else:
            ranked = mrmr_classif(X=X_cont, y=pd.Series(y),
                                   K=k_request, relevance='f', redundancy='c')
    except Exception as e:
        print(f'mRMR error: {e}. Falling back to RF importance.')
        rf_imp = get_rf_model(200, 5)
        rf_imp.fit(X_cont, y)
        ranked = X_cont.columns[
            np.argsort(rf_imp.feature_importances_)[::-1][:k_request]].tolist()

    selected, total_excluded = [], []
    for candidate in ranked:
        if len(selected) >= n_vars:
            break
        corr_too_high = False
        for sel_var in selected:
            r = abs(np.corrcoef(X_cont[candidate].values,
                                X_cont[sel_var].values)[0, 1])
            if r > SCALE_CORR_THRESHOLD:
                total_excluded.append((candidate, r, sel_var))
                corr_too_high = True
                break
        if not corr_too_high:
            selected.append(candidate)

    if len(selected) < n_vars:
        print(f'  Note: only {len(selected)} features pass correlation filter')
    if total_excluded:
        print(f'  Scale filter excluded {len(total_excluded)} redundant variables:')
        for var, r, ref in total_excluded[:10]:
            print(f'    {var} (r={r:.3f} with {ref})')
        if len(total_excluded) > 10:
            print(f'    ... and {len(total_excluded)-10} more')
    return selected[:n_vars]


def test_categorical_features(selected, X, y, cat_cols, sample_weights=None):
    if not cat_cols:
        print('  No categorical features to test')
        rf, metric = fit_rf(X[selected], y, sample_weight=sample_weights)
        return selected, rf, metric

    print(f'\nTesting {len(cat_cols)} categorical features:')
    rf_base, metric_base = fit_rf(X[selected], y, sample_weight=sample_weights)
    print(f'  Baseline {cfg["metric_name"]}: {metric_base:.4f}')
    current_features, current_metric, current_rf = list(selected), metric_base, rf_base
    added_cats = []

    for cat in cat_cols:
        test_features    = current_features + [cat]
        rf_test, metric_test = fit_rf(X[test_features], y,
                                       sample_weight=sample_weights)
        improvement = metric_test - current_metric
        added = improvement > CAT_IMP_THRESHOLD
        print(f'  {cat}: {cfg["metric_name"]} = {metric_test:.4f} '
              f'({improvement:+.4f}) {"✅ added" if added else "❌ skipped"}')
        if added:
            current_features, current_metric, current_rf = (
                test_features, metric_test, rf_test)
            added_cats.append(cat)

    if added_cats:
        print(f'\n  Added {len(added_cats)} categorical features: {added_cats}')
    else:
        print('\n  No categorical features improved the model')
    return current_features, current_rf, current_metric


def run_variable_selection(X, y, metric_full, metric_target,
                            sample_weights, cat_cols, cont_cols):
    X_cont    = X[cont_cols]
    n_vars    = min(INITIAL_N_VARS, len(cont_cols))
    iteration = 0

    while n_vars <= len(cont_cols):
        iteration += 1
        print(f'\n--- Iteration {iteration}: trying {n_vars} MRMR features ---')
        selected_cont = run_mrmr_with_scale_filter(X_cont, y, n_vars)
        print(f'Selected: {selected_cont}')
        rf, metric_sel = fit_rf(X[selected_cont], y, sample_weight=sample_weights)
        print(f'{cfg["metric_name"]}: {metric_sel:.4f} (target: {metric_target:.4f})')

        if metric_sel >= metric_target:
            print('✅ Meets target')
            imp = get_permutation_importance(rf, X[selected_cont], y)
            print(f'\nPermutation importance (top {min(20, len(imp))}):')
            print(imp.head(20).to_string(index=False))

            remaining, removed = remove_low_importance(selected_cont, imp)
            if removed:
                print(f'\nRemoving {len(removed)} low importance features:')
                for r in removed:
                    pct = imp[imp['feature'] == r]['importance_pct'].values
                    print(f'  {r}: {pct[0]:.3%}' if len(pct) > 0 else f'  {r}')
                rf2, metric_rem = fit_rf(X[remaining], y, sample_weight=sample_weights)
                print(f'{cfg["metric_name"]} after removal: {metric_rem:.4f}')
                if metric_rem >= metric_target:
                    print('✅ Still meets target')
                    selected_cont, rf, metric_sel = remaining, rf2, metric_rem
                else:
                    print('⚠️  Dropped below target — keeping original selection')
            else:
                print('No low importance features to remove')

            return test_categorical_features(selected_cont, X, y,
                                             cat_cols, sample_weights)
        else:
            print('❌ Below target')
            if n_vars >= len(cont_cols):
                print('Cannot add more features. Testing categoricals anyway...')
                return test_categorical_features(selected_cont, X, y,
                                                 cat_cols, sample_weights)
            else:
                n_vars = min(n_vars + INCREMENT, len(cont_cols))
                print(f'Adding {INCREMENT} → trying {n_vars}')

    return selected_cont, rf, metric_sel

# =============================================================================
# 7. MAIN
# =============================================================================

if __name__ == '__main__':
    print('=' * 60)
    print(f'SDM RF MODEL FITTING PIPELINE')
    print(f'Species: {SPECIES.upper()} | Model ID: {MODEL_ID} | Mode: {MODE.upper()}')
    print(f'Training assets: {TRAINING_ASSET_FOLDER}')
    print(f'Working dir:     {WORKING_DIR}')
    print(f'GEE classifier:  {GEE_CLASSIFIER_ID}')
    print('=' * 60)

    os.chdir(WORKING_DIR)

    # ── STEP 1: LOAD DATA ─────────────────────────────────────────────────────
    print('\nSTEP 1: Loading data')
    print('-' * 40)
    df = load_training_data()

    feature_cols = [c for c in df.columns if c not in META_COLS and c != TARGET_COL]
    X_full       = df[feature_cols].copy()
    y_raw        = df[TARGET_COL].values

    if MODE == 'multiclass':
        le = LabelEncoder()
        y  = le.fit_transform(y_raw)
        print(f'Classes: {le.classes_}')
    else:
        y, le = y_raw.astype(float), None

    print(f'Features before cleaning: {len(feature_cols)}')
    if MODE == 'binary':
        print(f'Presences: {int(y.sum())} ({y.mean():.1%}) | '
              f'Absences: {int((1-y).sum())}')

    # Drop non-numeric
    string_cols = X_full.select_dtypes(exclude=[np.number]).columns.tolist()
    if string_cols:
        print(f'Dropping non-numeric: {string_cols}')
        X_full = X_full.drop(columns=string_cols)

    X_full = X_full.astype(np.float64)

    # Drop all-sentinel columns
    sentinel_cols = [c for c in X_full.columns if (X_full[c] == -9999).all()]
    if sentinel_cols:
        print(f'Dropping {len(sentinel_cols)} all-sentinel columns')
        X_full = X_full.drop(columns=sentinel_cols)

    # Fill sentinels
    FILL_ZERO_COLS = ['lcluc_terra_ht', 'lcluc_wet_ht',
                      'lcluc_terra_short_pct', 'lcluc_wet_short_pct']
    sentinel_mask = X_full == -9999
    if sentinel_mask.any().any():
        n_sentinel = int(sentinel_mask.sum().sum())
        print(f'Replacing {n_sentinel} sentinels')
        for col in X_full.columns:
            mask = X_full[col] == -9999
            if not mask.any(): continue
            base = col
            for sfx in ('_250m','_1km','_5km','_10km'):
                if col.endswith(sfx): base = col[:-len(sfx)]; break
            if base in FILL_ZERO_COLS or col in FILL_ZERO_COLS:
                X_full.loc[mask, col] = 0.0
            else:
                X_full.loc[mask, col] = X_full.loc[~mask, col].median()

    # Impute NAs
    na_count = X_full.isna().sum().sum()
    if na_count > 0:
        print(f'Imputing {na_count} NAs with column median')
        X_full = X_full.fillna(X_full.median(numeric_only=True))

    # Remove near-zero variance
    nzv = X_full.var()[X_full.var() < 1e-6].index.tolist()
    if nzv:
        print(f'Removing {len(nzv)} near-zero variance features')
        X_full = X_full.drop(columns=nzv)

    feature_cols = X_full.columns.tolist()
    print(f'Features after cleaning: {len(feature_cols)}')

    cat_cols, cont_cols = identify_variable_types(X_full)
    print(f'Categorical: {len(cat_cols)} | Continuous: {len(cont_cols)}')
    sample_weights = get_sample_weights(y)

    # ── STEP 2: FULL MODEL BASELINE ───────────────────────────────────────────
    print(f'\nSTEP 2: Full model baseline')
    print('-' * 40)
    rf_full, metric_full = fit_rf(X_full, y, n_trees=300, node_size=5,
                                   sample_weight=sample_weights)
    metric_target = metric_full * ACCURACY_THRESH
    print(f'Full {cfg["metric_name"]}: {metric_full:.4f} | '
          f'Target (95%): {metric_target:.4f} | '
          f'Size: {get_model_size_mb(rf_full):.1f} MB')

    # ── STEP 3: VARIABLE SELECTION ────────────────────────────────────────────
    print(f'\nSTEP 3: Variable selection')
    print('-' * 40)
    selected_features, rf_selected, metric_selected = run_variable_selection(
        X_full, y, metric_full, metric_target,
        sample_weights, cat_cols, cont_cols)

    n_cat  = sum(f in cat_cols  for f in selected_features)
    n_cont = sum(f in cont_cols for f in selected_features)

    print('\n' + '=' * 60)
    print('VARIABLE SELECTION COMPLETE')
    print('=' * 60)

    from scipy import stats as _stats
    def _assoc(f, X, y):
        vals = X[f].values
        try:
            if f in cat_cols:
                if MODE == 'binary':
                    return None, '+' if vals[y==1].mean() > vals[y==0].mean() else '-'
                return None, '?'
            r, _ = _stats.pearsonr(vals, y.astype(float))
            return r, '+' if r > 0 else '-'
        except Exception:
            return None, '?'

    print(f'Selected {len(selected_features)} features '
          f'({n_cont} continuous, {n_cat} categorical):')
    print(f'  {"#":>2}  {"Type":4}  {"Dir":3}  {"r":>7}  Feature')
    print(f'  {"--"}  {"----"}  {"---"}  {"-------"}  {"------------------------------"}')
    for i, f in enumerate(selected_features, 1):
        ftype    = 'cat' if f in cat_cols else 'cont'
        r, dirn  = _assoc(f, X_full, y)
        r_str    = f'{r:+.3f}' if r is not None else '  N/A '
        print(f'  {i:2d}. [{ftype}]  {dirn:>3}  {r_str}  {f}')

    print(f'\n{cfg["metric_name"]}: {metric_selected:.4f} | '
          f'Full: {metric_full:.4f} | '
          f'Retention: {metric_selected/metric_full:.1%}')

    # ── STEP 4: MODEL SIZING ──────────────────────────────────────────────────
    print(f'\nSTEP 4: Model sizing')
    print('-' * 40)
    X_sel, sizing_results = X_full[selected_features], []
    for n_trees in SIZE_TREES:
        for node_size in SIZE_NODE_SIZES:
            rf = get_rf_model(n_trees, node_size)
            rf.fit(X_sel, y, sample_weight=sample_weights)
            metric = compute_oob_metric(rf, X_sel, y)
            size   = get_model_size_mb(rf)
            valid  = (metric >= metric_target) and (size <= MAX_MODEL_MB)
            print(f'Trees: {n_trees:3d}, node: {node_size:3d} → '
                  f'{cfg["metric_name"]}: {metric:.4f}, '
                  f'Size: {size:5.1f}MB {"✅" if valid else "❌"}')
            sizing_results.append({
                'n_trees': n_trees, 'node_size': node_size,
                'metric': metric, 'size_mb': size,
                'metric_drop': metric_full - metric,
                'meets_metric': metric >= metric_target,
                'meets_size': size <= MAX_MODEL_MB, 'valid': valid
            })

    results_df = pd.DataFrame(sizing_results)
    SIZING_TOLERANCE = 0.05
    valid_df = results_df[
        (results_df['metric'] >= metric_selected - SIZING_TOLERANCE) &
        (results_df['meets_size'])]
    if len(valid_df) == 0:
        print(f'⚠️  No config within tolerance — using best available')
        valid_df = results_df[results_df['meets_size']]
    optimal = valid_df.sort_values(['metric','n_trees'],
                                   ascending=[False,True]).iloc[0]
    print(f'\n✅ Optimal: {int(optimal.n_trees)} trees, '
          f'node {int(optimal.node_size)}, '
          f'{cfg["metric_name"]} {optimal.metric:.4f}, '
          f'{optimal.size_mb:.1f} MB')

    # ── STEP 5: FINAL MODEL ───────────────────────────────────────────────────
    print(f'\nSTEP 5: Final model')
    print('-' * 40)
    rf_final = get_rf_model(int(optimal.n_trees), int(optimal.node_size))
    rf_final.fit(X_sel, y, sample_weight=sample_weights)
    final_metric  = compute_oob_metric(rf_final, X_sel, y)
    final_size    = get_model_size_mb(rf_final)
    final_metrics = compute_additional_metrics(rf_final, X_sel, y)
    print(f'Final {cfg["metric_name"]}: {final_metric:.4f} | Size: {final_size:.1f} MB')
    print_metrics(final_metrics)

    # ── STEP 6: SAVE LOCAL OUTPUTS ────────────────────────────────────────────
    print(f'\nSTEP 6: Saving outputs')
    print('-' * 40)

    with open(OUTPUT_MODEL, 'wb') as f:
        pickle.dump(rf_final, f)
    print(f'Model:    {OUTPUT_MODEL}')

    if le is not None:
        enc_path = OUTPUT_MODEL.replace('.pkl', '_encoder.pkl')
        with open(enc_path, 'wb') as f:
            pickle.dump(le, f)
        print(f'Encoder:  {enc_path}')

    pd.DataFrame({
        'rank'   : range(1, len(selected_features) + 1),
        'feature': selected_features,
        'type'   : ['categorical' if f in cat_cols else 'continuous'
                    for f in selected_features]
    }).to_csv(OUTPUT_FEATURES, index=False)
    print(f'Features: {OUTPUT_FEATURES}')

    results_df.to_csv(OUTPUT_SIZING, index=False)
    print(f'Sizing:   {OUTPUT_SIZING}')

    pd.Series({
        'species'            : SPECIES,
        'model_id'           : MODEL_ID,
        'mode'               : MODE,
        'n_features_full'    : len(feature_cols),
        'n_features_selected': len(selected_features),
        'n_continuous'       : n_cont,
        'n_categorical'      : n_cat,
        'metric_name'        : cfg['metric_name'],
        'metric_full'        : metric_full,
        'metric_selected'    : metric_selected,
        'metric_final'       : final_metric,
        'metric_retention'   : final_metric / metric_full,
        'n_trees'            : int(optimal.n_trees),
        'node_size'          : int(optimal.node_size),
        'model_size_mb'      : final_size,
        'training_assets'    : TRAINING_ASSET_FOLDER,
        'gee_classifier'     : GEE_CLASSIFIER_ID,
        'gee_strings'        : GEE_STRINGS_ID,
        'gee_features'       : GEE_FEATURES_ID,
    }).to_csv(OUTPUT_SUMMARY)
    print(f'Summary:  {OUTPUT_SUMMARY}')

    # ── STEP 7: GEE EXPORT ────────────────────────────────────────────────────
    print(f'\nSTEP 7: GEE export')
    print('-' * 40)
    print(f'  Classifier: {GEE_CLASSIFIER_ID}')
    print(f'  Strings:    {GEE_STRINGS_ID}')
    print(f'  Features:   {GEE_FEATURES_ID}')

    try:
        from geemap import ml
        import types

        output_mode = {'regression': 'REGRESSION',
                       'binary': 'PROBABILITY',
                       'multiclass': 'MULTIPROBABILITY'}[MODE]

        _real_int = int
        def _safe_int(x):
            try: return _real_int(x)
            except (ValueError, TypeError): return _real_int(float(x))

        original_fn     = ml.tree_to_string
        patched_globals = dict(original_fn.__globals__)
        patched_globals['int'] = _safe_int
        ml.tree_to_string = types.FunctionType(
            original_fn.__code__, patched_globals,
            original_fn.__name__, original_fn.__defaults__,
            original_fn.__closure__)

        try:
            ee_strings = ml.rf_to_strings(rf_final, selected_features,
                                          processes=1, output_mode=output_mode)
        finally:
            ml.tree_to_string = original_fn

        # Classifier asset
        task = ee.batch.Export.classifier.toAsset(
            classifier=ml.strings_to_classifier(ee_strings),
            description=GEE_CLASSIFIER_ID.split('/')[-1],
            assetId=GEE_CLASSIFIER_ID)
        task.start()
        print(f'  Classifier export started')

        # Tree strings (preserves output mode for wall-to-wall)
        strings_fc = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point([0, 0]),
                       {'tree_index': i, 'tree_string': s, 'output_mode': output_mode})
            for i, s in enumerate(ee_strings)
        ])
        ee.batch.Export.table.toAsset(
            collection=strings_fc,
            description=f'{MODEL_ID}_RF_strings',
            assetId=GEE_STRINGS_ID).start()
        print(f'  Strings export started')

        # Selected features
        ee.batch.Export.table.toAsset(
            collection=ee.FeatureCollection([
                ee.Feature(ee.Geometry.Point([0, 0]), {
                    'selected_features': ','.join(selected_features),
                    'n_features'        : len(selected_features),
                    'species'           : SPECIES,
                    'model_id'          : MODEL_ID,
                    'mode'              : MODE,
                    'metric_name'       : cfg['metric_name'],
                    'final_metric'      : final_metric,
                })
            ]),
            description=f'{MODEL_ID}_selected_features',
            assetId=GEE_FEATURES_ID).start()
        print(f'  Features export started')

    except ImportError:
        print('geemap not installed — skipping GEE export')
    except Exception as e:
        print(f'GEE export error: {e}')

    # ── FINAL SUMMARY ─────────────────────────────────────────────────────────
    print('\n' + '=' * 60)
    print('PIPELINE COMPLETE')
    print('=' * 60)
    print(f'Species:       {SPECIES}')
    print(f'Model ID:      {MODEL_ID}')
    print(f'Features:      {len(selected_features)} ({n_cont} cont, {n_cat} cat)')
    print(f'Full {cfg["metric_name"]}:   {metric_full:.4f}')
    print(f'Final {cfg["metric_name"]}:  {final_metric:.4f} '
          f'({final_metric/metric_full:.1%} retention)')
    print(f'Trees:         {int(optimal.n_trees)} | '
          f'Node size: {int(optimal.node_size)} | '
          f'Size: {final_size:.1f} MB')
    print(f'\nLocal outputs:')
    print(f'  {OUTPUT_MODEL}')
    print(f'  {OUTPUT_FEATURES}')
    print(f'  {OUTPUT_SUMMARY}')
    print(f'\nGEE outputs:')
    print(f'  {GEE_CLASSIFIER_ID}')
    print(f'  {GEE_STRINGS_ID}')
    print(f'  {GEE_FEATURES_ID}')
    print('\nDetailed metrics:')
    print_metrics(final_metrics)
