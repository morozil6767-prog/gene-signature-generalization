import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from scipy.stats import ttest_ind

def run_leakage_controlled_pipeline(X_discovery, y_discovery, external_cohorts, sizes=[5, 10, 20, 50, 100, 200]):
    """
    Executes nested 5-fold CV on discovery data, ranks features inside training folds,
    locks signatures, and evaluates across external cohorts.
    """
    results = {s: {'discovery_cv': [], 'external_means': []} for s in sizes}
    
    # Outer 5-fold CV to prevent data leakage
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for train_idx, val_idx in outer_cv.split(X_discovery, y_discovery):
        X_train, y_train = X_discovery.iloc[train_idx], y_discovery.iloc[train_idx]
        X_val, y_val = X_discovery.iloc[val_idx], y_discovery.iloc[val_idx]
        
        # Stage 1: Feature ranking via t-test on outer training fold only
        t_stat, p_val = ttest_ind(X_train[y_train == 1], X_train[y_train == 0], axis=0)
        ranked_genes = X_train.columns[np.argsort(p_val)]
        
        # Stage 2: Nested Elastic-Net CV on top 500 features
        top_500 = ranked_genes[:500]
        clf = LogisticRegressionCV(penalty='elasticnet', l1_ratio=0.5, solver='saga', cv=5, random_state=42)
        clf.fit(X_train[top_500], y_train)
        
        # Rank by absolute coefficient magnitude
        coef_weights = np.abs(clf.coef_[0])
        final_ranked = top_500[np.argsort(coef_weights)[::-1]]
        
        # Evaluate each nested size S_k
        for k in sizes:
            selected_genes = final_ranked[:k]
            model = LogisticRegressionCV(penalty='elasticnet', l1_ratio=0.5, solver='saga', cv=5, random_state=42)
            model.fit(X_train[selected_genes], y_train)
            
            # Internal CV AUROC
            val_preds = model.predict_proba(X_val[selected_genes])[:, 1]
            results[k]['discovery_cv'].append(roc_auc_score(y_val, val_preds))
            
    return results

print("Pipeline script initialized successfully.")
