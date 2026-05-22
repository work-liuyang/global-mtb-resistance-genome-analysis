import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use ('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
from scipy import stats
from scipy.stats import ttest_rel
import itertools
import sys
from statsmodels.stats.multitest import multipletests

from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (f1_score, precision_score, recall_score, 
                            average_precision_score, precision_recall_curve)

from imblearn.ensemble import BalancedRandomForestClassifier

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'


class DrugFeatures:
    class INH_feat:
        AssowR_predominant = [
            "katG_p.Ser315Thr","n.1673425C>T"
        ]
        AssowR = [
            "katG_p.Ser315Thr","n.1673425C>T","katG_p.Ser315Asn","n.1674048G>A",
            "katG_frameshift_variant","n.1673432T>C","n.1673432T>A","katG_stop_lost"
        ]
        AssowR_and_Uncertain = [
            "katG_p.Ser315Thr","n.1673425C>T","katG_p.Ser315Asn","n.2726141C>T","n.1674048G>A",
            "n.2726145G>A","n.2726139C>T","n.2726112C>T","katG_p.Trp300Gly","n.2726136C>T",
            "inhA_p.Ser94Ala","katG_p.Asn138Asp","n.2156121T>G","n.1673423G>T","n.1673432T>C",
            "katG_p.Ala109Val","n.1673432T>A","n.2726119G>A","katG_p.Ser315Gly","katG_frameshift_variant",
            "katG_stop_lost","katG_p.Tyr337Cys","katG_p.Ser315Arg","katG_p.Tyr155Ser","katG_p.Gly234Glu",
            "katG_p.Gly285Asp","katG_p.Met126Ile","katG_p.Thr394Ala","katG_p.Trp91Arg","katG_p.Gln127Pro",
            "inhA_p.Ile21Thr","inhA_p.Ile194Thr","n.2156118C>T","katG_p.Ser140Asn","katG_p.Ala106Val",
            "katG_p.Thr308Pro","katG_p.Trp191Arg"
        ]

    class RIF_feat:
        AssowR_predominant = [
            "rpoB_p.Ser450Leu","rpoB_p.Asp435Val"
        ]
        AssowR = [
            "rpoB_p.Ser450Leu","rpoB_p.His445Tyr","rpoB_p.Leu452Pro","rpoB_p.His445Asp","rpoB_p.Asp435Val",
            "rpoB_p.Leu430Pro","rpoB_p.His445Leu","rpoB_p.Ser450Trp","rpoB_p.Asp435Tyr","rpoB_p.His445Arg",
            "rpoB_p.His445Asn","rpoB_p.Val170Phe","rpoB_p.Ser441Leu","rpoB_p.Ile491Phe","rpoB_p.His445Ser",
            "rpoB_p.Ser450Phe","rpoB_p.His445Cys","rpoB_p.His445Gly","rpoB_p.Gln432Pro","rpoB_p.His445Gln",
            "rpoB_p.Asn438del","rpoB_p.Asp435Gly","rpoB_p.Leu430Arg","rpoB_p.Met434Val","rpoB_p.Asn437Asp",
            "rpoB_p.Asp435Ala","rpoB_p.Ala451Val","rpoB_p.Thr427Ala"
        ]
        AssowR_and_Uncertain = [
            "rpoB_p.Ser450Leu","rpoB_p.His445Tyr","rpoB_p.His445Asp","rpoB_p.Leu452Pro","rpoB_p.Asp435Val",
            "rpoB_p.Leu430Pro","rpoB_p.His445Leu","rpoB_p.Ser450Trp","rpoB_p.Asp435Tyr","rpoB_p.His445Arg",
            "rpoB_p.His445Asn","rpoB_p.Val170Phe","rpoB_p.Ser441Leu","rpoB_p.Ile491Phe","rpoB_p.Ser450Phe",
            "rpoB_p.His445Ser","rpoB_p.His445Cys","rpoB_p.His445Gly","rpoB_p.Gln432Pro","rpoB_p.His445Gln",
            "rpoB_p.Asn438del","rpoB_p.Asp435Gly","rpoB_p.Leu430Arg","rpoB_p.Thr400Ala","rpoB_p.Ile480Val",
            "rpoC_p.Leu527Val","rpoC_p.Ile491Thr","rpoC_p.Gly332Arg","rpoB_p.Asp435Ala","rpoB_p.Met434Val",
            "rpoB_p.Ser493Leu","rpoB_p.Asn437Asp","rpoB_p.Ala451Val","rpoC_p.Asp485Asn","rpoB_p.Ala286Val",
            "rpoC_p.Asp747Ala","rpoC_p.Phe452Ser","rpoB_p.Pro454His","rpoB_p.Thr427Ala"
        ]



    class EMB_feat:
        AssowR_predominant = [
            "embB_p.Met306Val","embB_p.Met306Ile","embB_p.Gln497Arg","n.4243221C>T"
        ]
        AssowR = [
            "embB_p.Met306Val","embB_p.Met306Ile","embB_p.Gln497Arg","embB_p.Gly406Ala","embB_p.Gly406Asp",
            "embB_p.Asp354Ala","n.4243221C>T","embB_p.Gln497Lys","embB_p.Tyr319Ser","embB_p.Gly406Ser",
            "embB_p.Asp328Tyr","embB_p.Tyr319Cys","embB_p.Met306Leu","embB_p.Gly406Cys"
        ]
        AssowR_and_Uncertain = [
            "embB_p.Met306Val","embB_p.Met306Ile","embB_p.Gln497Arg","embB_p.Gly406Ala","embB_p.Gly406Asp",
            "embB_p.Asp354Ala","n.4243221C>T","embB_p.Gln497Lys","embB_p.Tyr319Ser","embB_p.Gly406Ser",
            "n.4243217C>G","embB_p.Asp1024Asn","embB_p.His1002Arg","embB_p.Gln497Pro","embB_p.Tyr319Cys",
            "embB_p.Asp328Tyr","embB_p.Met306Leu","embB_p.Ser380Asn","n.4243217C>T","n.4243222C>A","embB_p.Gly406Cys","embB_p.Ala356Val"
        ]

    class ETO_feat:
        AssowR_predominant = [
            "n.1673425C>T","ethA_frameshift_variant","n.1674048G>A"
        ]
        AssowR = [
            "n.1673425C>T","ethA_frameshift_variant","n.1674048G>A","ethA_stop_gained","inhA_p.Ser94Ala","n.1673432T>C"
        ]
        AssowR_and_Uncertain = [
            "n.1673425C>T","ethA_frameshift_variant","n.1674048G>A","ethA_stop_gained","inhA_p.Ser94Ala","n.1673432T>C",
            "inhA_p.Ile21Thr","n.1673423G>T","inhA_p.Ile194Thr"
        ]

    class MFX_feat:
        AssowR_predominant = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro"
        ]
        AssowR = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro","gyrB_p.Asp461Asn","gyrB_p.Glu501Asp","gyrB_p.Ala504Val","gyrB_p.Asn499Asp"
        ]
        AssowR_and_Uncertain = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro","gyrA_p.Asp89Asn","gyrB_p.Asp461Asn","gyrB_p.Glu501Asp","gyrB_p.Thr500Asn",
            "gyrB_p.Asn499Thr","gyrB_p.Ser447Phe","gyrB_p.Ala504Val","gyrB_p.Asn499Lys","gyrB_p.Asp461His","gyrB_p.Asn499Asp",
            "gyrA_p.Ala74Ser"
        ]

    class LFX_feat:
        AssowR_predominant = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro"
        ]
        AssowR = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro","gyrB_p.Asp461Asn","gyrB_p.Glu501Asp","gyrB_p.Ala504Val","gyrB_p.Asn499Asp"
        ]
        AssowR_and_Uncertain = [
            "gyrA_p.Ala90Val","gyrA_p.Ser91Pro","gyrA_p.Asp89Asn","gyrB_p.Asp461Asn","gyrB_p.Glu501Asp","gyrB_p.Thr500Asn",
            "gyrB_p.Asn499Thr","gyrB_p.Ser447Phe","gyrB_p.Ala504Val","gyrB_p.Asn499Lys","gyrB_p.Asp461His","gyrB_p.Asn499Asp",
            "gyrA_p.Ala74Ser"
        ]

    class AMK_feat:
        AssowR_predominant = [
            "n.1473246A>G","n.1473329G>T","n.2715346G>A"
        ]
        AssowR = [
            "n.1473246A>G","n.1473329G>T","n.2715346G>A"
        ]
        AssowR_and_Uncertain = [
            "n.1473246A>G","n.1473329G>T","n.2715346G>A"
        ]

    class KAN_feat:
        AssowR_predominant = [
            "n.1473246A>G","n.2715342C>T"
        ]
        AssowR = [
            "n.1473246A>G","n.1473329G>T","n.1473369C>A","n.2715342C>T","n.2715346G>A"
        ]
        AssowR_and_Uncertain = [
            "n.1473246A>G","n.1473329G>T","n.1473369C>A","n.2715342C>T","n.2715346G>A"
        ]

    drug_mapping = {
        'Isoniazid': INH_feat, 'Rifampicin': RIF_feat, 'Ethambutol': EMB_feat, 'Ethionamide': ETO_feat,
        'Moxifloxacin': MFX_feat, 'Levofloxacin': LFX_feat, 'Amikacin': AMK_feat, 'Kanamycin': KAN_feat
    }
    
    @classmethod
    def get_features(cls, drug_name, strategy):
        print(drug_name)
        if drug_name not in cls.drug_mapping:
            raise ValueError(f"{drug_name} error")
        feat_class = cls.drug_mapping[drug_name]
        if not hasattr(feat_class, strategy):
            raise ValueError(f"{drug_name} : {strategy} error")
        return getattr(feat_class, strategy)
    
    @classmethod
    def get_all_features_for_drug(cls, drug_name):
        feat_class = cls.drug_mapping[drug_name]
        all_feats = []
        for s in ['AssowR_predominant', 'AssowR', 'AssowR_and_Uncertain']:
            all_feats.extend(getattr(feat_class, s))
        return sorted(list(set(all_feats)))

class SNPClassifierPipeline:
    def __init__(self, snp_matrix_file, outdir_base, random_state=42):
        self.snp_matrix_file = snp_matrix_file
        self.outdir_base = os.path.abspath(outdir_base)
        self.random_state = random_state
        self.snp_data = None
        self.snp_names = None

        self.all_fold_metrics = {}
        self.all_summary_results = {}
        self.all_ttest_results = {}
        self.all_feature_importance = {}

        os.makedirs(self.outdir_base, exist_ok=True)
        self._load_snp_matrix()

    def _load_snp_matrix(self):
        print(f"\nloading data: {self.snp_matrix_file}")
        self.snp_data = pd.read_csv(self.snp_matrix_file, index_col=0, sep='\t')
        self.snp_names = self.snp_data.index.tolist()
        print(f"SNP matrix shape: {self.snp_data.shape[0]} SNPs × {self.snp_data.shape[1]} samples")

    def _load_group_file(self, group_file):
        print(f"\nloading groups: {group_file}")
        group_df = pd.read_csv(group_file, sep='\t', names=['sample', 'group'])
        group_df['group'] = group_df['group'].str.strip().str.lower()
        group_df['y'] = group_df['group'].apply(lambda x: 0 if 'sensi' in x else 1)
        return group_df.set_index('sample')

    def _filter_data(self, drug, strategy, group_df):
        selected = DrugFeatures.get_features(drug, strategy)
        exist = [f for f in selected if f in self.snp_names]
        X = self.snp_data.loc[exist].T
        common = X.index.intersection(group_df.index)
        X = X.loc[common]
        y = group_df.loc[common, 'y']
        print(f"Samples: {len(y)} | Sens:{sum(y==0)} Res{sum(y==1)} | Features:{len(exist)}")
        return X, y, exist

    def _cross_validate(self, X, y, cv_folds=10):
        skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
        fold_metrics = {'AUPRC':[],'F1':[],'Precision':[],'Recall':[]}
        imp_list = []

        for fold, (tr, va) in enumerate(skf.split(X, y), 1):
            model = BalancedRandomForestClassifier(n_estimators=200, random_state=self.random_state, n_jobs=-1)
            model.fit(X.iloc[tr], y.iloc[tr])
            imp_list.append(model.feature_importances_)
            y_pred = model.predict(X.iloc[va])
            y_prob = model.predict_proba(X.iloc[va])[:, 1]

            fold_metrics['AUPRC'].append(average_precision_score(y.iloc[va], y_prob))
            fold_metrics['F1'].append(f1_score(y.iloc[va], y_pred, zero_division=0))
            fold_metrics['Precision'].append(precision_score(y.iloc[va], y_pred, zero_division=0))
            fold_metrics['Recall'].append(recall_score(y.iloc[va], y_pred, zero_division=0))

        summary = {}
        for k, v in fold_metrics.items():
            mean_val = np.mean(v)
            std_val = np.std(v)
            summary[k] = {'mean±std': f"{mean_val:.4f}±{std_val:.4f}", 'folds': v}

        avg_imp = dict(zip(X.columns, np.mean(imp_list, axis=0)))
        return summary, fold_metrics, avg_imp

    def _train_holdout(self, X, y, drug, strategy, outdir):
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=self.random_state,stratify=y
            
        )
        model = BalancedRandomForestClassifier(n_estimators=200, random_state=self.random_state, n_jobs=-1)
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test)

        ap = average_precision_score(y_test, y_prob)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        pre = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred, zero_division=0)

        save_path = os.path.join(outdir, f"{drug}_{strategy}_pr_curve.pdf")
        self.plot_pr_curve(y_test, y_prob, save_path)

        return {
            'AUPRC': f"{ap:.4f}",
            'F1': f"{f1:.4f}",
            'Precision': f"{pre:.4f}",
            'Recall': f"{rec:.4f}"
        }

    def plot_pr_curve(self, y_true, y_pred_proba, save_path):

        def center_dot(x):
            return str(round(x,1)).replace(".", "·")

        precision, recall, thresholds = precision_recall_curve(y_true, y_pred_proba)
        avg_precision = average_precision_score(y_true, y_pred_proba)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
        ix = np.argmax(f1_scores[:-1])

        plt.figure(figsize=(10, 8))
        plt.plot(recall, precision, color='darkgreen', lw=2,
                 label=f'PR curve (AUPRC = {avg_precision:.4f})')
        plt.plot(recall[ix], precision[ix], 'ro', markersize=5,
                 label=f'Best threshold (F1={f1_scores[ix]:.3f}, Precision={precision[ix]:.3f}, Recall={recall[ix]:.3f})')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        ax = plt.gca()
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: center_dot(x)))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: center_dot(x)))
        plt.xlabel('Recall', fontsize=9)
        plt.ylabel('Precision', fontsize=9)
        plt.title('Precision-recall - balanced random forest', fontsize=9)
        plt.legend(loc="upper right", fontsize=9)
        plt.tick_params(axis='both', labelsize=9) 

        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

    def run_drug_analysis(self, drug_name, group_file, cv_folds=10):
        print("\n" + "="*80)
        print(f"analyzing: {drug_name}")
        print("="*80)

        drug_dir = os.path.join(self.outdir_base, drug_name)
        os.makedirs(drug_dir, exist_ok=True)
        group_df = self._load_group_file(group_file)
        strategies = ['AssowR_predominant', 'AssowR', 'AssowR_and_Uncertain']

        cv_results = {}
        fold_results = {}
        holdout_results = {}
        feat_imp_results = {}

        for strategy in strategies:
            print(f"\nstrategy: {strategy}")
            X, y, _ = self._filter_data(drug_name, strategy, group_df)
            cv, fold, imp = self._cross_validate(X, y, cv_folds)
            hold = self._train_holdout(X, y, drug_name, strategy, drug_dir)

            cv_results[strategy] = cv
            fold_results[strategy] = fold
            holdout_results[strategy] = hold
            feat_imp_results[strategy] = imp

            
            print(f"Overview: AUPRC {hold['AUPRC']} | {cv['AUPRC']['mean±std']} | F1 {cv['F1']['mean±std']}")

        self.all_fold_metrics[drug_name] = fold_results
        self.all_summary_results[drug_name] = (cv_results, holdout_results)
        self.all_feature_importance[drug_name] = feat_imp_results

        ttest = self._paired_ttest(fold_results)
        self.all_ttest_results[drug_name] = ttest


    def _paired_ttest(self, fold_metrics):
        print("\nPair t test (with FDR BH correction + NA handling)...")
        metrics = ['AUPRC', 'F1', 'Recall', 'Precision']
        strategies = [s for s in ['AssowR_predominant','AssowR','AssowR_and_Uncertain'] if fold_metrics[s]]
        pairs = list(itertools.combinations(strategies, 2)) 
        pairs.remove(('AssowR_predominant', 'AssowR_and_Uncertain'))
        
        rows = []
        for m in metrics:
            for a, b in pairs:
                d1 = fold_metrics[a][m]
                d2 = fold_metrics[b][m]
                
                try:
                    t_stat, p_val = ttest_rel(d1, d2)
                except:
                    t_stat, p_val = np.nan, np.nan
                
                rows.append([
                    m, a, b, t_stat, p_val,
                    f"{np.mean(d1):.4f}±{np.std(d1):.4f}",
                    f"{np.mean(d2):.4f}±{np.std(d2):.4f}"
                ])

        df = pd.DataFrame(rows, columns=['metric','group1','group2','t','p','s1_mean±std','s2_mean±std'])

        valid_p = df['p'].dropna()
        p_corrected = np.full(len(df), np.nan) 
        
        if len(valid_p) > 0:
            rejected, corrected, _, _ = multipletests(valid_p, method='fdr_bh')
            p_corrected[valid_p.index] = corrected

        df['p_corrected'] = p_corrected


        def get_sig(p):
            if pd.isna(p):
                return 'ns'
            elif p < 0.001:
                return '***'
            elif p < 0.01:
                return '**'
            elif p < 0.05:
                return '*'
            else:
                return 'ns'

        df['significance'] = df['p_corrected'].apply(get_sig)
        return df


    def save_all_results(self):
        print("\nSaving all results...")

        with pd.ExcelWriter(os.path.join(self.outdir_base, 'summary_results.xlsx'), engine='openpyxl') as w:
            for drug, (cv, hold) in self.all_summary_results.items():
                rows = []
                data = {
                    "Metrics": ["AUPRC", "F1 Score", "Precision", "Recall"],
                    "Model(AssowR_predominant)": [
                        hold['AssowR_predominant']['AUPRC'],
                        hold['AssowR_predominant']['F1'],
                        hold['AssowR_predominant']['Precision'],
                        hold['AssowR_predominant']['Recall']
                    ],
                    "Cross-validation(AssowR_predominant)": [
                        cv['AssowR_predominant']['AUPRC']['mean±std'],
                        cv['AssowR_predominant']['F1']['mean±std'],
                        cv['AssowR_predominant']['Precision']['mean±std'],
                        cv['AssowR_predominant']['Recall']['mean±std']
                    ],
                    "Model(AssowR)": [
                        hold['AssowR']['AUPRC'],
                        hold['AssowR']['F1'],
                        hold['AssowR']['Precision'],
                        hold['AssowR']['Recall']
                    ],
                    "Cross-validation(AssowR)": [
                        cv['AssowR']['AUPRC']['mean±std'],
                        cv['AssowR']['F1']['mean±std'],
                        cv['AssowR']['Precision']['mean±std'],
                        cv['AssowR']['Recall']['mean±std']
                    ],
                    "Model(AssowR_and_Uncertain)": [
                        hold['AssowR_and_Uncertain']['AUPRC'],
                        hold['AssowR_and_Uncertain']['F1'],
                        hold['AssowR_and_Uncertain']['Precision'],
                        hold['AssowR_and_Uncertain']['Recall']
                    ],
                    "Cross-validation(AssowR_and_Uncertain)": [
                        cv['AssowR_and_Uncertain']['AUPRC']['mean±std'],
                        cv['AssowR_and_Uncertain']['F1']['mean±std'],
                        cv['AssowR_and_Uncertain']['Precision']['mean±std'],
                        cv['AssowR_and_Uncertain']['Recall']['mean±std']
                    ]
                }
                df_summary = pd.DataFrame(data)
                df_summary.to_excel(w, sheet_name=drug, index=False)

        with pd.ExcelWriter(os.path.join(self.outdir_base, 'fold_metrics.xlsx'), engine='openpyxl') as w:
            for drug, fold in self.all_fold_metrics.items():
                rows = []
                for s in ['AssowR_predominant','AssowR','AssowR_and_Uncertain']:
                    for i in range(len(fold[s]['F1'])):
                        rows.append([s, i+1,
                                     fold[s]['AUPRC'][i],
                                     fold[s]['F1'][i],
                                     fold[s]['Precision'][i],
                                     fold[s]['Recall'][i]])
                pd.DataFrame(rows, columns=['strategy','fold','AUPRC','F1','Precision','Recall']).to_excel(w, sheet_name=drug, index=False)

        with pd.ExcelWriter(os.path.join(self.outdir_base, 't_test_results.xlsx'), engine='openpyxl') as w:
            for drug, df in self.all_ttest_results.items():
                df.to_excel(w, sheet_name=drug, index=False)

        with pd.ExcelWriter(os.path.join(self.outdir_base, 'feature_importance.xlsx'), engine='openpyxl') as w:
            for drug, imp_dict in self.all_feature_importance.items():
                all_feats = DrugFeatures.get_all_features_for_drug(drug)
                df = pd.DataFrame({'Feature': all_feats})
                for s in ['AssowR_predominant','AssowR','AssowR_and_Uncertain']:
                    df[s] = df['Feature'].map(imp_dict.get(s, {})).fillna('-')
                    df[s] = df[s].apply(lambda x: f"{x:.6f}" if isinstance(x, (int, float)) else x)
                
                df['sort_col'] = df['AssowR_and_Uncertain'].replace('-', '0').astype(float)
                df_sorted = df.sort_values('sort_col', ascending=False).drop('sort_col', axis=1)
                
                df_sorted.to_excel(w, sheet_name=drug, index=False)

        print("OK!")

def main():
    seed = int(sys.argv[1])

    pipeline = SNPClassifierPipeline(snp_matrix_file="AllSamples.FullSnps.matrix",outdir_base=f"./Seed_{seed}",random_state=seed)

    pipeline.run_drug_analysis("Isoniazid", "INH.group", cv_folds=10)
    pipeline.run_drug_analysis("Rifampicin", "RIF.group", cv_folds=10)
    pipeline.run_drug_analysis("Ethambutol", "EMB.group", cv_folds=10)
    pipeline.run_drug_analysis("Ethionamide", "ETO.group", cv_folds=10)
    pipeline.run_drug_analysis("Moxifloxacin", "MFX.group", cv_folds=10)
    pipeline.run_drug_analysis("Levofloxacin", "LFX.group", cv_folds=10)
    pipeline.run_drug_analysis("Kanamycin", "KAN.group", cv_folds=10)
    pipeline.run_drug_analysis("Amikacin", "AMK.group", cv_folds=10)

    pipeline.save_all_results()

    print("\nAll done!")

if __name__ == "__main__":
    main()
