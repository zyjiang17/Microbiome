#!/usr/bin/env python3
"""Secondary weight analysis; weight unit unknown, no disease or breed covariate."""
from pathlib import Path
import os
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
input_dir=Path(os.environ.get('HUMANN_INPUT_DIR',str(p.parent)))
A=pd.read_csv(input_dir/'qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv',sep='\t')
R=pd.read_csv(p/'results/pathway_age_associations_sensitivity.csv')
use=R.pathway.tolist()
m=A.weight.notna().values
x=A.loc[m,use].to_numpy(float)
log=np.log1p(x)
clr=log-log.mean(axis=1,keepdims=True)
wt=A.loc[m,'weight'].values;age=A.loc[m,'age'].values
rs=[]
for i,c in enumerate(use):
 r,pp=stats.spearmanr(wt,x[:,i]);rc,pc=stats.spearmanr(wt,clr[:,i]);
 rs.append({'pathway':c,'rho_weight_cpm':r,'p_weight_cpm':pp,'rho_weight_clr':rc,'p_weight_clr':pc})
out=pd.DataFrame(rs)
out['q_weight_cpm']=multipletests(out.p_weight_cpm,method='fdr_bh')[1]
out['q_weight_clr']=multipletests(out.p_weight_clr,method='fdr_bh')[1]
D=sm.add_constant(np.column_stack([age,np.log1p(wt)]))
for c,mat in [('cpm',log),('clr',clr)]:
 weightcoef=[];pval=[]
 for i in range(len(use)):
  fit=sm.OLS(mat[:,i],D).fit(cov_type='HC3')
  weightcoef.append(fit.params[2]);pval.append(fit.pvalues[2])
 out[f'age_adjusted_weight_coef_{c}']=weightcoef
 out[f'age_adjusted_weight_p_{c}']=pval
 out[f'age_adjusted_weight_q_{c}']=multipletests(pval,method='fdr_bh')[1]
out=out.merge(R[['pathway','rho_cpm','q_cpm','coverage_ge_0.5_fraction']],on='pathway')
out.sort_values('q_weight_cpm').to_csv(p/'results/pathway_weight_associations_secondary.csv',index=False)
fig,ax=plt.subplots(figsize=(5.3,5.1))
ax.scatter(out.rho_cpm,out.rho_weight_cpm,s=16,alpha=.65)
ax.axhline(0,linewidth=.8,color='0.5');ax.axvline(0,linewidth=.8,color='0.5')
ax.set(xlabel='Spearman rho: pathway CPM vs age (n=695)',ylabel='Spearman rho: pathway CPM vs weight (n=626)',xlim=(-.7,.7),ylim=(-.7,.7))
fig.tight_layout();fig.savefig(p/'figures/12_age_vs_weight_pathway_correlations.png',dpi=240);plt.close(fig)
print('Raw CPM significant weight:',(out.q_weight_cpm<.05).sum(),'; CLR weight:',(out.q_weight_clr<.05).sum(),'; age-adjusted weight CPM:',(out.age_adjusted_weight_q_cpm<.05).sum(),'; CLR:',(out.age_adjusted_weight_q_clr<.05).sum())
print('age-weight spearman',stats.spearmanr(age,wt))
