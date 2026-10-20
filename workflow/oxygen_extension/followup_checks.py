#!/usr/bin/env python3
"""Audit annotation coverage, multirank aggregation, and weight-sensitive Neisseria trends."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr,rankdata
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent; W=R.parent; O=R/'outputs'
a=pd.read_csv(O/'species_annotation_review_queue.csv'); e=pd.read_csv(R/'reviewed_evidence.csv').set_index('species')
for k in ['oxygen_category','evidence_type','source_url','evidence_note','periodontal_association']:
 a[k]=a.species.map(e[k]).fillna(a[k].fillna(''))
a['confidence']=a.species.map(e.confidence).fillna('unreviewed')
a.to_csv(O/'species_annotation_current.csv',index=False)
coverage=[]; adjusted=[]; groups=[]; order_checks=[]
for cohort,filename in [('qc3','qc3_duplicates_cluster1_and_minus1_removed.csv'),('qc4','qc4_duplicates_cluster_1_minus_1_removed_age_ge_1.csv')]:
 d=pd.read_csv(W/'inputs'/filename); x=d[a.clade]; total=x.sum(axis=1)
 labels=a.set_index('clade').oxygen_category
 known=labels.index[labels!='unknown']; reviewed=x[known].sum(axis=1)
 coverage.append({'cohort':cohort,'reviewed_species_count':len(known),'mean_reviewed_original_abundance':reviewed.mean(),'mean_reviewed_fraction_of_direct_species':(reviewed/total).mean(),'coverage_age_rho':spearmanr(d.age,reviewed/total).statistic})
 for name,v in [('direct_genus',d[[c for c in d if c.endswith('|g__Neisseria')]].sum(axis=1)),('sum_direct_species',x[[c for c in x if '|g__Neisseria|' in c]].sum(axis=1))]:
  rho,p=spearmanr(d.age,v); groups.append({'cohort':cohort,'feature':name,'rho':rho,'p':p})
 mask=d.age.notna() & d.weight.notna(); dm=d.loc[mask]
 for c in [c for c in x if '|g__Neisseria|' in c]:
  if x.loc[mask,c].nunique()<2: continue
  y=rankdata(x.loc[mask,c]); age=rankdata(dm.age); wt=rankdata(dm.weight)
  unadj=sm.OLS(y,sm.add_constant(age)).fit(cov_type='HC3'); adj=sm.OLS(y,sm.add_constant(np.column_stack([age,wt]))).fit(cov_type='HC3')
  adjusted.append({'cohort':cohort,'species':c.split('s__')[-1],'n_same_complete_cases':len(dm),'rank_age_coefficient_unadjusted':unadj.params[1],'p_unadjusted':unadj.pvalues[1],'rank_age_coefficient_weight_adjusted':adj.params[1],'p_adjusted':adj.pvalues[1]})
 # Compare frozen rank-summing formula with a direct-order calculation.
 features=[c for c in d if c not in ['sample_id','age','weight']]
 def order(c): return next((v[3:] for v in c.split('|') if v.startswith('o__')),None)
 mult=[c for c in features if order(c)]; direct=[c for c in features if c.split('|')[-1].startswith('o__')]
 def collapse(cols):
  z=d[cols].T.groupby({c:order(c) for c in cols}).sum().T
  return z.div(z.sum(axis=1),axis=0)
 old,new=collapse(mult).align(collapse(direct),join='outer',axis=1,fill_value=0)
 delta=old-new
 order_checks.append({'cohort':cohort,'multirank_input_columns':len(mult),'direct_order_columns':len(direct),'mean_total_variation_distance':float((abs(delta).sum(axis=1)/2).mean()),'max_sample_total_variation_distance':float((abs(delta).sum(axis=1)/2).max()),'largest_absolute_age_group_mean_difference':float(delta.groupby(pd.cut(d.age,[0,2,7,np.inf],include_lowest=True),observed=True).mean().abs().max().max())})
 new.assign(sample_id=d.sample_id,age=d.age).to_csv(O/f'{cohort}_direct_order_composition.csv',index=False)
pd.DataFrame(coverage).to_csv(O/'annotation_coverage.csv',index=False)
f=pd.DataFrame(adjusted)
for co in f.cohort.unique():
 m=f.cohort==co; f.loc[m,'q_bh_adjusted_within_neisseria']=multipletests(f.loc[m,'p_adjusted'],method='fdr_bh')[1]
f.to_csv(O/'neisseria_weight_sensitivity.csv',index=False)
pd.DataFrame(groups).to_csv(O/'neisseria_genus_checks.csv',index=False)
pd.DataFrame(order_checks).to_csv(O/'order_aggregation_audit.csv',index=False)
print(pd.DataFrame(coverage).to_string(index=False)); print(pd.DataFrame(order_checks).to_string(index=False)); print(pd.DataFrame(groups).to_string(index=False))
