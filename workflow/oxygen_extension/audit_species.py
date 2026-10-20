#!/usr/bin/env python3
"""Inventory direct species and test exploratory age associations, without inferring oxygen labels.
Run with the baseline environment: python workflow/oxygen_extension/audit_species.py
BH correction includes all nonconstant direct-species tests within each cohort/denominator.
"""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent/'outputs'
OUT.mkdir(exist_ok=True)
records=[]; manifest=[]
for cohort,filename in [('qc3','qc3_duplicates_cluster1_and_minus1_removed.csv'),('qc4','qc4_duplicates_cluster_1_minus_1_removed_age_ge_1.csv')]:
 path=ROOT/'inputs'/filename
 d=pd.read_csv(path)
 assert d.sample_id.notna().all() and not d.sample_id.duplicated().any()
 cols=[c for c in d if c.split('|')[-1].startswith('s__')]
 x=d[cols].apply(pd.to_numeric,errors='raise')
 assert x.notna().all().all() and (x>=0).all().all()
 total=x.sum(axis=1)
 assert (total>0).all()
 cov=pd.DataFrame({'sample_id':d.sample_id,'age':d.age,'weight':d.weight,'direct_species_total':total})
 cov.to_csv(OUT/f'{cohort}_denominator_audit.csv',index=False)
 rr,pp=spearmanr(d.age,total,nan_policy='omit')
 manifest.append({'cohort':cohort,'samples':len(d),'direct_species_features':len(cols),'species_total_mean':total.mean(),'species_total_min':total.min(),'species_total_age_rho':rr,'species_total_age_p':pp,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 for denom,values in [('original_abundance',x),('within_direct_species',x.div(total,axis=0))]:
  rows=[]
  for c in cols:
   y=values[c]; valid=d.age.notna() & y.notna()
   rho,p=spearmanr(d.loc[valid,'age'],y[valid]) if y[valid].nunique()>1 else (np.nan,np.nan)
   rows.append({'cohort':cohort,'denominator':denom,'clade':c,'species':c.split('s__')[-1], 'n':int(valid.sum()),'prevalence':float((x[c]>0).mean()),'mean_abundance':y.mean(),'rho':rho,'p':p})
  f=pd.DataFrame(rows); mask=f.p.notna(); f.loc[mask,'q_bh']=multipletests(f.loc[mask,'p'],method='fdr_bh')[1]; records.append(f)
 # Both direct genus and sum of species: check unresolved genus abundance.
 neis=[c for c in cols if '|g__Neisseria|' in c]
 g=[c for c in d if c.endswith('|g__Neisseria')]
 z=cov[['sample_id','age','weight']].copy(); z['neisseria_species_sum']=x[neis].sum(axis=1)
 if g: z['neisseria_direct_genus']=d[g].sum(axis=1); z['genus_minus_species']=z.neisseria_direct_genus-z.neisseria_species_sum
 z.to_csv(OUT/f'{cohort}_neisseria_abundance.csv',index=False)
 if cohort=='qc3':
  inv=pd.DataFrame({'clade':cols,'species':[c.split('s__')[-1] for c in cols],'mean_original_abundance':x.mean().values,'prevalence':(x>0).mean().values})
  inv=inv.sort_values('mean_original_abundance',ascending=False); inv['cumulative_share_of_mean_species_abundance']=inv.mean_original_abundance.cumsum()/inv.mean_original_abundance.sum()
  inv['oxygen_category']='unknown'; inv['evidence_type']='unreviewed'; inv['source_url']=''; inv['evidence_note']=''; inv['periodontal_association']='not_reviewed'
  inv.to_csv(OUT/'species_annotation_review_queue.csv',index=False)
  fig,axs=plt.subplots(3,3,figsize=(12,9),sharex=True)
  for ax,c in zip(axs.flat,neis):
   ax.scatter(d.age,x[c]*100,s=5,alpha=.25); ax.set_title(c.split('s__')[-1].replace('_',' '),fontsize=9); ax.set_ylabel('Original abundance (%)'); ax.set_xlabel('Age (years)')
  fig.tight_layout(); fig.savefig(OUT/'neisseria_species_age.png',dpi=180); plt.close(fig)
f=pd.concat(records,ignore_index=True); f.to_csv(OUT/'all_species_age_associations.csv',index=False)
f[f.clade.str.contains(r'\|g__Neisseria\|')].to_csv(OUT/'neisseria_age_associations.csv',index=False)
(OUT/'input_audit.json').write_text(json.dumps(manifest,indent=2)+'\n')
(OUT/'environment.txt').write_text(subprocess.check_output([__import__('sys').executable,'-m','pip','freeze'],text=True))
print(json.dumps(manifest,indent=2))
print(f[(f.cohort=='qc3')&(f.denominator=='original_abundance')&f.clade.str.contains(r'\|g__Neisseria\|')][['species','mean_abundance','prevalence','rho','q_bh']].to_string(index=False))
