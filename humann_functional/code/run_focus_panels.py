#!/usr/bin/env python3
"""Hypothesis-led panels connecting age-associated pathways to taxonomic narrative.
These are community metagenomic functional associations, NOT periodontal outcomes.
"""
from pathlib import Path
import os
import pandas as pd,numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
p=Path(__file__).resolve().parent
input_dir=Path(os.environ.get('HUMANN_INPUT_DIR',str(p.parent)))
A=pd.read_csv(input_dir/'qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv',sep='\t').rename(columns={'Unnamed: 0':'sample_id'})
R=pd.read_csv(p/'results/pathway_age_associations_sensitivity.csv').set_index('pathway')
C=pd.read_csv(input_dir/'qc1_pathcoverage_n_reads_ge_2_mil.tsv',sep='\t').rename(columns={'Unnamed: 0':'sample_id'})
keys=[
'PWY-5030: L-histidine degradation III',
'P164-PWY: purine nucleobases degradation I (anaerobic)',
'P163-PWY: L-lysine fermentation to acetate and butanoate',
'DENITRIFICATION-PWY: nitrate reduction I (denitrification)',
'FUCCAT-PWY: fucose degradation',
'FERMENTATION-PWY: mixed acid fermentation',
]
assert all(z in A for z in keys)
assert all(z in C for z in keys)
# Pathways hypothesized to align with older dog taxa vs younger dog taxa; don't infer all anaerobic functions move same way.
fig,axes=plt.subplots(2,3,figsize=(13.5,7.0),sharex=True)
rows=[]
for ax,k in zip(axes.flat,keys):
 y=np.log10(A[k].to_numpy(float)+1);x=A.age.values
 rho=R.loc[k,'rho_cpm'];q=R.loc[k,'q_cpm'];clrr=R.loc[k,'rho_clr'];adult=R.loc[k,'rho_cpm_age_ge_1']
 ax.scatter(x,y,s=10,alpha=.2)
 sm=lowess(y,x,frac=.38);ax.plot(sm[:,0],sm[:,1],lw=2)
 ax.set_title(k.split(':',1)[1].strip(),fontsize=10)
 ax.text(.03,.97,f'ρ={rho:+.2f}  q={q:.1e}',transform=ax.transAxes,va='top',fontsize=9)
 ax.set_ylabel('log10(CPM + 1)')
 rows.append({'pathway':k,'rho_cpm':rho,'q_cpm':q,'rho_clr':clrr,
              'rho_age_ge_1':adult,'coverage_mean':float(C[k].mean()),
              'n_strat_cols':sum(z.startswith(k+'|') for z in A.columns),
              'mean_cpm_age_under1':float(A.loc[A.age<1,k].mean()),
              'mean_cpm_age_7plus':float(A.loc[A.age>=7,k].mean())})
for ax in axes[-1,:]: ax.set_xlabel('Age (years)')
fig.tight_layout();fig.savefig(p/'figures/08_selected_metabolic_hypothesis_pathways.png',dpi=240);plt.close(fig)
pd.DataFrame(rows).to_csv(p/'results/selected_metabolic_pathways_summary.csv',index=False)

# Top contributors per selected pathway. Each contribution proportion is explicitly among assigned stratified quantities,
# NOT among total pathway abundance; output also includes fraction of overall abundance explained by species-stratified.
parts=[];sdata=[]
for pathway in keys:
 names=[c for c in A if c.startswith(pathway+'|') and '|g__' in c]
 if not names:continue
 groups=A[names].sum(axis=0).sort_values(ascending=False)
 total=groups.sum()
 total_unstrat=A[pathway].sum()
 frac_cov=total/total_unstrat if total_unstrat>0 else np.nan
 for c,ab in groups.head(8).items():
  parts.append({'pathway':pathway,'species':c.split('|',1)[1],
                'cumulative_assigned_cpm':ab, 'share_of_assigned_taxonomic_contributions':ab/total if total else np.nan,
                'overall_taxonomic_attribution_to_community_abundance':frac_cov})
 sdata.append({'pathway':pathway,'all_taxonomically_assigned_fraction_of_total_community_cpm':frac_cov,
    'n_species_columns':len(names), 'total_community_cpm':total_unstrat})
pd.DataFrame(parts).to_csv(p/'results/selected_pathways_taxonomic_contributors.csv',index=False)
pd.DataFrame(sdata).to_csv(p/'results/selected_pathways_attribution_fraction.csv',index=False)

# Within-community pathway age trends and coverage detection (side-by-side) as a compact diagnostic.
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
for kk in [keys[0],keys[1],keys[4]]:
  xx=A.age.values
  y=np.log10(A[kk].values+1)
  ys=lowess(y,xx,frac=.4)
  axs[0].plot(ys[:,0],ys[:,1],label=kk.split(':',1)[-1][:34])
  yy=C[kk].values
  cc=lowess(yy,xx,frac=.4)
  axs[1].plot(cc[:,0],cc[:,1],label=kk.split(':',1)[-1][:34])
axs[0].set(xlabel='Age (years)',ylabel='LOWESS log10(CPM+1)')
axs[1].set(xlabel='Age (years)',ylabel='LOWESS HUMAnN pathway coverage')
axs[1].legend(fontsize=8,loc='lower right',frameon=False)
fig.tight_layout();fig.savefig(p/'figures/09_selected_pathway_abundance_vs_coverage.png',dpi=230);plt.close(fig)
print(pd.DataFrame(rows).to_string(index=False))
print('\nAttribution fractions:\n',pd.DataFrame(sdata).to_string(index=False))

# Coverage-supported shortlist: exploratory rule, not a standard HUMAnN clinical cutoff.
full=R.reset_index()
short=full[(full.supported_all_four)&(full['coverage_ge_0.5_fraction']>=.30)].copy()
short=short.sort_values('q_cpm')
short.to_csv(p/'results/coverage_supported_priority_pathways.csv',index=False)

strong=['HISDEG-PWY: L-histidine degradation I',
        'PWY-7791: UMP biosynthesis III',
        'PWY-7977: L-methionine biosynthesis IV',
        'VALSYN-PWY: L-valine biosynthesis',
        'FASYN-ELONG-PWY: fatty acid elongation -- saturated',
        'PWY-7663: gondoate biosynthesis (anaerobic)']
assert all(x in A for x in strong)
fig,axs=plt.subplots(2,3,figsize=(13.5,7.2),sharex=True)
for ax,k in zip(axs.flat,strong):
 y=np.log10(A[k].values+1)
 ax.scatter(A.age,y,s=10,alpha=.18)
 w=lowess(y,A.age,frac=.37)
 ax.plot(w[:,0],w[:,1],lw=2)
 ax.set_title(k.split(':',1)[-1].strip(),fontsize=10)
 ax.text(.03,.95,f"ρ={R.loc[k,'rho_cpm']:+.2f}; ≥0.5 coverage in {R.loc[k,'coverage_ge_0.5_fraction']:.0%}",transform=ax.transAxes,va='top',fontsize=8)
 ax.set_ylabel('log10(CPM+1)')
for ax in axs[-1,:]:ax.set_xlabel('Age (years)')
fig.tight_layout();fig.savefig(p/'figures/10_coverage_supported_age_pathways.png',dpi=240);plt.close(fig)

# Figure of age-group pathway changes, standardized within feature (avoids comparing pathway CPM magnitudes).
keep=['HISDEG-PWY: L-histidine degradation I','PWY-7791: UMP biosynthesis III',
      'PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)',
      'PWY-7977: L-methionine biosynthesis IV','VALSYN-PWY: L-valine biosynthesis',
      'FASYN-ELONG-PWY: fatty acid elongation -- saturated',
      'PWY-7663: gondoate biosynthesis (anaerobic)']
group=pd.cut(A.age,[-np.inf,1,3,7,np.inf],right=False,labels=['<1','1–<3','3–<7','≥7'])
V=np.log1p(A[keep]); V=(V-V.mean())/V.std(ddof=0)
V['age_group']=group
G=V.groupby('age_group',observed=True).mean().T
fig,ax=plt.subplots(figsize=(8.1,4.4));im=ax.imshow(G.values,aspect='auto',cmap='coolwarm',vmin=-1,vmax=1)
ax.set_xticks(np.arange(4),G.columns.tolist());ax.set_yticks(np.arange(len(keep)),[z.split(':',1)[1].strip()[:58] for z in keep],fontsize=9)
cb=fig.colorbar(im,ax=ax);cb.set_label('Mean within-pathway standardized log(CPM + 1)')
fig.tight_layout();fig.savefig(p/'figures/11_age_group_supported_pathway_heatmap.png',dpi=230);plt.close(fig)

# Decompose genus contribution for coverage-supported histidine degradation I.
target='HISDEG-PWY: L-histidine degradation I'
subset=[v for v in A if v.startswith(target+'|') and '|g__' in v]
sp=A[subset].sum()
sp.index=[v.split('|',1)[1] for v in subset]
sp=sp.sort_values(ascending=False)
denom=A[target].sum()
spout=[]
for spname,quantity in sp.items():
 spout.append({'pathway':target,'species':spname,
 'total_stratified_cpm':quantity,'fraction_among_all_assigned_taxa':quantity/sp.sum(),
 'fraction_of_total_community_cpm':quantity/denom})
pd.DataFrame(spout).to_csv(p/'results/histidine_I_species_contributors.csv',index=False)
print('\nPrioritized coverage-supported pathway count:',len(short),'positive',int((short.rho_cpm>0).sum()))
