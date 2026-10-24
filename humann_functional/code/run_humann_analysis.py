#!/usr/bin/env python3
"""Exploratory community MetaCyc age trends and technical QC of HUMAnN tables.
This pipeline does NOT infer periodontal disease status from age or taxonomy.
Only unstratified community-level pathways enter the per-pathway tests.
"""
from pathlib import Path
import argparse, json, re, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

PTH=Path(__file__).resolve().parent
DEFAULT='/mnt/data/'
parser=argparse.ArgumentParser()
parser.add_argument('--abundance',default=DEFAULT+'qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv')
parser.add_argument('--coverage',default=DEFAULT+'qc1_pathcoverage_n_reads_ge_2_mil.tsv')
parser.add_argument('--output',default=str(PTH))
a=parser.parse_args()
out=Path(a.output);(out/'results').mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(parents=True,exist_ok=True)
A=pd.read_csv(a.abundance,sep='\t');C=pd.read_csv(a.coverage,sep='\t')
for d in (A,C): d.rename(columns={d.columns[0]:'sample_id'},inplace=True)
assert A['sample_id'].is_unique and C['sample_id'].is_unique
assert set(A.sample_id)==set(C.sample_id)
C=C.set_index('sample_id').loc[A.sample_id].reset_index()
assert np.allclose(A.age,C.age,equal_nan=True)
assert np.allclose(A.weight,C.weight,equal_nan=True)
met=['sample_id','age','weight']
cols=[x for x in A.columns if x not in met and '|' not in x and not x.startswith('UN')]
covcols=[x for x in C.columns if x not in met and '|' not in x and not x.startswith('UN')]
assert set(cols)==set(covcols),('not common',len(set(cols)-set(covcols)))
assert not A[cols].isna().any().any()
assert not C[cols].isna().any().any()
ids=A.sample_id
X=A[cols].astype(float)
XC=C[cols].astype(float)
age=A.age.to_numpy(float)
weight=A.weight.to_numpy(float)
prev=(X>0).mean(axis=0)
use=prev[prev>=0.10].index.tolist()
U=X[use]
cu=C[use]
print('Read',len(A),'samples',len(cols),'pathways','tested',len(use))

# QC information and sample manifests.
sample=pd.DataFrame({'sample_id':ids,'age':age,'weight':weight,'age_ge_1':age>=1,'weight_recorded':np.isfinite(weight),
  'pathways_detected':(X>0).sum(axis=1), 'pathways_coverage_ge_0.5':(XC>=.5).sum(axis=1),
  'community_pathway_cpm_sum':X.sum(axis=1)})
sample.to_csv(out/'results/sample_manifest_695.csv',index=False)
qc_summary={'n_samples_abundance':len(A),'n_samples_coverage':len(C),'n_unique_ids':ids.nunique(),
 'n_with_age':int(np.isfinite(age).sum()),'n_with_weight':int(np.isfinite(weight).sum()),
 'n_age_below_1':int((age<1).sum()),'n_age_ge_1':int((age>=1).sum()),
 'n_unstratified_pathways':len(cols),'n_tested_prevalence_ge_10pct':len(use),
 'n_unstratified_coverage_features_including_unmapped_unintegrated':len([x for x in C if x not in met and '|' not in x]),
 'max_abs_cpm_sum_minus_1m':float((sample.community_pathway_cpm_sum-1e6).abs().max()),
 'coverage_out_of_range_count':int(((XC<0)|(XC>1)).sum().sum()),
 'missing_pd_severity':True,'notes':'No read-count manifest provided; 695 rows cannot be reconciled to reported n=625.'}
(out/'results/qc_summary.json').write_text(json.dumps(qc_summary,indent=2))

# Age comparisons: community CPM; log-ratio sensitivity (CLR with pseudocount 1 CPM)
# and restriction to adult dogs (>=1 year). These methods are complementary, not definitive.
cp=np.log(U.to_numpy()+1.0)
clr=cp-cp.mean(axis=1,keepdims=True)
clr_df=pd.DataFrame(clr,columns=use)
adult=age>=1
wtmask=np.isfinite(weight)

def corr_columns(mat,xx,mask):
    arr=mat.loc[mask].to_numpy(dtype=float)
    val=xx[mask]
    # rank correlation through scipy per column, ties retained
    if np.isnan(arr).any():raise ValueError('NaN in features')
    rho=np.empty(arr.shape[1]);p=np.empty(arr.shape[1]);
    for i in range(arr.shape[1]):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r,pv=stats.spearmanr(val,arr[:,i]);
        rho[i]=r if np.isfinite(r) else np.nan
        p[i]=pv if np.isfinite(pv) else 1.0
    q=multipletests(p,method='fdr_bh')[1]
    return rho,p,q
mask=np.ones(len(age),dtype=bool)
rr,pp,qq=corr_columns(U,age,mask)
rc,pc,qc=corr_columns(clr_df,age,mask)
ra,pa,qa=corr_columns(U,age,adult)
rac,pac,qac=corr_columns(clr_df,age,adult)

# Robust HC3 OLS for age + log(weight), complete records; model is exploratory linear model.
model_agecoef=[];model_p=[];model_coef_clr=[];model_p_clr=[]
D=sm.add_constant(np.column_stack([age[wtmask],np.log1p(weight[wtmask])]))
for i,col in enumerate(use):
    y=np.log1p(U.loc[wtmask,col].to_numpy(dtype=float))
    fit=sm.OLS(y,D).fit(cov_type='HC3')
    model_agecoef.append(fit.params[1]);model_p.append(fit.pvalues[1])
    yc=clr_df.loc[wtmask,col].to_numpy(dtype=float)
    f2=sm.OLS(yc,D).fit(cov_type='HC3')
    model_coef_clr.append(f2.params[1]);model_p_clr.append(f2.pvalues[1])
model_q=multipletests(model_p,method='fdr_bh')[1]
model_q_clr=multipletests(model_p_clr,method='fdr_bh')[1]

# high-low age contrasts (per-group CPM, no log-ratio or confidence intervals here).
group = pd.cut(age,[-np.inf,1,3,7,np.inf],labels=['<1','1–<3','3–<7','>=7'],right=False)
group_n={str(k):int((group==k).sum()) for k in group.categories}
medians=U.groupby(group,observed=True).median()
means=U.groupby(group,observed=True).mean()
medians.T.reset_index(names='pathway').to_csv(out/'results/pathway_age_group_medians_cpm.csv',index=False)
means.T.reset_index(names='pathway').to_csv(out/'results/pathway_age_group_means_cpm.csv',index=False)

R=pd.DataFrame({'pathway':use,'metacyc_id':[s.split(':',1)[0] for s in use],
 'prevalence':prev.reindex(use).values,'rho_cpm':rr,'p_cpm':pp,'q_cpm':qq,
 'rho_clr':rc,'p_clr':pc,'q_clr':qc,
 'rho_cpm_age_ge_1':ra,'q_cpm_age_ge_1':qa,
 'rho_clr_age_ge_1':rac,'q_clr_age_ge_1':qac,
 'weight_adjusted_age_coef_log1p_cpm':model_agecoef,'weight_adjusted_p_cpm':model_p,'weight_adjusted_q_cpm':model_q,
 'weight_adjusted_age_coef_clr':model_coef_clr,'weight_adjusted_p_clr':model_p_clr,'weight_adjusted_q_clr':model_q_clr,
 'median_cpm_under_1':medians.loc['<1'].reindex(use).values,
 'median_cpm_age_7plus':medians.loc['>=7'].reindex(use).values,
 'mean_cpm_under_1':means.loc['<1'].reindex(use).values,
 'mean_cpm_age_7plus':means.loc['>=7'].reindex(use).values,
 'coverage_mean':cu.mean().values,'coverage_ge_0.5_fraction':(cu>=.5).mean().values,
 'coverage_ge_0.9_fraction':(cu>=.9).mean().values})
R['raw_and_clr_significant']=(R.q_cpm<.05)&(R.q_clr<.05)&(np.sign(R.rho_cpm)==np.sign(R.rho_clr))
R['also_significant_in_age_ge_1']=(R.q_cpm_age_ge_1<.05)&(np.sign(R.rho_cpm)==np.sign(R.rho_cpm_age_ge_1))
R['also_weight_adjusted']=(R.weight_adjusted_q_cpm<.05)&(np.sign(R.rho_cpm)==np.sign(R.weight_adjusted_age_coef_log1p_cpm))
R['supported_all_four']=(R.raw_and_clr_significant&R.also_significant_in_age_ge_1&R.also_weight_adjusted)
R['mean_fold_age_7plus_vs_under_1_plus1']=(R.mean_cpm_age_7plus+1)/(R.mean_cpm_under_1+1)
R.sort_values('q_cpm',inplace=True)
R.to_csv(out/'results/pathway_age_associations_sensitivity.csv',index=False)

# Keyword-only themes are provisional, and not official MetaCyc hierarchy annotations.
# Assign mutually exclusive categories for presentation but retain all matched tags.
themes=[
 ('Sulfur and sulfur-containing amino acids',r'sulfur|sulfate|sulfite|sulfide|methionine|cysteine|taurine|homocysteine|thiosulfate|seleno'),
 ('Fermentation and anaerobic energy',r'ferment|anaerob|denitrif|nitrate reduction|nitrite reduction|acetogen|methanogen'),
 ('Respiration and TCA',r'respirat|tricarboxylic|citric acid cycle|tca cycle|glyoxylate cycle|electron transfer|oxidative phosphorylation'),
 ('Cell envelope and biofilm-associated',r'peptidoglycan|lipopolysaccharide|lipooligosaccharide|lipoteichoic|teichoic|capsul|biofilm|colanic|exopolysaccharide|cell wall|cell envelope'),
 ('Vitamins and cofactors',r'folate|biotin|thiamine|riboflavin|vitamin|coenzyme|nad\b|flavin|pantothen|heme biosyn|tetrahydrofolate|cobalamin|pyridox'),
 ('Amino acid degradation',r'(amino acid|alanine|arginine|asparagine|aspartate|glutamate|glutamine|glycine|histidine|isoleucine|leucine|lysine|ornithine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine).*?(degrad|catabol|ferment)|(?:degrad|catabol).*(amino acid|alanine|arginine|asparagine|aspartate|glutamate|glutamine|glycine|histidine|isoleucine|leucine|lysine|ornithine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine)'),
 ('Carbohydrate metabolism',r'glycolys|gluconeogen|glycogen|starch|sucrose|lactose|fructose|glucose|fucose|galactose|glucuron|mannose|pentose|hexose|xylose|ribose|arabinose|glycan|polysaccharide|cellulose|pectin'),
 ('Lipids and fatty acids',r'fatty acid|lipid|phospholipid|beta-oxidation|butyrate|propionat|acetate'),
 ('Nucleotide metabolism',r'purine|pyrimidine|nucleotide|adenosine|guanine|guanosine|uracil|uridi|thymid|cytid'),
 ('Amino acid biosynthesis',r'amino acid|alanine|arginine|asparagine|aspartate|glutamate|glutamine|glycine|histidine|isoleucine|leucine|lysine|ornithine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine'),
 ('Other degradation',r'degrad|catabol|detox|salvage'),
 ('Other biosynthesis',r'biosynth|synthes'),
]
cat=[]; tags=[]; for_output=[]
for col in cols:
 nm=col.split(':',1)[-1].lower().strip()
 matched=[label for label,expr in themes if re.search(expr,nm)]
 cat.append(matched[0] if matched else 'Other / unclassified by keywords')
 tags.append('; '.join(matched) if matched else '')
 for_output.append({'pathway':col,'metacyc_id':col.split(':',1)[0], 'name':nm,
                    'keyword_primary_category':cat[-1],'all_keyword_tags':tags[-1],
                    'annotation_method':'HEURISTIC LABEL ONLY (not MetaCyc hierarchy)'})
catdf=pd.DataFrame(for_output)
catdf.to_csv(out/'results/pathway_category_keyword_audit.csv',index=False)
R=R.merge(catdf[['pathway','keyword_primary_category']],on='pathway',how='left')
R.sort_values('q_cpm').to_csv(out/'results/pathway_age_associations_sensitivity.csv',index=False)

# Group totals are disjoint because every pathway is assigned a single primary category.
cat_mat=pd.DataFrame(index=X.index)
for label in catdf.keyword_primary_category.unique():
    cat_mat[label]=X[catdf.loc[catdf.keyword_primary_category==label,'pathway'].tolist()].sum(axis=1)
cat_mat['age']=age;cat_mat['sample_id']=ids
cat_mat.to_csv(out/'results/sample_keyword_functional_categories_cpm.csv',index=False)
cat_anal=[]
for t in cat_mat.columns[:-2]:
    rho,p=stats.spearmanr(age,cat_mat[t].values)
    cat_anal.append({'category':t,'n_pathways':int((catdf.keyword_primary_category==t).sum()),'rho_age_cpm':rho,'p':p,
        'mean_cpm':cat_mat[t].mean()})
cat_anal=pd.DataFrame(cat_anal)
cat_anal['q']=multipletests(cat_anal.p,method='fdr_bh')[1]
cat_anal.sort_values('q').to_csv(out/'results/keyword_category_age_associations.csv',index=False)

# Abundance vs coverage: show prevalence and coverage observations, not a percentage of complete pathways.
covcheck=[]
for col in cols:
 z=XC[col].values; v=X[col].values
 covcheck.append({'pathway':col,'nonzero_abundance_fraction':float(np.mean(v>0)),
  'coverage_nonzero_fraction':float(np.mean(z>0)),
  'coverage_ge_0.5_fraction':float(np.mean(z>=0.5)),
  'coverage_ge_0.9_fraction':float(np.mean(z>=0.9)),
  'median_coverage_when_abundant':float(np.median(z[v>0])) if np.any(v>0) else np.nan,
  'abundance_coverage_spearman':float(stats.spearmanr(v,z).statistic) if len(np.unique(v))>1 and len(np.unique(z))>1 else np.nan})
pd.DataFrame(covcheck).to_csv(out/'results/abundance_coverage_qc.csv',index=False)

# Species attribution diagnostic for select robust pathways. Species observations NOT disease diagnosis.
# Select top 8 robust by abs rho, aiming for biological themes and at least 3 strat columns.
strat=[s for s in A.columns if '|' in s]
select=R[R.supported_all_four].sort_values('rho_cpm',key=lambda s:s.abs(),ascending=False)
selected=[]
for t in ['Sulfur and sulfur-containing amino acids','Fermentation and anaerobic energy','Cell envelope and biofilm-associated',
          'Carbohydrate metabolism','Amino acid degradation','Respiration and TCA','Vitamins and cofactors']:
    sub=select[select.keyword_primary_category==t]
    if len(sub): selected.append(sub.iloc[0].pathway)
selected+=select.pathway.head(6).tolist()
selected=list(dict.fromkeys(selected))[:10]
spp=[]
for col in selected:
 found=[v for v in strat if v.startswith(col+'|') and '|g__' in v]
 if not found:continue
 vals=A[found].to_numpy(float)
 attrib=vals.sum(axis=0)
 denom=attrib.sum()
 ratio=np.divide(vals.sum(axis=1),X[col].values,where=X[col].values>0,out=np.full(len(A),np.nan))
 for c,quantity in sorted(zip(found,attrib),key=lambda x:-x[1])[:12]:
   spp.append({'pathway':col,'taxon':c.split('|',1)[1],
     'total_stratified_cpm':quantity,'fraction_of_total_stratified_attribution':quantity/denom if denom>0 else np.nan,
     'median_stratified_over_unstratified_ratio_across_detected_samples':np.nanmedian(ratio)})
pd.DataFrame(spp).to_csv(out/'results/top_species_contributors_selected_pathways.csv',index=False)

# Figures (title-free for manuscript use; text exported to report)
plt.rcParams.update({'font.size':10, 'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(6.5,4.5))
ax.hist(age,bins=np.arange(0,18.5,.5),edgecolor='white',alpha=.9)
ax.set_xlabel('Age (years; as recorded)');ax.set_ylabel('Number of dogs'); fig.tight_layout()
fig.savefig(out/'figures/01_age_distribution.png',dpi=250);plt.close(fig)

fig,ax=plt.subplots(figsize=(7,5))
qval=np.clip(R.q_cpm.values,1e-300,1)
ax.scatter(R.rho_cpm,-np.log10(qval),s=17,alpha=.65)
ax.axhline(-np.log10(.05),color='0.4',linestyle='--',linewidth=1)
ax.set_xlabel('Spearman rho: age vs pathway CPM');ax.set_ylabel('-log10(BH-adjusted p)');fig.tight_layout()
fig.savefig(out/'figures/02_pathway_age_volcano.png',dpi=250);plt.close(fig)

fig,ax=plt.subplots(figsize=(5.2,5.2))
ax.scatter(R.rho_cpm,R.rho_clr,s=17,alpha=.6)
ax.axhline(0,color='0.6');ax.axvline(0,color='0.6');ax.plot([-1,1],[-1,1],linestyle='--',color='0.3',linewidth=.8)
ax.set(xlabel='Age correlation (raw CPM)',ylabel='Age correlation (CLR sensitivity)',xlim=(-.85,.85),ylim=(-.85,.85));fig.tight_layout()
fig.savefig(out/'figures/03_raw_vs_clr_age_correlations.png',dpi=250);plt.close(fig)

# Representative robust pathway curves: top 6 by |rho|, require >=10% prevalence.
rob=R[R.supported_all_four].sort_values('rho_cpm',key=lambda s:s.abs(),ascending=False)
chosen=list(rob.pathway.head(6));
if len(chosen)<6:chosen=list(R.pathway.head(6))
fig,axs=plt.subplots(2,3,figsize=(14,7.2),sharex=True)
for ax,col in zip(axs.flat,chosen):
 y=np.log10(X[col].values+1)
 ax.scatter(age,y,s=7,alpha=.19)
 xo=np.linspace(0,np.ceil(age.max()),140)
 smooth=sm.nonparametric.lowess(y,age,frac=.35,return_sorted=True)
 ax.plot(smooth[:,0],smooth[:,1],lw=2)
 info=R[R.pathway==col].iloc[0];short=col.split(':',1)[-1].strip()
 if len(short)>53:short=short[:50]+'…'
 ax.set_title(short,fontsize=9)
 ax.text(.02,.95,f"rho={info.rho_cpm:+.2f}; q={info.q_cpm:.1e}",ha='left',va='top',transform=ax.transAxes,fontsize=8)
 ax.set_ylabel('log10(CPM + 1)')
ax.set_xlabel('Age (years)')
for ax in axs[-1,:]: ax.set_xlabel('Age (years)')
fig.tight_layout();fig.savefig(out/'figures/04_top_robust_pathways_vs_age.png',dpi=220);plt.close(fig)

# Category stacked bars, purely a provisional keyword classification.
cat_no_meta=cat_mat.drop(columns=['age','sample_id']).copy();cat_no_meta['age_group']=group
G=cat_no_meta.groupby('age_group',observed=True).mean()
G=G.reindex(columns=G.mean(axis=0).sort_values(ascending=False).index)
ax=G.plot(kind='bar',stacked=True,figsize=(11,5.4),width=.73)
ax.set(xlabel='Dog age group (years)',ylabel='Mean functional category CPM');ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize=8,frameon=False)
plt.xticks(rotation=0);plt.tight_layout();plt.savefig(out/'figures/05_provisional_functional_categories_by_age.png',dpi=220);plt.close()

# Ordination: PCA of CLR-transformed profile, no read-depth covariate available.
transform=clr_df.to_numpy()
standard=StandardScaler(with_std=False).fit_transform(transform)
pca=PCA(n_components=3,random_state=1);z=pca.fit_transform(standard)
fig,ax=plt.subplots(figsize=(6.7,5.1));pts=ax.scatter(z[:,0],z[:,1],c=age, cmap='viridis',s=13,alpha=.8)
ax.set(xlabel=f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% explained)',ylabel=f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% explained)')
fig.colorbar(pts,ax=ax,label='Age (years)');fig.tight_layout();fig.savefig(out/'figures/06_clr_pathway_pca_by_age.png',dpi=220);plt.close(fig)
pd.DataFrame({'sample_id':ids,'age':age,'weight':weight,'PC1':z[:,0],'PC2':z[:,1],'PC3':z[:,2]}).to_csv(out/'results/pathway_clr_pca_scores.csv',index=False)

fig,ax=plt.subplots(figsize=(6.2,4.4))
ax.scatter(age,sample.pathways_detected,s=12,alpha=.4)
s=sm.nonparametric.lowess(sample.pathways_detected,age,frac=.32)
ax.plot(s[:,0],s[:,1],linewidth=2)
ax.set(xlabel='Age (years)',ylabel='Number of detected community-level pathways');fig.tight_layout()
fig.savefig(out/'figures/07_pathway_detection_vs_age.png',dpi=220);plt.close(fig)

# Publication-cautionary summary.
r_rich,p_rich=stats.spearmanr(age,sample.pathways_detected)
pca_age,pca_p=stats.spearmanr(age,z[:,0])
summary={'qc':qc_summary,'age_group_counts':group_n,'n_tested':len(R),
  'cpm_fdr_lt05':int((R.q_cpm<.05).sum()),
  'clr_fdr_lt05':int((R.q_clr<.05).sum()),
  'age_ge_1_cpm_fdr_lt05':int((R.q_cpm_age_ge_1<.05).sum()),
  'age_ge_1_clr_fdr_lt05':int((R.q_clr_age_ge_1<.05).sum()),
  'weight_adjusted_cpm_fdr_lt05':int((R.weight_adjusted_q_cpm<.05).sum()),
  'weight_adjusted_clr_fdr_lt05':int((R.weight_adjusted_q_clr<.05).sum()),
  'consistent_cpm_clr':int(R.raw_and_clr_significant.sum()),
  'consistent_cpm_adult':int(((R.q_cpm<.05)&R.also_significant_in_age_ge_1).sum()),
  'supported_all_four':int(R.supported_all_four.sum()),
  'pathway_detection_vs_age_spearman_rho':float(r_rich),
  'pathway_detection_vs_age_spearman_p':float(p_rich),
  'pca_pc1_age_spearman_rho':float(pca_age),
  'pca_pc1_age_spearman_p':float(pca_p),
  'pca_explained_variance_ratio_first_3':pca.explained_variance_ratio_.tolist(),
  'provisional_categories':cat_anal.sort_values('q').round(4).to_dict('records'),
  'top_age_associations':R[['pathway','rho_cpm','q_cpm','rho_clr','q_clr','rho_cpm_age_ge_1','q_cpm_age_ge_1','supported_all_four']].head(12).to_dict('records')}
(out/'results/analysis_summary.json').write_text(json.dumps(summary,indent=2,default=str))
print('RESULTS:',json.dumps({k:v for k,v in summary.items() if k in ['age_group_counts','n_tested','cpm_fdr_lt05','clr_fdr_lt05','age_ge_1_cpm_fdr_lt05','age_ge_1_clr_fdr_lt05','weight_adjusted_cpm_fdr_lt05','consistent_cpm_clr','supported_all_four','pathway_detection_vs_age_spearman_rho','pca_pc1_age_spearman_rho']},indent=2))
