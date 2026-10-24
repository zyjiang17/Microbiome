#!/usr/bin/env python3
"""Compare HUMAnN sample IDs against pre-existing MetaPhlAn manuscript QC cohorts.
This is NOT data harmonization; never guess/collapse nonidentical sample names.
"""
from pathlib import Path
import argparse
import pandas as pd
p=argparse.ArgumentParser()
p.add_argument('--input-dir',required=True)
p.add_argument('--taxa-dir',required=True)
p.add_argument('--out-dir',required=True)
a=p.parse_args()
out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
h=pd.read_csv(Path(a.input_dir)/'qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv',sep='\t',usecols=[0])
H=set(h.iloc[:,0].astype(str).str.strip())
files={
 'QC1':'qc1_duplicates_removed_only.csv',
 'QC2':'qc2_duplicates_removed_prevalence_1pct.csv',
 'QC3':'qc3_duplicates_cluster1_and_minus1_removed.csv',
 'QC4':'qc4_duplicates_cluster_1_minus_1_removed_age_ge_1.csv',
}
results=[]
for cohort,fname in files.items():
 path=Path(a.taxa_dir)/fname
 if not path.exists():
  results.append({'cohort':cohort,'file':fname,'status':'missing','taxonomic_n':None,'humann_n':len(H),'intersection':None})
  continue
 t=pd.read_csv(path,usecols=[0]); T=set(t.iloc[:,0].astype(str).str.strip())
 results.append({'cohort':cohort,'file':fname,'status':'OK','taxonomic_n':len(T),'humann_n':len(H),'intersection':len(H&T),'humann_only':len(H-T),'taxonomic_only':len(T-H)})
 pd.DataFrame({'sample_id':sorted(H|T),'in_humann':[x in H for x in sorted(H|T)],'in_taxonomy':[x in T for x in sorted(H|T)]}).to_csv(out/f'{cohort.lower()}_sample_overlap.csv',index=False)
df=pd.DataFrame(results);df.to_csv(out/'taxonomy_humann_overlap_summary.csv',index=False)
print(df.to_string(index=False))
if 'intersection' in df and not df['intersection'].fillna(0).any():
 print('WARNING: Zero exact-ID overlap. Sample naming may differ; do not merge by approximate suffix without validating mapping.')
