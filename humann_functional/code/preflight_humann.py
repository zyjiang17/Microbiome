#!/usr/bin/env python3
"""Check inputs before cluster analysis; never silently drop samples."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--input-dir',required=True);p.add_argument('--out',required=True);a=p.parse_args()
x=Path(a.input_dir)
names=['qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv','qc1_pathcoverage_n_reads_ge_2_mil.tsv']
tables=[pd.read_csv(x/f,sep='\t') for f in names]
for t in tables:t.rename(columns={t.columns[0]:'sample_id'},inplace=True)
A,C=tables
if A.sample_id.duplicated().any() or C.sample_id.duplicated().any():raise ValueError('Duplicate sample IDs')
if set(A.sample_id)!=set(C.sample_id):raise ValueError('Abundance and coverage sample IDs differ')
meta=['sample_id','age','weight']
for m in ('age','weight'):
 if m not in A or m not in C:raise ValueError('Missing '+m)
base=[v for v in A.columns if v not in meta and '|' not in v and not v.startswith('UN')]
coverage=[v for v in C.columns if v not in meta and '|' not in v and not v.startswith('UN')]
if set(base)!=set(coverage):raise ValueError('Unstratified pathway IDs mismatch')
if not A.age.notna().all():raise ValueError('Missing ages')
n=len(A)
info={'n_samples':n,'n_unique_ids':A.sample_id.nunique(),'n_pathways':len(base),'n_age_available':int(A.age.notna().sum()),'n_weight_available':int(A.weight.notna().sum()),'sample_ids_identical':True,'n_stratified_abundance':int(sum('|' in v for v in A.columns)),
'warning': ('Observed '+str(n)+' dogs, not 625 claimed in the Oct 15 email; do not call this the verified >=2-million-read cohort.' if n!=625 else None),
'no_read_counts_in_inputs':True,'periodontal_disease_labels_absent':True}
Path(a.out).write_text(json.dumps(info,indent=2))
print('PREFLIGHT',json.dumps(info))
