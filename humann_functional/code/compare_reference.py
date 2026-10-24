#!/usr/bin/env python3
"""Compare reproducibility headline values without hiding mismatches."""
from pathlib import Path
import argparse,json,math
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--actual',required=True);p.add_argument('--output',required=True);a=p.parse_args()
ref=json.loads(Path(a.reference).read_text());act=json.loads(Path(a.actual).read_text())
keys=['n_tested','cpm_fdr_lt05','clr_fdr_lt05','age_ge_1_cpm_fdr_lt05','weight_adjusted_cpm_fdr_lt05','supported_all_four','pathway_detection_vs_age_spearman_rho','pca_pc1_age_spearman_rho']
comp={k:{'expected':ref.get(k),'actual':act.get(k),'matches':(math.isclose(ref[k],act[k],rel_tol=1e-6,abs_tol=1e-6) if isinstance(ref.get(k),(int,float)) and isinstance(act.get(k),(int,float)) else ref.get(k)==act.get(k))} for k in keys}
Path(a.output).write_text(json.dumps(comp,indent=2))
if not all(v['matches'] for v in comp.values()):print('WARNING: some reference values differ; see',a.output)
else: print('Reference headline statistics reproduced')
