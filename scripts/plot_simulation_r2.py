#!/usr/bin/env python3
import argparse,json,hashlib
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
from scipy.stats import t
VERSION='2.1.0-audited-random-label'; BASES=('A','T','G','C')
ORDER=[('random',.5,.5,'random\nfragmentation'),('b055_050',.55,.5,'0.55/0.50'),('b070_050',.7,.5,'0.70/0.50'),('b080_050',.8,.5,'0.80/0.50'),('b090_050',.9,.5,'0.90/0.50'),('b100_050',1.,.5,'1.00/0.50'),('b100_100',1.,1.,'1.00/1.00')]
ap=argparse.ArgumentParser(); ap.add_argument('--input-csv',type=Path,required=True); ap.add_argument('--outdir',type=Path,required=True); a=ap.parse_args(); df=pd.read_csv(a.input_csv)
req={'genome','condition','bias_5prime','bias_3prime','replicate','base','r2'}
if req-set(df): raise ValueError(f'Missing columns: {sorted(req-set(df))}')
rep=df.groupby(['genome','condition','bias_5prime','bias_3prime','replicate'],as_index=False).agg(mean_r2_across_bases=('r2','mean'))
rows=[]
for keys,g in rep.groupby(['genome','condition','bias_5prime','bias_3prime'],sort=False):
 v=g.mean_r2_across_bases.to_numpy(float); n=len(v); mean=float(v.mean()); sd=float(v.std(ddof=1)) if n>1 else 0.; hw=float(t.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 else 0.
 rows.append(dict(genome=keys[0],condition=keys[1],bias_5prime=keys[2],bias_3prime=keys[3],n_replicates=n,mean_r2=mean,sd_between_replicates=sd,ci95_lower=max(0,mean-hw),ci95_upper=min(1,mean+hw)))
s=pd.DataFrame(rows); a.outdir.mkdir(parents=True,exist_ok=True); rep.to_csv(a.outdir/'simulation_replicate_means.csv',index=False); s.to_csv(a.outdir/'simulation_condition_summary.csv',index=False)
for genome in ('pseudomonas','human'):
 part=s[s.genome==genome]; rr=[]
 for cond,b5,b3,label in ORDER:
  hit=part[(part.condition==cond)&np.isclose(part.bias_5prime,b5)&np.isclose(part.bias_3prime,b3)]
  if len(hit)!=1: raise ValueError(f'{genome}: expected one row for {cond}, got {len(hit)}')
  rr.append((hit.iloc[0],label))
 vals=np.array([r.mean_r2 for r,_ in rr]); lo=np.array([r.ci95_lower for r,_ in rr]); hi=np.array([r.ci95_upper for r,_ in rr]); x=np.arange(7)
 fig,ax=plt.subplots(figsize=(9,5.5)); ax.errorbar(x,vals,yerr=np.vstack([vals-lo,hi-vals]),marker='o',linewidth=2,capsize=4); ax.set_xticks(x,[lab for _,lab in rr]); ax.set_xlabel("5' bias / 3' bias"); ax.set_ylabel('Mean period-3 R² across bases'); ax.set_ylim(0,1); ax.set_title('Human genome' if genome=='human' else 'Pseudomonas aeruginosa'); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(a.outdir/f'simulation_r2_{genome}.png',dpi=400,bbox_inches='tight'); fig.savefig(a.outdir/f'simulation_r2_{genome}.pdf',bbox_inches='tight'); plt.close(fig)
manifest={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'summary_rule':'Mean R2 across A/T/G/C within each replicate; mean and 95% t CI across independent replicate means.','random_fragmentation_definition':'bias_5prime=0.50 and bias_3prime=0.50; all A/C/G/T boundaries accepted equally'}
(a.outdir/'simulation_plot_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); print(f'Wrote plots to {a.outdir}')
