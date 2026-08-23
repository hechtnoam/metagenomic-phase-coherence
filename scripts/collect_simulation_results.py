#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
BASES=("A","T","G","C")
CONDS={"random":(.5,.5),"b055_050":(.55,.5),"b070_050":(.7,.5),"b080_050":(.8,.5),"b090_050":(.9,.5),"b100_050":(1.,.5),"b100_100":(1.,1.)}
ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--read-length',type=int,default=59); a=ap.parse_args(); rows=[]; missing=[]
for genome in ('human','pseudomonas'):
 for cond,(b5,b3) in CONDS.items():
  for rep in range(1,11):
   d=a.root/genome/cond/f'replicate_{rep}'; meta=d/'simulated_reads.fasta.metadata.json'; p3=d/'period3'/'simulated_reads'/'stats'/'period3'/'kmer'/f'L{a.read_length}'/'k1_period3.json'
   if not meta.is_file() or not p3.is_file(): missing.append(str(d)); continue
   m=json.loads(meta.read_text()); j=json.loads(p3.read_text())
   if abs(float(m['bias_5prime'])-b5)>1e-12 or abs(float(m['bias_3prime'])-b3)>1e-12 or int(m['seed'])!=rep: raise ValueError(f'Parameter mismatch: {d}')
   fits={x['base'].upper():x for x in j['fits']}
   for base in BASES:
    x=fits[base]; rows.append(dict(genome=genome,condition=cond,bias_5prime=b5,bias_3prime=b3,replicate=rep,seed=rep,base=base,r2=float(x['r2']),amplitude=float(x['amplitude']),p_value=float(x['p_value']),p_value_adjusted=float(x['p_value_adjusted']),read_count=int(j['read_count']),simulation_attempts=int(m['attempts']),simulation_acceptance_fraction=float(m['acceptance_fraction'])))
if missing: raise RuntimeError('Missing/incomplete runs:\n'+'\n'.join(missing))
if len(rows)!=560: raise RuntimeError(f'Expected 560 rows, got {len(rows)}')
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f'Wrote {len(rows)} rows to {a.output}')