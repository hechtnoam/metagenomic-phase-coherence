#!/usr/bin/env python3
"""
Synthetic known-truth and iid-null validation for six_state_em.py.

Generates synthetic FASTA datasets with known latent six-state assignments, runs
the repository EM as a subprocess, and evaluates parameter/state recovery while
explicitly resolving the model's cyclic-row and strand-complement symmetries.

Validation conditions:
  positive : sequences generated from a known triplet-specific 3x4 matrix.
  iid_null : sequences generated from the positive model's overall base
             composition, with no position-class dependence.

The script creates independent replicate datasets, runs 20-restart hard-EM for
each, and writes machine-readable summaries.
"""
from __future__ import annotations
import argparse, csv, itertools, json, math, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

BASES = ("A","C","G","T")
B2I = {b:i for i,b in enumerate(BASES)}
STATES = ("+0","+1","+2","-0","-1","-2")
COMP = str.maketrans("ACGT","TGCA")
DEFAULT_P = np.array([
    [0.13,0.30,0.39,0.18],
    [0.28,0.17,0.24,0.31],
    [0.14,0.33,0.24,0.29],
], dtype=float)

def rc(s): return s.translate(COMP)[::-1]

def generate_from_state(p, length, state, rng):
    sign, off = state[0], int(state[1])
    # EM +o scores observed row r against model row (r+o) mod 3 because
    # np.roll(logp,-o)[r] = logp[(r+o)%3].
    chars=[]
    for i in range(length):
        model_row=(i+off)%3
        chars.append(rng.choice(BASES,p=p[model_row]))
    s="".join(chars)
    # A '-' state is represented by a reverse-complemented molecule; the EM
    # scores its reverse complement against the same underlying P.
    return rc(s) if sign=="-" else s

def generate_iid(freq, length, rng):
    return "".join(rng.choice(BASES,size=length,p=freq))

def write_dataset(path, truth_path, condition, p, n, length, seed):
    rng=np.random.default_rng(seed)
    truth=[]
    path.parent.mkdir(parents=True,exist_ok=True)
    overall=p.mean(axis=0); overall/=overall.sum()
    with path.open("w") as f:
        for i in range(n):
            state=STATES[i%6] if condition=="positive" else STATES[i%6]
            if condition=="positive":
                seq=generate_from_state(p,length,state,rng)
            else:
                seq=generate_iid(overall,length,rng)
            f.write(f">read_{i}|true_state={state}\n{seq}\n")
            truth.append((f"read_{i}",state))
    with truth_path.open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["read_id","true_state"]); w.writerows(truth)

def read_probs(path):
    df=pd.read_csv(path,sep="\t")
    L=int(df["length"].iloc[0])
    p=np.zeros((3,4))
    for _,r in df.iterrows():
        p[int(r.position_class)-1,B2I[str(r.base)]]=float(r.prob)
    return L,p

def counts_mod3(seq):
    c=np.zeros((3,4),dtype=np.int32)
    for i,b in enumerate(seq): c[i%3,B2I[b]]+=1
    return c

def score_states(seq,p):
    cp=counts_mod3(seq); cm=counts_mod3(rc(seq)); logp=np.log(p)
    scores=[]
    for c in (cp,cm):
        for o in range(3):
            scores.append(float(np.sum(c*np.roll(logp,-o,axis=0))))
    return np.array(scores)

def read_fasta(path):
    h=None; chunks=[]
    with path.open() as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None: yield h,"".join(chunks)
                h=line[1:].split("|")[0]; chunks=[]
            else: chunks.append(line)
    if h is not None: yield h,"".join(chunks)

def all_equivalent_truth_matrices(p,length):
    # Match six_state_em canonicalization transformations.
    col_comp=np.array([3,2,1,0])
    candidates=[]
    for rot in range(3):
        candidates.append((np.roll(p,-rot,axis=0),("rotation",rot,False)))
    mod=length%3
    row_src=[((mod-1)-q)%3 for q in range(3)]
    prc=p[row_src][:,col_comp]
    for rot in range(3):
        candidates.append((np.roll(prc,-rot,axis=0),("full6",rot,True)))
    return candidates

def matrix_recovery(est,true,length):
    best=None
    for equiv,desc in all_equivalent_truth_matrices(true,length):
        d=est-equiv
        rmse=float(np.sqrt(np.mean(d*d))); mx=float(np.max(np.abs(d)))
        if best is None or rmse<best[0]: best=(rmse,mx,desc,equiv)
    return best

def infer_states_from_est(fasta,p):
    pred={}
    for rid,seq in read_fasta(fasta):
        pred[rid]=STATES[int(np.argmax(score_states(seq,p)))]
    return pred

def normalize_state_label(x):
    """Normalize truth/prediction state labels to one of +0,+1,+2,-0,-1,-2."""
    s = str(x).strip()
    # Preserve explicit signed labels.
    if s in STATES:
        return s
    # Pandas may parse +0,+1,+2 as numeric 0,1,2.
    if s in {"0", "0.0", "1", "1.0", "2", "2.0"}:
        return f"+{int(float(s))}"
    raise ValueError(f"Invalid six-state label: {x!r}")

def best_state_accuracy(truth,pred):
    # State labels have a 6-element symmetry group. For validation, maximize
    # accuracy over all bijections induced by the 3 rotations and sign flip.
    true=[normalize_state_label(truth[k]) for k in truth]
    pr=[normalize_state_label(pred[k]) for k in truth]
    maps=[]
    for flip in (False,True):
        for rot in range(3):
            m={}
            for s in STATES:
                sign=s[0]; o=int(s[1])
                if flip: sign="-" if sign=="+" else "+"
                m[s]=f"{sign}{(o-rot)%3}"
            maps.append(m)
    best=(0,None)
    for m in maps:
        acc=np.mean([m[t]==q for t,q in zip(true,pr)])
        if acc>best[0]: best=(float(acc),m)
    return best

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--em-script",type=Path,default=Path("scripts/six_state_em.py"))
    ap.add_argument("--outdir",type=Path,default=Path("results/em_synthetic_validation"))
    ap.add_argument("--replicates",type=int,default=10)
    ap.add_argument("--reads",type=int,default=60000)
    ap.add_argument("--length",type=int,default=60)
    ap.add_argument("--restarts",type=int,default=20)
    ap.add_argument("--threads",type=int,default=10)
    ap.add_argument("--seed",type=int,default=1000)
    ap.add_argument("--max-iter",type=int,default=100)
    ap.add_argument("--tol",type=float,default=1e-8)
    ap.add_argument("--alpha",type=float,default=0.5)
    a=ap.parse_args()
    if not a.em_script.is_file(): raise FileNotFoundError(a.em_script)
    a.outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for cond in ("positive","iid_null"):
        for rep in range(1,a.replicates+1):
            seed=a.seed + (0 if cond=="positive" else 100000) + rep
            rd=a.outdir/cond/f"replicate_{rep}"
            fasta=rd/"synthetic.fasta"; truthp=rd/"truth.tsv"; emout=rd/"em"
            rd.mkdir(parents=True,exist_ok=True)
            write_dataset(fasta,truthp,cond,DEFAULT_P,a.reads,a.length,seed)
            cmd=[sys.executable,str(a.em_script),"--input",str(fasta),"--input-type","fasta",
                 "--out",str(emout),"--lengths",str(a.length),"--dedup","none",
                 "--alpha",str(a.alpha),"--max-iter",str(a.max_iter),"--tol",str(a.tol),
                 "--seed",str(seed),"--restarts",str(a.restarts),"--threads",str(a.threads),
                 "--canonicalize","full6"]
            subprocess.run(cmd,check=True)
            _,est=read_probs(emout/"base_probs_by_length.tsv")
            rmse,mx,desc,_=matrix_recovery(est,DEFAULT_P,a.length)
            truth_df = pd.read_csv(
                truthp,
                sep="\t",
                dtype={"read_id": str, "true_state": str},
                keep_default_na=False,
            )
            truth = dict(zip(truth_df["read_id"], truth_df["true_state"]))
            pred=infer_states_from_est(fasta,est)
            state_acc,_=best_state_accuracy(truth,pred)
            restart = pd.read_csv(emout / "restart_diagnostics_by_length.tsv", sep="\t")
            # six_state_em.py v4.3 writes the final complete-data hard-assignment
            # log score as final_LL (and delta_LL_from_best relative to the winner).
            if "final_LL" not in restart.columns:
                raise ValueError(
                    "restart_diagnostics_by_length.tsv is missing required column "
                    f"'final_LL'; found: {list(restart.columns)}"
                )
            best = float(restart["final_LL"].max())
            near = int(np.sum(np.abs(restart["final_LL"] - best) <= 1e-6))
            # Null structure: deviation among position rows.
            row_spread=float(np.sqrt(np.mean((est-est.mean(axis=0,keepdims=True))**2)))
            rows.append(dict(condition=cond,replicate=rep,seed=seed,n_reads=a.reads,
                             length=a.length,matrix_rmse_to_positive_truth=rmse,
                             matrix_max_abs_error_to_positive_truth=mx,
                             state_accuracy_after_symmetry=state_acc,
                             estimated_position_row_spread=row_spread,
                             restarts=a.restarts,restarts_at_best_objective=near,
                             best_objective=best))
            print(f"[{cond} rep {rep}] RMSE={rmse:.5g} state_acc={state_acc:.4f} row_spread={row_spread:.5g}")
    df=pd.DataFrame(rows)
    df.to_csv(a.outdir/"validation_summary.tsv",sep="\t",index=False)
    summary=df.groupby("condition").agg(
        n_replicates=("replicate","count"),
        mean_matrix_rmse=("matrix_rmse_to_positive_truth","mean"),
        mean_state_accuracy=("state_accuracy_after_symmetry","mean"),
        mean_position_row_spread=("estimated_position_row_spread","mean"),
        mean_restarts_at_best=("restarts_at_best_objective","mean"),
    ).reset_index()
    summary.to_csv(a.outdir/"validation_condition_summary.tsv",sep="\t",index=False)
    manifest={"script":"validate_six_state_em.py","conditions":["positive","iid_null"],
              "true_matrix":{"bases":list(BASES),"rows":DEFAULT_P.tolist()},
              "state_proportions":"exactly balanced across six states in positive datasets",
              "null":"iid bases drawn from mean composition of positive matrix; no position-class dependence",
              "args":vars(a)}
    # stringify Paths
    manifest["args"]={k:str(v) if isinstance(v,Path) else v for k,v in manifest["args"].items()}
    (a.outdir/"validation_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("\nWrote:",a.outdir/"validation_summary.tsv")
    print("Wrote:",a.outdir/"validation_condition_summary.tsv")

if __name__=="__main__":
    main()
