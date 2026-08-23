#!/usr/bin/env python3
"""Create same-ID, same-length frame-normalized and matched-control FASTAs."""
import argparse, hashlib, json
from pathlib import Path

VERSION = "1.0.0"

def canonical_id(header: str) -> str:
    """Ordinal|read_name key shared by normalized and matched-control headers."""
    parts = header.split("|")
    if len(parts) < 2:
        return header
    return f"{parts[0]}|{parts[1]}"

def records(path):
    h = None; seq = []
    with Path(path).open() as f:
        for raw in f:
            s = raw.strip()
            if not s: continue
            if s.startswith(">"):
                if h is not None: yield h, "".join(seq).upper()
                h = canonical_id(s[1:].split()[0]); seq = []
            else:
                if h is None: raise ValueError(f"Sequence before header in {path}")
                seq.append(s)
    if h is not None: yield h, "".join(seq).upper()

def sha(path):
    x = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""): x.update(b)
    return x.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--normalized", required=True, type=Path)
    ap.add_argument("--unnormalized", required=True, type=Path)
    ap.add_argument("--length", required=True, type=int)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--prefix")
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    prefix = a.prefix or f"L{a.length}"
    no = a.outdir / f"{prefix}.frame_normalized.fasta"
    co = a.outdir / f"{prefix}.matched_unnormalized.fasta"
    mo = a.outdir / f"{prefix}.matched_comparison.metadata.json"

    ctrl = {}
    for rid, seq in records(a.unnormalized):
        if rid in ctrl: raise ValueError(f"Duplicate control ID: {rid}")
        ctrl[rid] = seq

    seen=set(); n=missing=short=dup=bad=0
    with no.open("w") as nf, co.open("w") as cf:
        for rid, norm in records(a.normalized):
            if len(norm) != a.length: continue
            if rid in seen: dup += 1; continue
            seen.add(rid)
            raw = ctrl.get(rid)
            if raw is None: missing += 1; continue
            if len(raw) < a.length: short += 1; continue
            raw = raw[:a.length]
            if (set(norm) | set(raw)) - set("ACGT"): bad += 1; continue
            nf.write(f">{rid}\n{norm}\n")
            cf.write(f">{rid}\n{raw}\n")
            n += 1

    ni=[r for r,_ in records(no)]; ci=[r for r,_ in records(co)]
    if not n or ni != ci: raise RuntimeError("Matched-output validation failed")
    meta = {
      "script": Path(__file__).name, "script_version": VERSION,
      "normalized_input": str(a.normalized),
      "matched_unnormalized_input": str(a.unnormalized),
      "selected_normalized_length": a.length,
      "control_definition": "Same accepted read IDs; first L bases in coding orientation before phase-dependent trimming.",
      "counts": {"matched_pairs_written": n, "missing_control": missing,
                 "control_shorter_than_requested_length": short,
                 "duplicate_normalized_id_excluded": dup, "non_acgt_excluded": bad},
      "validation": {"identical_id_order": True, "normalized_records": len(ni),
                     "control_records": len(ci), "all_output_sequences_length": a.length},
      "outputs": {"frame_normalized": str(no), "matched_unnormalized": str(co)},
      "sha256": {"normalized_input": sha(a.normalized), "matched_unnormalized_input": sha(a.unnormalized),
                 "frame_normalized_output": sha(no), "matched_unnormalized_output": sha(co)}
    }
    mo.write_text(json.dumps(meta, indent=2)+"\n")
    print(f"Matched pairs written: {n:,}")
    print(no); print(co); print(mo)

if __name__ == "__main__": main()