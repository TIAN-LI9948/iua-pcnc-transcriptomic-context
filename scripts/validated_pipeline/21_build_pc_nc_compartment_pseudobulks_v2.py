#!/usr/bin/env python
"""Build v2 donor pseudobulks from the verified H5AD without changing it."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

_context_root = os.environ.get("IUA_CONTEXT_ROOT")
if not _context_root:
    raise RuntimeError("Set IUA_CONTEXT_ROOT to the audited context-validation workspace before running this provenance script.")
ROOT = Path(_context_root).expanduser().resolve()
H5AD = ROOT / "raw" / "PC_NC_JCI2026" / "natproendometrium.h5ad"
BRANCHES = {
    "nonciliated_epithelial": {"glandular epithelial cell of endometrium", "luminal epithelial cell of endometrium"},
    "glandular_epithelial": {"glandular epithelial cell of endometrium"},
    "broad_epithelial_audit": {"glandular epithelial cell of endometrium", "luminal epithelial cell of endometrium", "ciliated epithelial cell"},
}

def decode(x): return x.decode() if isinstance(x, bytes) else str(x)
def stamp(): return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(4*1024*1024), b""): h.update(c)
    return h.hexdigest()
def category(h5, key):
    node=h5[f"obs/{key}"]; cats=[decode(x) for x in node["categories"][...]]
    return np.asarray([cats[int(i)] if int(i)>=0 else "NA" for i in node["codes"][...]], dtype=object)

def main():
    out_i, out_m = ROOT/"interim", ROOT/"metadata"
    expected=[]
    for branch in BRANCHES:
        expected.extend([out_i/f"PC_NC_{branch}_pseudobulk_counts_v2.tsv.gz",out_m/f"PC_NC_{branch}_pseudobulk_metadata_v2.tsv"])
    if any(x.exists() for x in expected):
        raise RuntimeError("Refusing to overwrite a v2 pseudobulk checkpoint.")
    with h5py.File(H5AD,"r") as h5:
        donor,cycle,ctype=category(h5,"donor_id"),category(h5,"CycleType"),category(h5,"cell_type")
        genes=np.asarray([decode(x) for x in h5["raw/var/gene_symbols"][...]],dtype=object)
        donors=sorted(pd.unique(donor).tolist())
        donor_index={x:i for i,x in enumerate(donors)}
        donor_cycle=pd.DataFrame({"donor_id":donor,"CycleType":cycle}).drop_duplicates().set_index("donor_id").loc[donors,"CycleType"].to_dict()
        raw=h5["raw/X"]; indptr, indices, data=raw["indptr"],raw["indices"],raw["data"]
        aggregate={b:np.zeros((len(genes),len(donors)),dtype=np.int64) for b in BRANCHES}
        nuclei={b:np.zeros(len(donors),dtype=np.int64) for b in BRANCHES}
        # Each row is allocated only to its predeclared branch(es); no outcome-based selection.
        for row in range(len(donor)):
            d=donor_index[donor[row]]; label=ctype[row]
            selected=[b for b, labels in BRANCHES.items() if label in labels]
            if not selected: continue
            start,end=int(indptr[row]),int(indptr[row+1])
            if start==end:
                for b in selected: nuclei[b][d]+=1
                continue
            idx=indices[start:end]; vals=data[start:end]
            if np.any(vals<0) or np.any(np.abs(vals-np.rint(vals))>1e-8): raise ValueError("Count semantics changed since phase B.")
            ivals=vals.astype(np.int64)
            for b in selected:
                np.add.at(aggregate[b][:,d],idx,ivals); nuclei[b][d]+=1
    qa=[]
    for b in BRANCHES:
        meta=pd.DataFrame({"pseudobulk_id":donors,"donor_id":donors,"CycleType":[donor_cycle[d] for d in donors],"nuclei":nuclei[b],"compartment":b})
        if (meta["nuclei"]<=0).any(): raise ValueError(f"{b} has donors with zero retained nuclei.")
        if meta.groupby("CycleType")["donor_id"].nunique().to_dict() != {"Controlled":7,"Natural":9}: raise ValueError(f"{b} lost a condition group.")
        table=pd.DataFrame(aggregate[b],columns=donors); table.insert(0,"gene",genes)
        table.to_csv(out_i/f"PC_NC_{b}_pseudobulk_counts_v2.tsv.gz",sep="\t",index=False,compression="gzip",quoting=csv.QUOTE_MINIMAL)
        meta.to_csv(out_m/f"PC_NC_{b}_pseudobulk_metadata_v2.tsv",sep="\t",index=False,quoting=csv.QUOTE_MINIMAL)
        qa.append({"branch":b,"labels":";".join(sorted(BRANCHES[b])),"donors":len(meta),"Controlled":int((meta.CycleType=="Controlled").sum()),"Natural":int((meta.CycleType=="Natural").sum()),"total_nuclei":int(meta.nuclei.sum()),"minimum_nuclei":int(meta.nuclei.min()),"all_counts_nonnegative_integer":bool(np.all(aggregate[b]>=0)),"h5ad_sha256":sha256(H5AD)})
    pd.DataFrame(qa).to_csv(out_m/"PC_NC_compartment_pseudobulk_qa_v2.tsv",sep="\t",index=False)
    print(json.dumps({"status":"complete","branches":qa,"completed_at":stamp()}))
if __name__=="__main__": main()
