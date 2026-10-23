"""Descriptive post-freeze audit: balance in final paired complete-case train rows.

This does not change the matching, include test outcomes, or rerun endpoints.
The same original training-background scales are used for both null types.
"""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import json
from pathlib import Path
import numpy as np
import pandas as pd

root=Path("/root/autodl-tmp/arid1a_baf_r15/strong_null")
primary=pd.read_csv(root/"r15_primary_complex_results.csv")
keep=set(primary.source_complex)
results=pd.read_csv(root/"r15_pair_results.csv")
results=results.loc[results.source_complex.isin(keep)&results.status.eq("ok")]
catalog=pd.read_csv(root/"r15_frozen_complexes.csv").set_index("source_complex")
nulls=pd.read_csv(root/"r15_null_gene_manifest.csv").set_index(["source_complex","split","kind","null_rep"])
matches=pd.read_csv(root/"r15_gene_matching_balance.csv")
maps={(cid,int(s),int(rep)):dict(zip(d.real_gene,d.matched_gene)) for (cid,s,rep),d in matches.loc[matches.source_complex.isin(keep)].groupby(["source_complex","split","null_rep"])}
splits=pd.read_csv(root/"r15_patient_splits.csv")
trains={int(s):d.loc[d.partition.eq("train"),"row"].to_numpy() for s,d in splits.groupby("split")}
arr=np.load(root/"dependency.npy",mmap_mode="r")
original=np.load(root/"training_features.npy",mmap_mode="r")
genes=json.loads((root/"genes.json").read_text())
gi={g:i for i,g in enumerate(genes)}
scales={}
for split in range(10):
    feat=original[split]
    eligible=(feat[:,4]>=150)&(feat[:,1]>1e-8)&np.isfinite(feat[:,:4]).all(axis=1)
    scale=np.std(feat[eligible,:4],axis=0)
    scale[scale<1e-12]=1
    scales[split]=scale

def features(a):
    return np.column_stack([a.mean(axis=0),a.std(axis=0),(a<-1).mean(axis=0),np.zeros(a.shape[1])])

rows=[]
feature_names=["mean","sd_ddof0","fraction_observed_below_minus1","missing_fraction"]
for r in results.itertuples():
    rg=catalog.loc[r.source_complex,"genes"].split(";")
    ng=nulls.loc[(r.source_complex,r.split,r.kind,r.null_rep),"genes"].split(";")
    ri=[gi[g] for g in rg]
    ni=[gi[g] for g in ng]
    train=trains[r.split]
    block=arr[np.ix_(train,ri+ni)]
    block=block[np.isfinite(block).all(axis=1)]
    assert len(block)==r.train_complete
    rf,nf=features(block[:,:len(ri)]),features(block[:,len(ri):])
    delta=(nf.mean(axis=0)-rf.mean(axis=0))/scales[r.split]
    row={"source_complex":r.source_complex,"split":r.split,"kind":r.kind,"null_rep":r.null_rep,"n_train_complete":len(block),"k_total":len(ri)}
    for j,f in enumerate(feature_names):
        row[f"module_mean_standardized_delta_{f}"]=float(delta[j])
    if r.kind=="strong":
        mapping=maps[(r.source_complex,r.split,r.null_rep)]
        target_order=[ng.index(mapping[g]) for g in rg]
        gene_delta=(nf[target_order]-rf)/scales[r.split]
        for j,f in enumerate(feature_names):
            row[f"gene_pair_mean_abs_delta_{f}"]=float(abs(gene_delta[:,j]).mean())
            row[f"gene_pair_max_abs_delta_{f}"]=float(abs(gene_delta[:,j]).max())
            row[f"gene_pair_above_original_caliper_{f}"]=int((abs(gene_delta[:,j])>.5+1e-10).sum())
    rows.append(row)

df=pd.DataFrame(rows)
df.to_csv(root/"r15_final_training_balance_pairs.csv",index=False)
summary=[]
for kind,d in df.groupby("kind"):
    for f in feature_names:
        x=d[f"module_mean_standardized_delta_{f}"].to_numpy()
        row={"kind":kind,"feature":f,"n_primary_pairs":len(d),"mean_signed_module_delta":float(x.mean()),"mean_abs_module_delta":float(abs(x).mean()),"p95_abs_module_delta":float(np.quantile(abs(x),.95)),"max_abs_module_delta":float(abs(x).max())}
        if kind=="strong":
            row["mean_within_module_gene_pair_abs_delta"]=float(d[f"gene_pair_mean_abs_delta_{f}"].mean())
            row["max_gene_pair_abs_delta"]=float(d[f"gene_pair_max_abs_delta_{f}"].max())
            row["n_gene_pairs_above_original_caliper"]=int(d[f"gene_pair_above_original_caliper_{f}"].sum())
            row["n_gene_pairs_evaluated"]=int(d.k_total.sum())
        summary.append(row)
out=pd.DataFrame(summary)
out.to_csv(root/"r15_final_training_balance_summary.csv",index=False)
print(out.to_string(index=False),flush=True)
print("Complete-case missing fraction is mechanically zero in both arms; this is not independent evidence of matching quality.",flush=True)
