import pandas as pd
import numpy as np
from scipy import stats as sstats

# ------------------------------------------------------------------
# COMPLEX-LEVEL (module) co-dependency: instead of individual genes,
# compute a mean-dependency "module score" for each SWI/SNF sub-complex
# per cell line, then correlate ARID1A dependency with each module.
# Modules defined A PRIORI by canonical complex biology (Mashtalir 2018),
# NOT by the earlier result, so this is not circular. ARID1A itself is
# excluded from every module. Predicts (if the collateral-dependency is
# core-complex-specific): core+cBAF strongly positive, PBAF ~0, ncBAF ~0/neg.
# ------------------------------------------------------------------

MODULES = {
    "Core_shared":   ["SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCE1","ACTL6A","ACTB"],
    "cBAF_specific": ["ARID1B","DPF1","DPF2","DPF3","BCL7A","BCL7B","BCL7C","SS18","SS18L1"],
    "PBAF_specific": ["PBRM1","ARID2","BRD7","PHF10"],
    "ncBAF_specific":["BRD9","BICRA","BICRAL"],
    "Alt_ATPase":    ["SMARCA2"],   # BRM, alternative to SMARCA4 -- known SMARCA4-synthetic-lethal, expect NOT co-dependent
}

models = pd.read_csv("model_list.csv", low_memory=False)
lung_ids = set(models.loc[models["tissue"].str.contains("Lung", case=False, na=False), "model_id"])

fname = "Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_fold_change_values_20250624.tsv"
with open(fname, encoding="utf-8") as f:
    f.readline(); row2 = f.readline().rstrip("\n").split("\t"); row3 = f.readline().rstrip("\n").split("\t")
model_ids_row = row2[3:]; sources_row = row3[3:]
df = pd.read_csv(fname, sep="\t", skiprows=5, header=None, low_memory=False)
ncols = df.shape[1] - 3
genes = df.iloc[:, 1].astype(str).values
mat = df.iloc[:, 3:3+ncols].apply(pd.to_numeric, errors="coerce")
model_ids_row = model_ids_row[:ncols]; sources_row = sources_row[:ncols]

src_rank = np.array([0 if s=="Broad" else 1 for s in sources_row])
colinfo = pd.DataFrame({"idx": range(ncols), "model_id": model_ids_row, "src_rank": src_rank})
keep = np.sort(colinfo.sort_values(["model_id","src_rank"]).drop_duplicates("model_id",keep="first")["idx"].values)
mat = mat.iloc[:, keep]
kept_ids = np.array(model_ids_row)[keep]
gdf = pd.DataFrame(mat.values, index=genes).groupby(level=0).mean()
gdf.columns = kept_ids

def analyze(col_ids, label):
    sel = [c for c in gdf.columns if c in col_ids]
    sub = gdf[sel]
    arid = sub.loc["ARID1A"] if "ARID1A" in sub.index else None
    print(f"\n===== [{label}] complex-module co-dependency with ARID1A ({len(sel)} lines) =====")
    rows = []
    for mod, mgenes in MODULES.items():
        present = [g for g in mgenes if g in sub.index]
        if not present:
            continue
        modscore = sub.loc[present].mean(axis=0, skipna=True)  # per-line mean dependency of the module
        pair = pd.concat([arid, modscore], axis=1).dropna()
        pair.columns = ["ARID1A","module"]
        if len(pair) < 10:
            continue
        rho, p = sstats.spearmanr(pair["ARID1A"], pair["module"])
        rows.append({"module": mod, "n_subunits": len(present), "subunits": ";".join(present),
                    "n_lines": len(pair), "rho": rho, "P": p})
        print(f"  {mod:<15} ({len(present)} subunits) rho={rho:+.3f} P={p:.2e}")
    out = pd.DataFrame(rows)
    out["tissue"] = label
    return out

res = pd.concat([analyze(lung_ids, "lung"), analyze(set(gdf.columns), "pancancer")], ignore_index=True)
res.to_csv("complex_module_codependency.csv", index=False)
print("\nSaved complex_module_codependency.csv")
