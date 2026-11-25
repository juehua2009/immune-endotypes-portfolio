import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import nnls
import scanpy as sc
import decoupler as dc

DATA = Path.home() / "Documents/immune-endotypes/data"
MOD = Path.home() / "Documents/github/immune-endotypes-portfolio/modules/05_cross_cohort"
RES = MOD / "results"

SIG8 = ["IL11", "IL13RA2", "CHI3L1", "MMP3", "MMP10", "IL6", "CXCL1", "CXCL5"]
T = 0.107229
import json
UNI = set(json.load(open(RES / "universe_human.json"))["genes"])
rng = np.random.default_rng(20261123)

def score_iaf_v3(a):
    genes = [g for g in a.var_names if g in UNI]
    sub = a[:, genes].copy()
    net = pd.DataFrame([{"source": "SIG8", "target": g, "weight": 1.0} for g in SIG8])
    dc.mt.aucell(data=sub, net=net, tmin=1, verbose=False)
    s = sub.obsm["score_aucell"]["SIG8"].values
    return s >= T

# ---------- load + label (v3 ruler replaces ruler_tool SIG9) ----------
sm = sc.read_h5ad(DATA / "smillie_uc/adata_qc.h5ad")
print("Smillie:", sm.shape)
FIB_CLUSTERS = ["WNT2B+ Fos-lo 1", "WNT2B+ Fos-lo 2", "WNT2B+ Fos-hi",
                "WNT5B+ 1", "WNT5B+ 2", "RSPO3+",
                "Inflammatory Fibroblasts", "Myofibroblasts"]
missing = [c for c in FIB_CLUSTERS if c not in set(sm.obs["Cluster"].astype(str))]
assert not missing, f"fibroblast clusters missing: {missing}"
sm_fib = sm[sm.obs["Cluster"].astype(str).isin(FIB_CLUSTERS)].copy()
sm_fib.obs["ref_group"] = np.where(score_iaf_v3(sm_fib), "IAF", "Fibroblast_other")
print("Smillie fibroblast reference:", sm_fib.n_obs, "|", sm_fib.obs["ref_group"].value_counts().to_dict())

def compartment(cluster):
    c = str(cluster)
    if any(k in c for k in ["Plasma", "Follicular"]) or c.startswith("B "): return "B_plasma"
    if any(k in c for k in ["CD4", "CD8", "Treg", "NK", "ILC"]): return "T_NK"
    if any(k in c for k in ["Macrophage", "DC", "Mast", "Monocyte", "Neutrophil"]): return "Myeloid"
    if any(k in c for k in ["Endothelial", "Venule", "Microvascular", "Lymphatic"]): return "Endothelial"
    if any(k in c for k in ["Enterocyte", "TA", "Goblet", "Stem", "Best4", "Tuft",
                            "Enteroendocrine", "Secretory", "Cycling"]): return "Epithelial"
    return "Other"
sm_other = sm[~sm.obs["Cluster"].astype(str).isin(FIB_CLUSTERS)].copy()
sm_other.obs["ref_group"] = [compartment(c) for c in sm_other.obs["Cluster"]]
del sm

ma_fib = sc.read_h5ad(DATA / "martin_2019/adata_martin_fibroblast.h5ad")
print("Martin fibroblasts:", ma_fib.shape)
ma_fib.obs["ref_group"] = np.where(score_iaf_v3(ma_fib), "IAF", "Fibroblast_other")
print("Martin fibroblast reference:", ma_fib.obs["ref_group"].value_counts().to_dict())

# ---------- donor split (same seed/order as nb12) ----------
def donor_col(obs):
    for c in ["Subject", "subject", "donor", "Donor", "patient", "Patient", "sample_id"]:
        if c in obs.columns: return c
    return None

def lognorm(a):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a

sm_fib = lognorm(sm_fib); ma_fib = lognorm(ma_fib); sm_other = lognorm(sm_other)
for a in (sm_fib, ma_fib):
    dc_ = donor_col(a.obs)
    donors = a.obs[dc_].astype(str).unique()
    rng.shuffle(donors)
    a.obs["donor_split"] = np.where(a.obs[dc_].astype(str).isin(donors[:max(1, int(.7 * len(donors)))]),
                                    "reference", "heldout")
    print("donor split:", a.obs["donor_split"].value_counts().to_dict(), "| donor col:", dc_)

# ---------- reference machinery (verbatim from nb12) ----------
def harmonize(objs):
    common = objs[0].var_names
    for a in objs[1:]:
        common = common.intersection(a.var_names)
    print("common genes across reference objects:", len(common))
    return [a[:, common].copy() for a in objs]

def group_profiles(objs, group_col="ref_group"):
    genes = objs[0].var_names
    sums, counts = {}, {}
    for a in objs:
        g = a.obs[group_col].astype(str)
        for grp in g.unique():
            Xg = a.X[g.values == grp].tocsr().copy()
            Xg.data = np.expm1(Xg.data)
            sums[grp] = sums.get(grp, 0) + np.asarray(Xg.sum(axis=0)).ravel()
            counts[grp] = counts.get(grp, 0) + int((g.values == grp).sum())
    return pd.DataFrame({grp: sums[grp] / counts[grp] for grp in sorted(sums)}, index=genes)

def detection_fractions(objs, group_col="ref_group"):
    genes = objs[0].var_names
    hits, counts = {}, {}
    for a in objs:
        g = a.obs[group_col].astype(str)
        for grp in g.unique():
            Xg = a.X[g.values == grp]
            hits[grp] = hits.get(grp, 0) + np.asarray((Xg > 0).sum(axis=0)).ravel()
            counts[grp] = counts.get(grp, 0) + int((g.values == grp).sum())
    return pd.DataFrame({grp: hits[grp] / counts[grp] for grp in sorted(hits)}, index=genes)

def select_markers(prof, detfrac, topn=60, min_frac=0.10):
    markers = {}
    for grp in prof.columns:
        lfc = np.log2((prof[grp] + 0.01) / (prof.drop(columns=[grp]).max(axis=1) + 0.01))
        cand = lfc[detfrac[grp] >= min_frac]
        markers[grp] = list(cand.sort_values(ascending=False).head(topn).index)
    return markers

def dwls(bulk_vec, S, iters=200, damp=0.1):
    x, _ = nnls(S.values, bulk_vec)
    if x.sum() == 0: return pd.Series(0.0, index=S.columns)
    x = x / x.sum()
    for _ in range(iters):
        w = 1.0 / (S.values @ x + damp)
        Sw = S.values * np.sqrt(w)[:, None]; bw = bulk_vec * np.sqrt(w)
        xn, _ = nnls(Sw, bw)
        if xn.sum() == 0: break
        xn = xn / xn.sum()
        if np.max(np.abs(xn - x)) < 1e-5: x = xn; break
        x = xn
    return pd.Series(x, index=S.columns)

def nnls_plain(bulk_vec, S):
    x, _ = nnls(S.values, bulk_vec)
    return pd.Series(x / x.sum() if x.sum() > 0 else x, index=S.columns)

ref_objs = harmonize([sm_fib[sm_fib.obs["donor_split"] == "reference"],
                      ma_fib[ma_fib.obs["donor_split"] == "reference"], sm_other])
prof = group_profiles(ref_objs)
detfrac = detection_fractions(ref_objs)
markers = select_markers(prof, detfrac)
all_markers = sorted(set(sum(markers.values(), [])))
S_ref = prof.loc[all_markers]
print("reference groups:", list(prof.columns))
print("markers:", {k: len(v) for k, v in markers.items()}, "| union:", len(all_markers))
S_ref.to_csv(RES / "ruler_v3_reference_signatures.csv")

# ---------- the gate (verbatim design from nb12) ----------
def titrated_pseudobulk(fib_a, other_a, target_frac, n_fib=800, n_other=3000, seed=0):
    r = np.random.default_rng(seed)
    iaf_idx = np.where(fib_a.obs["ref_group"].values == "IAF")[0]
    oth_idx = np.where(fib_a.obs["ref_group"].values == "Fibroblast_other")[0]
    n_iaf = int(round(target_frac * n_fib)); n_oth = n_fib - n_iaf
    fi = np.concatenate([
        r.choice(iaf_idx, size=n_iaf, replace=len(iaf_idx) < n_iaf),
        r.choice(oth_idx, size=n_oth, replace=len(oth_idx) < n_oth)])
    oi = r.choice(other_a.n_obs, size=min(n_other, other_a.n_obs), replace=False)
    Xf = fib_a.layers["counts"] if "counts" in fib_a.layers else fib_a.X
    Xo = other_a.layers["counts"] if "counts" in other_a.layers else other_a.X
    sf = pd.Series(np.asarray(Xf[fi].sum(axis=0)).ravel(), index=fib_a.var_names)
    so = pd.Series(np.asarray(Xo[oi].sum(axis=0)).ravel(), index=other_a.var_names)
    bulk = sf.add(so, fill_value=0)
    bulk = bulk / bulk.sum() * 1e4
    return bulk, float(fib_a.obs["ref_group"].iloc[fi].eq("IAF").mean())

LEVELS = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50]
sm_ho = sm_fib[sm_fib.obs["donor_split"] == "heldout"]
rows = []
for f in LEVELS:
    for s in range(4):
        bv, true_f = titrated_pseudobulk(sm_ho, sm_other, f, seed=int(f * 1000) + s)
        for solver_name, solver in [("dwls", dwls), ("nnls", nnls_plain)]:
            est = solver(bv.reindex(S_ref.index).fillna(0).values, S_ref)
            denom = est["IAF"] + est["Fibroblast_other"]
            rows.append({"solver": solver_name, "target": f, "seed": s,
                         "true_iaf_frac": true_f,
                         "est_iaf_frac": float(est["IAF"] / denom) if denom > 0 else 0.0})
gate = pd.DataFrame(rows)
gate.to_csv(RES / "ruler_v3_deconv_gate.csv", index=False)
for solver_name, gsub in gate.groupby("solver"):
    rho = stats.spearmanr(gsub.true_iaf_frac, gsub.est_iaf_frac).statistic
    mae = float((gsub.est_iaf_frac - gsub.true_iaf_frac).abs().median())
    print(f"{solver_name}: rho = {rho:.3f} | median |err| = {mae * 100:.1f} pp | GATE {'PASS' if (rho >= 0.80 and mae <= 0.10) else 'FAIL'} (pass = rho>=0.80 and |err|<=10pp)")
print("STEP8 DONE")
