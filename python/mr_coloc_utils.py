"""Shared: 1000G EUR bed LD utils, MR methods, coloc.abf (Giambartolomei 2014)."""
import numpy as np
from scipy import stats

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
NREF = 503
ROWB = (NREF + 3) // 4

_lut = None
def _get_lut():
    global _lut
    if _lut is None:
        _lut = np.full((256, 4), np.nan)
        for b in range(256):
            for j in range(4):
                c = (b >> (2 * j)) & 3
                _lut[b, j] = {0: 0.0, 2: 1.0, 3: 2.0, 1: np.nan}[c]
    return _lut

def load_bim_index():
    import pandas as pd, os
    f = os.path.join(BASE, "data", "gwas", "bim_full_index.csv")
    import os.path as op
    if op.exists(f):
        return pd.read_csv(f)
    rows = []
    with open(os.path.join(BASE, "g1000_eur", "g1000_eur.bim"), "r", buffering=1 << 20) as fh:
        for i, line in enumerate(fh):
            p = line.split(None, 3)
            if len(p) >= 4:
                rows.append((p[1], i))
    df = pd.DataFrame(rows, columns=["rsid", "bim_row"])
    df.to_csv(f, index=False)
    return df

def genotypes_for(rsids):
    """Return dict rsid -> dosage vector (NREF,), NaN=missing."""
    import os
    bim = load_bim_index()
    sub = bim[bim["rsid"].isin(set(rsids))]
    bed = np.memmap(os.path.join(BASE, "g1000_eur", "g1000_eur.bed"), dtype=np.uint8, mode="r")
    lut = _get_lut()
    out = {}
    for rs, row in zip(sub["rsid"], sub["bim_row"]):
        o = 3 + int(row) * ROWB
        raw = np.frombuffer(bed[o:o + ROWB], dtype=np.uint8)
        out[rs] = lut[raw].reshape(-1)[:NREF]
    return out

def ld_r2(g1, g2):
    ok = ~(np.isnan(g1) | np.isnan(g2))
    if ok.sum() < 50: return 0.0
    a, b = g1[ok], g2[ok]
    if a.std() == 0 or b.std() == 0: return 0.0
    return float(np.corrcoef(a, b)[0, 1] ** 2)

def clump(df, pcol="p_exposure", chr_col="chr", pos_col="pos", r2_thr=0.1, win=250000):
    """df sorted by p; greedy clump using 1000G EUR genotypes."""
    df = df.sort_values(pcol).reset_index(drop=True)
    rs = df["SNP"].tolist()
    geno = genotypes_for(rs)
    keep_idx = []
    for i in range(len(df)):
        r = rs[i]
        if r not in geno: continue
        drop = False
        for j in keep_idx:
            rj = rs[j]
            if df.loc[j, chr_col] == df.loc[i, chr_col] and abs(df.loc[j, pos_col] - df.loc[i, pos_col]) <= win:
                if ld_r2(geno[r], geno[rj]) > r2_thr:
                    drop = True; break
        if not drop:
            keep_idx.append(i)
    return df.loc[keep_idx].reset_index(drop=True)

# ---------------- MR ----------------
def wald_ratio(bx, sx, by, sy):
    b = by / bx
    se = abs(sy / bx)  # first-order
    p = 2 * stats.norm.sf(abs(b / se))
    return b, se, p

def ivw(bx, sx, by, sy, multiplicative=True):
    w = 1.0 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se_fe = 1.0 / np.sqrt(np.sum(w * bx ** 2))
    q = np.sum(w * (by - b * bx) ** 2)
    dfree = len(bx) - 1
    q_p = stats.chi2.sf(q, dfree) if dfree > 0 else np.nan
    if multiplicative and dfree > 0:
        phi = max(q / dfree, 1.0)
        se = se_fe * np.sqrt(phi)
    else:
        se = se_fe
    p = 2 * stats.norm.sf(abs(b / se))
    return dict(b=b, se=se, p=p, Q=q, Q_df=dfree, Q_p=q_p)

def mr_egger(bx, sx, by, sy):
    w = 1.0 / sy ** 2
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ by)
    resid = by - X @ beta
    dfree = len(bx) - 2
    q = float(np.sum(w * resid ** 2))
    phi = max(q / dfree, 1e-8) if dfree > 0 else 1.0
    covb = phi * np.linalg.inv(XtWX)
    a, b = beta
    se_a, se_b = np.sqrt(covb[0, 0]), np.sqrt(covb[1, 1])
    p_a = 2 * stats.t.sf(abs(a / se_a), max(dfree, 1))
    p_b = 2 * stats.t.sf(abs(b / se_b), max(dfree, 1))
    return dict(intercept=a, intercept_se=se_a, intercept_p=p_a,
                b=b, se=se_b, p=p_b, Q=q, Q_df=dfree)

# ---------------- coloc.abf ----------------
def _lbf(beta, se, W):
    V = se ** 2
    r = W ** 2 / (W ** 2 + V)
    z = beta / se
    return 0.5 * np.log1p(-r) + 0.5 * z ** 2 * r

def _logsumexp(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))

def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, W1=0.15, W2=0.2):
    """Giambartolomei et al. 2014. W1: quantitative (0.15); W2: cc (0.2)."""
    l1 = _lbf(beta1, se1, W1)
    l2 = _lbf(beta2, se2, W2)
    n = len(l1)
    lp = np.full(5, -np.inf)
    lp[0] = 0.0                                    # H0
    lp[1] = np.log(p1) + _logsumexp(l1)            # H1
    lp[2] = np.log(p2) + _logsumexp(l2)            # H2
    s1, s2 = _logsumexp(l1), _logsumexp(l2)
    both = _logsumexp(l1 + l2)
    # H3: p1*p2 * (sum_i!=j e^{l1_i+l2_j}) = p1p2 (e^{s1+s2} - e^{both})
    a = s1 + s2
    diff = np.exp(a - both) - 1.0
    lp[3] = np.log(p1) + np.log(p2) + both + np.log(max(diff, 1e-300))
    lp[4] = np.log(p12) + both                     # H4
    pp = np.exp(lp - _logsumexp(lp))
    return dict(nsnps=n, PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3], PP_H4=pp[4])

def _coloc_selftest():
    rng = np.random.default_rng(1)
    n = 200
    b1 = rng.normal(0, 0.05, n); s1 = np.full(n, 0.05)
    b2 = rng.normal(0, 0.05, n); s2 = np.full(n, 0.05)
    r = coloc_abf(b1, s1, b2, s2)
    tot = r["PP_H0"] + r["PP_H1"] + r["PP_H2"] + r["PP_H3"] + r["PP_H4"]
    assert abs(tot - 1) < 1e-9, tot
    # null data -> H0 dominant
    assert r["PP_H0"] > 0.9, r
    # strong shared signal at SNP k -> H4 dominant
    b1[50] = 0.5; b2[50] = 0.5
    r2 = coloc_abf(b1, s1, b2, s2)
    assert r2["PP_H4"] > 0.9, r2
    # strong signal only in trait1 -> H1 dominant
    b1x = b1.copy(); b2x = rng.normal(0, 0.05, n)
    r3 = coloc_abf(b1x, s1, b2x, s2)
    assert r3["PP_H1"] > 0.8, r3
    return True
