"""Radial velocities of four DA white dwarfs between SDSS-V visits (data/rv_variable_da_sources.csv): per-visit velocities from the
SDSS-V DR20 mwmVisit spectra, per-exposure velocities from the BOSS pipeline spec files, SDSS and DESI DR1 velocities from SPARCL,
and the ZTF periodograms.

Line velocities (H-alpha 6564.61, H-beta 4862.68, H-gamma 4341.69, H-delta 4102.89 A, vacuum):
- Each spectrum is cut to +-4500 km/s around each line (velocity x = c ln(lambda/lambda0)); a linear pseudo-continuum is fitted
  outside |x| > 3300 km/s and divided out for the template only.
- Template per line: a three-Gaussian absorption profile with a common centre times a linear continuum, fitted to the object's own
  coadd of all its spectra of that route (rebuilt after aligning the spectra when they differ by more than 60 km/s).
- Velocity per spectrum and line: on a grid of -1500 to +1500 km/s (10 km/s) the model [1, x, P(x - v), x P(x - v)] is fitted by
  linear least squares over |x| < 2600 km/s (native pixels); parabola refinement of chi2(v); error from the curvature times
  sqrt(max(1, chi2_r)); one 5-sigma residual clip. Lines are kept with dchi2 (no line minus best) > 25, error below 150 km/s
  (200 km/s for BOSS exposures), and not at the grid edge.
- Calibrated line errors: the formal error times a factor interpolated in log(error) from pair-difference tests
  (data/rv_variable_da_calibration.json; store route: all pairs of visits of the 7,203 DA white dwarfs of the screen; BOSS route:
  exposure pairs within one visit of 150 random DA white dwarfs).
- Consensus per spectrum: the largest internally consistent subset of at least two lines (inter-line chi2 p > 1e-3 and spread
  below max(150 km/s, 4 sigma)), highest weight among equal sizes; spectra with two disjoint consistent subsets more than
  100 km/s apart are dropped. Consensus error = 1/sqrt(sum w), times sqrt(chi2_r) of the lines when > 1, then rescaled with the
  consensus-level pair calibration.
Store route (SDSS-V mwmVisit, Astra 0.8.1): visits with in_stack = True only, XCSAO velocity shift removed (lambda = grid x
(1 + v_xcsao/c)); visits with in_stack = False are not used. Field-night zero point: the median residual of all screened objects
on the same field and MJD when there are at least 5 (data/rv_variable_da_zero_points.csv, from the screen); otherwise 0 with
an error of 10 km/s. Visits with S/N < 8 are not used.
BOSS route: spec files of run2d v6_2_1 (data.sdss.org, spectra/daily/full), one per visit listed in the mwmVisit file (in_stack
True or False); per exposure (MJD_EXP_* extensions) FLUX, IVAR with AND_MASK = 0, LOGLAM. Sky-line check: Gaussian centroids
of [O I] 5578.89 and 6302.05 A (vacuum) in the SKY column, as a velocity, plus HELIO_RV (sky_offset_kms; near zero when the
heliocentric wavelength scale is right); exposures with |sky_offset_kms| > 20 km/s are dropped. Between visits an error floor (data/rv_variable_da_calibration.json) is added in quadrature.
SDSS and DESI DR1: spectra within 3 arcsec from SPARCL, each measured with its own template (single spectrum).
Per object and route: chi2 of a constant velocity, log10 p, the largest pair difference and its significance.
ZTF (data/ztf_<gaia_dr3>.csv; IRSA, 2.5 arcsec, ZTF DR24): catflags = 0; per filter flux normalised to the median; mjd +
exptime/2 converted to BJD_TDB; Lomb-Scargle of the joint g, r, i fluxes over 0.1-200 c/d with 1-5 c/d +-0.03 c/d masked;
Baluev false-alarm probability of the highest peak.
Usage: python rv_variable_da.py  -> ../tables/rv_variable_da.csv, ../tables/rv_variable_da_velocities.csv"""
import os, json, time, itertools, requests, numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from scipy.optimize import least_squares, curve_fit
from scipy.stats import chi2 as CHI2
from sdssv import fetch, CACHE, C

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(R, "data")
CAL = json.load(open(os.path.join(D, "rv_variable_da_calibration.json")))
ZP = pd.read_csv(os.path.join(D, "rv_variable_da_zero_points.csv"))
LINES = {"Ha": 6564.61, "Hb": 4862.68, "Hg": 4341.69, "Hd": 4102.89}; LN = list(LINES)
VG = np.arange(-1500, 1500.1, 10.0); XW = 2600.0; GRID = np.arange(-4400, 4400.1, 46.0)
BOSS = "https://data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/spectra/daily/full"
SKY = [5578.887, 6302.046]
PALOMAR = EarthLocation.from_geodetic(lon=-116.8597 * u.deg, lat=33.3563 * u.deg, height=1712 * u.m)


# ---- line-profile velocities
def prof(x, p):
    a1, s1, a2, s2, a3, s3 = p
    return a1 * np.exp(-0.5 * (x / s1) ** 2) + a2 * np.exp(-0.5 * (x / s2) ** 2) + a3 * np.exp(-0.5 * (x / s3) ** 2)


def windows(w, f, iv):
    iv = np.where(np.isfinite(f) & np.isfinite(iv), iv, 0); out = {}
    for n, l in LINES.items():
        x = C * np.log(w / l); m = np.abs(x) < 4500; out[n] = (x[m], np.nan_to_num(f[m]), iv[m])
    return out


def norm_spec(x, f, iv):
    o = (np.abs(x) > 3300) & (iv > 0)
    if o.sum() < 6:
        return None
    c = np.polyval(np.polyfit(x[o], f[o], 1, w=np.sqrt(iv[o])), x)
    if np.any(c <= 0) or np.median(c) <= 0:
        return None
    return f / c, iv * c ** 2


def coadd(V, L, vshift):
    num = np.zeros_like(GRID); den = np.zeros_like(GRID)
    for d, vs in zip(V, vshift):
        x, f, iv = d[L]; r = norm_spec(x, f, iv)
        if r is None:
            continue
        fn, ivn = r; xx = x - vs; ok = ivn > 0
        if ok.sum() < 20:
            continue
        fi = np.interp(GRID, xx[ok], fn[ok], left=np.nan, right=np.nan); ii = np.interp(GRID, xx[ok], ivn[ok], left=0, right=0)
        g = np.isfinite(fi) & (ii > 0); num[g] += fi[g] * ii[g] / 1.5; den[g] += ii[g] / 1.5
    ok = den > 0
    return np.where(ok, num / np.where(ok, den, 1), np.nan), den


def fit_template(fc, den, warm=None):
    ok = np.isfinite(fc) & (den > 0) & (np.abs(GRID) < 4300)
    if ok.sum() < 100:
        return None
    x, y, w = GRID[ok], fc[ok], np.sqrt(den[ok])
    res = lambda q: (y - (q[0] + q[1] * x * 1e-4) * (1 - prof(x - q[2], q[3:]))) * w
    starts = [warm] if warm is not None else [[1.0, 0.0, 0.0, 0.15, s1, 0.15, s2, 0.15, s3] for s1, s2, s3 in ((150, 700, 2500), (80, 400, 1500), (250, 1100, 3500))]
    lo = [0.2, -5, -800, 0, 30, 0, 250, 0, 900]; hi = [5, 5, 800, 1.2, 450, 1.2, 1600, 1.2, 6000]; best = None
    for q0 in starts:
        q0 = np.clip(q0, np.array(lo) + 1e-6, np.array(hi) - 1e-6)
        try:
            r = least_squares(res, q0, bounds=(lo, hi), max_nfev=400)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return None if best is None else dict(p=best.x[3:], q=best.x)


def spectrum_rv(x, f, iv, p):
    m = (np.abs(x) < XW) & (iv > 0) & np.isfinite(f)
    if m.sum() < 30:
        return None
    x, f, iv = x[m], f[m], iv[m]; keep = np.ones(len(x), bool); out = None
    for it in range(2):
        xx, ff, ww = x[keep], f[keep], iv[keep]; P = prof(xx[None, :] - VG[:, None], p); one = np.ones_like(xx)
        B = np.stack([np.broadcast_to(one, P.shape), np.broadcast_to(xx * 1e-3, P.shape), P, P * xx * 1e-3], axis=1)
        A = np.einsum("vkn,vln,n->vkl", B, B, ww); bb = np.einsum("vkn,n,n->vk", B, ww, ff)
        coef = np.linalg.solve(A, bb[..., None])[..., 0]; model = np.einsum("vk,vkn->vn", coef, B); chi = np.sum((ff[None, :] - model) ** 2 * ww[None, :], axis=1)
        B0 = np.stack([one, xx * 1e-3]); c0 = np.linalg.solve((B0 * ww) @ B0.T, (B0 * ww) @ ff); chi0 = np.sum((ff - c0 @ B0) ** 2 * ww)
        i = int(np.argmin(chi)); n = len(xx)
        if i == 0 or i == len(VG) - 1:
            return dict(v=np.nan, e=np.nan, dchi2=chi0 - chi[i], edge=1)
        a, b, _ = np.polyfit(VG[i - 1:i + 2], chi[i - 1:i + 2], 2)
        if a <= 0:
            return dict(v=np.nan, e=np.nan, dchi2=chi0 - chi[i], edge=0)
        cr = chi[i] / max(n - 5, 1); out = dict(v=-b / (2 * a), e=np.sqrt(1 / a) * np.sqrt(max(cr, 1.0)), dchi2=chi0 - chi[i], edge=0)
        if it == 0:
            bad = np.abs((ff - model[i]) * np.sqrt(ww)) > 5
            if bad.sum() == 0:
                break
            idx = np.where(keep)[0]; keep[idx[bad]] = False
            if keep.sum() < 30:
                break
    return out


def measure(V):
    n = len(V); vs = np.zeros(n); tmpl = {}
    for it in range(2):
        res = {}; prev = tmpl; tmpl = {}
        for L in LN:
            t = fit_template(*coadd(V, L, vs), warm=(prev[L]["q"] if prev.get(L) is not None else None)); tmpl[L] = t
            if t is not None:
                for j, d in enumerate(V):
                    res[(j, L)] = spectrum_rv(*d[L], t["p"])
        new = np.zeros(n)
        for j in range(n):
            vv = [(r["v"], r["e"]) for L in LN if (r := res.get((j, L))) is not None and np.isfinite(r["v"]) and r["dchi2"] > 25]
            if vv:
                a = np.array(vv); w = 1 / a[:, 1] ** 2; new[j] = np.sum(w * a[:, 0]) / w.sum()
        if it == 0 and np.all(np.abs(new - np.median(new)) < 60):
            break
        vs = new
    return res


# ---- calibration and consensus
def apply_scale(e, bins):
    B = np.array(bins, float); r = B[:, 1] / B[:, 0]; w = np.maximum.accumulate(r * B[:, 0]); r = w / B[:, 0]
    return e * np.interp(np.log(np.sqrt(2) * e), np.log(B[:, 0]), r)


def consensus(g):
    """g: DataFrame indexed by line with columns v, ec. Returns (velocity, error, lines) or None."""
    avail = [L for L in LN if L in g.index]
    for k in range(len(avail), 1, -1):
        cands = []
        for sub in itertools.combinations(avail, k):
            v = g.loc[list(sub), "v"].values; e = g.loc[list(sub), "ec"].values; w = 1 / e ** 2; mu = np.sum(w * v) / w.sum()
            if CHI2.sf(np.sum(w * (v - mu) ** 2), k - 1) > 1e-3 and np.ptp(v) < max(150, 4 * np.sqrt(e.max() ** 2 + e.min() ** 2)):
                cands.append((w.sum(), list(sub), mu))
        if cands:
            cands.sort(key=lambda z: -z[0]); best = cands[0]
            if any(not set(c[1]) & set(best[1]) and abs(c[2] - best[2]) > 100 for c in cands[1:]):
                return None
            v = g.loc[best[1], "v"].values; e = g.loc[best[1], "ec"].values; w = 1 / e ** 2; mu = np.sum(w * v) / w.sum()
            c2 = np.sum(w * (v - mu) ** 2)
            return mu, 1 / np.sqrt(w.sum()), np.sqrt(max(c2 / (k - 1), 1.0)), "".join(best[1])
    return None


def line_table(res, n, bins, emax):
    rows = []
    for j in range(n):
        for L in LN:
            r = res.get((j, L))
            if r is None or not np.isfinite(r["v"]) or r["edge"] or r["dchi2"] <= 25 or not (0 < r["e"] < emax):
                continue
            rows.append(dict(j=j, line=L, v=r["v"], e=r["e"], ec=float(apply_scale(np.array([r["e"]]), bins[L])[0])))
    return pd.DataFrame(rows, columns=["j", "line", "v", "e", "ec"])


def stats(v, e):
    w = 1 / e ** 2; mu = np.sum(w * v) / w.sum(); c2 = float(np.sum(w * (v - mu) ** 2)); n = len(v)
    dv = np.abs(v[:, None] - v[None, :]); sg = dv / np.hypot(e[:, None], e[None, :]); i, j = np.unravel_index(np.argmax(sg), sg.shape)
    return dict(n=n, v_mean=mu, chi2=c2, logp=float(np.log10(max(CHI2.sf(c2, n - 1), 1e-300))), dv_pair=float(dv[i, j]), sig_pair=float(sg[i, j]))


# ---- store route
def store_route(sdss_id):
    meta, V = [], []
    with fits.open(fetch(sdss_id, "visit")) as h:
        for i in (1, 2):
            if h[i].data is None or len(h[i].data) == 0:
                continue
            hd = h[i].header; wg = 10 ** (hd["CRVAL"] + hd["CDELT"] * np.arange(hd["NPIXELS"])); tel = h[i].name.split("/")[1]
            for r in h[i].data:
                ins = bool(r["in_stack"]); v = float(r["xcsao_v_rad"])
                m = dict(tel=tel, mjd=int(r["mjd"]), fieldid=int(r["fieldid"]), catalogid=int(r["catalogid"]), in_stack=ins, xcsao_v=v, snr=float(r["snr"]),
                         n_exp=int(r["n_exp"]), mjd_mid=(int(r["tai_beg"]) + int(r["tai_end"])) / 2 / 86400)
                meta.append(m)
                if ins and np.isfinite(v):
                    V.append((m, windows(wg * (1 + v / C), np.asarray(r["flux"], float), np.asarray(r["ivar"], float))))
    res = measure([x[1] for x in V]); Lt = line_table(res, len(V), CAL["store_line_bins"], 150); rows = []
    for j, (m, _) in enumerate(V):
        g = Lt[Lt.j == j].set_index("line"); rec = dict(route="SDSS-V visit (store)", mjd=m["mjd"], fieldid=m["fieldid"], mjd_mid_tai=round(m["mjd_mid"], 5), snr=round(m["snr"], 1),
                                                       in_stack=m["in_stack"], xcsao_v_kms=round(m["xcsao_v"], 1), n_exp=m["n_exp"])
        for L in LN:
            rec[f"v_{L}"] = g.loc[L, "v"] if L in g.index else np.nan; rec[f"e_{L}"] = g.loc[L, "ec"] if L in g.index else np.nan
        cs = consensus(g) if m["snr"] >= 8 else None
        if cs is not None:
            mu, e0, s, lines = cs; z = ZP[(ZP.fieldid == m["fieldid"]) & (ZP.mjd == m["mjd"])]
            zp = float(z.zp.iloc[0]) if len(z) and z.n_zp.iloc[0] >= 5 else 0.0
            zpe = float(z.zp_e.iloc[0]) if len(z) and z.n_zp.iloc[0] >= 5 and np.isfinite(z.zp_e.iloc[0]) else CAL["store_zp_e_without_zero_point"]
            ez0 = np.hypot(e0 * s, zpe)
            rec.update(lines=lines, v_consensus=mu, zero_point=zp, e_zero_point=zpe, v=mu - zp, e=float(apply_scale(np.array([ez0]), CAL["store_consensus_bins"])[0]))
        rows.append(rec)
    return meta, rows


# ---- BOSS route
def boss_file(fid, mjd, cat):
    f = f"{int(fid):06d}"; name = f"spec-{f}-{int(mjd)}-{int(cat)}.fits"; p = os.path.join(CACHE, name); os.makedirs(CACHE, exist_ok=True)
    if os.path.exists(p) and os.path.getsize(p) > 100000:
        return p
    for k in range(3):
        try:
            q = requests.get(f"{BOSS}/{f[:3]}XXX/{f}/{int(mjd)}/{name}", timeout=180)
            if q.status_code == 200 and q.content[:6] == b"SIMPLE" and len(q.content) > 100000:
                open(p, "wb").write(q.content); return p
        except requests.RequestException:
            pass
        time.sleep(5)
    return None


def _g(x, a, mu, s, c0, c1):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + c0 + c1 * (x - mu)


def sky_offset(w, sky, hrv):
    out = []
    for l in SKY:
        m = np.abs(w - l) < 8
        if m.sum() < 7:
            continue
        x, y = w[m], sky[m]
        try:
            p, _ = curve_fit(_g, x, y, p0=[y.max() - np.median(y), x[np.argmax(y)], 1.2, np.median(y), 0.0], maxfev=4000)
            if p[0] > 0 and 0.5 < p[2] < 3 and abs(p[1] - l) < 3:
                out.append(C * (p[1] / l - 1))
        except RuntimeError:
            pass
    return (np.mean(out) + hrv) if out else np.nan


def boss_route(meta):
    E, V = [], []
    for m in meta:
        p = boss_file(m["fieldid"], m["mjd"], m["catalogid"])
        if p is None:
            print("BOSS spec file not retrieved:", m["fieldid"], m["mjd"], m["catalogid"]); continue
        with fits.open(p) as h:
            for x in h:
                if not x.name.startswith("MJD_EXP"):
                    continue
                d = x.data; hd = x.header; w = 10 ** d["LOGLAM"].astype(float); hrv = float(hd.get("HELIO_RV", np.nan))
                V.append(windows(w, d["FLUX"].astype(float), d["IVAR"].astype(float) * (d["AND_MASK"] == 0)))
                E.append(dict(route="BOSS exposure", mjd=m["mjd"], fieldid=m["fieldid"], exposure=x.name, in_stack=m["in_stack"],
                              mjd_mid_tai=round((hd["TAI-BEG"] + hd["EXPTIME"] / 2) / 86400, 5), helio_rv_kms=round(hrv, 2),
                              sky_offset_kms=round(sky_offset(w, d["SKY"].astype(float), hrv), 1)))
    res = measure(V); Lt = line_table(res, len(V), CAL["boss_line_bins"], 200)
    for j, rec in enumerate(E):
        g = Lt[Lt.j == j].set_index("line")
        for L in LN:
            rec[f"v_{L}"] = g.loc[L, "v"] if L in g.index else np.nan; rec[f"e_{L}"] = g.loc[L, "ec"] if L in g.index else np.nan
        if abs(rec["sky_offset_kms"]) > 20:
            rec["lines"] = "dropped (sky-line offset)"; continue
        cs = consensus(g)
        if cs is not None:
            mu, e0, s, lines = cs; e = float(apply_scale(np.array([e0]), CAL["boss_consensus_within_bins"])[0])
            rec.update(lines=lines, v_consensus=mu, v=mu, e_within_visit=e, e=float(np.hypot(e, CAL["boss_visit_floor_kms"])))
    return E


# ---- SDSS / DESI (SPARCL)
def sparcl_rv(ra, dec, r_as=3.0):
    from sparcl.client import SparclClient
    for k in range(5):  # the client's connection check times out often; retry
        try:
            cl = SparclClient(); break
        except Exception:
            if k == 4:
                raise
            time.sleep(10)
    d = r_as / 3600
    f = cl.find(outfields=["sparcl_id"], constraints={"ra": [ra - d / np.cos(np.radians(dec)), ra + d / np.cos(np.radians(dec))], "dec": [dec - d, dec + d]}, limit=50)
    ids = [x["sparcl_id"] for x in f.records]; out = []
    if not ids:
        return out
    for x in cl.retrieve(uuid_list=ids, include=["flux", "ivar", "wavelength", "data_release", "specid", "dateobs"]).records:
        dt = x["dateobs"]; dt = json.loads(dt) if isinstance(dt, str) and dt.startswith("[") else dt; dt = str(dt[0] if isinstance(dt, (list, tuple)) else dt)
        res = measure([windows(np.array(x["wavelength"], float), np.array(x["flux"], float), np.array(x["ivar"], float))])
        rec = dict(route=f"{x['data_release']} (SPARCL)", mjd=round(Time(dt[:19].replace(" ", "T")).mjd, 3), exposure=f"specid {x['specid']}")
        for L in LN:
            r = res.get((0, L)); ok = r is not None and np.isfinite(r["v"]) and r["dchi2"] > 25 and not r["edge"]
            rec[f"v_{L}"] = r["v"] if ok else np.nan; rec[f"e_{L}"] = r["e"] if ok else np.nan
        out.append(rec)
    return out


# ---- ZTF
def ztf(gaia, ra, dec):
    p = os.path.join(D, f"ztf_{gaia}.csv")
    if not os.path.exists(p):
        return dict(ztf="no ZTF light curve (outside the ZTF sky)")
    d = pd.read_csv(p); d = d[d.catflags == 0]; c0 = SkyCoord(ra * u.deg, dec * u.deg); tt, yy, ee = [], [], []
    for b in ("zg", "zr", "zi"):
        x = d[d.filtercode == b]
        if len(x) < 20:
            continue
        t = Time(x.mjd.values + x.exptime.values / 2 / 86400, format="mjd", scale="utc", location=PALOMAR)
        tt.append((t.tdb + t.light_travel_time(c0)).jd); yy.append(10 ** (-0.4 * (x.mag.values - np.median(x.mag))) - 1); ee.append(x.magerr.values * 0.921)
    f = np.linspace(0.1, 200, 1000000); f = f[np.min(np.abs(f[:, None] - np.arange(1, 6)[None, :]), axis=1) > 0.03]
    t, y, e = map(np.concatenate, (tt, yy, ee)); L = LombScargle(t, y, e); pw = L.power(f); i = int(np.argmax(pw)); mp = L.model_parameters(f[i])
    return dict(ztf_n=len(t), ztf_peak_cd=round(float(f[i]), 5), ztf_peak_fap=float(L.false_alarm_probability(pw[i], minimum_frequency=0.1, maximum_frequency=200)),
                ztf_peak_semi_amplitude=round(float(np.hypot(mp[1], mp[2])), 4))


if __name__ == "__main__":
    SRC = pd.read_csv(os.path.join(D, "rv_variable_da_sources.csv"), dtype={"gaia_dr3": str, "sdss_id": str})
    VEL, OBJ = [], []
    for s in SRC.itertuples():
        meta, sv = store_route(s.sdss_id); be = boss_route(meta); xe = sparcl_rv(s.ra_deg, s.dec_deg)
        rows = pd.DataFrame(sv + be + xe); rows.insert(0, "gaia_dr3", s.gaia_dr3); rows.insert(1, "name", s.name); VEL.append(rows)
        o = dict(gaia_dr3=s.gaia_dr3, name=s.name, sdss_id=s.sdss_id, **{k: getattr(s, k) for k in SRC.columns if k not in ("gaia_dr3", "name", "sdss_id")},
                 n_sdssv_visits=len(meta), n_sdssv_visits_not_in_stack=sum(not m["in_stack"] for m in meta))
        for lab, q in (("store", rows[rows.route.str.startswith("SDSS-V")]), ("boss", rows[rows.route == "BOSS exposure"])):
            q = q[np.isfinite(q.get("v", pd.Series(dtype=float)))] if "v" in q else q.iloc[0:0]
            if len(q) >= 2:
                st = stats(q.v.values.astype(float), q.e.values.astype(float)); o.update({f"{lab}_{k}": v for k, v in st.items()})
        o.update(ztf(s.gaia_dr3, s.ra_deg, s.dec_deg)); OBJ.append(o)
    V = pd.concat(VEL, ignore_index=True)
    for c in V.columns:
        if V[c].dtype.kind == "f" and c not in ("mjd", "mjd_mid_tai", "helio_rv_kms"):
            V[c] = V[c].round(1)
    order = ["gaia_dr3", "name", "route", "mjd", "exposure", "fieldid", "mjd_mid_tai", "snr", "in_stack", "n_exp", "xcsao_v_kms", "helio_rv_kms", "sky_offset_kms"] + \
        [f"{a}_{L}" for L in LN for a in ("v", "e")] + ["lines", "v_consensus", "zero_point", "e_zero_point", "e_within_visit", "v", "e"]
    V = V[[c for c in order if c in V.columns]]
    V.to_csv(os.path.join(R, "tables", "rv_variable_da_velocities.csv"), index=False)
    O = pd.DataFrame(OBJ)
    for c in O.columns:
        if O[c].dtype.kind == "f" and not c.endswith(("fap",)):
            O[c] = O[c].round(2)
    O.to_csv(os.path.join(R, "tables", "rv_variable_da.csv"), index=False)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
    print(V[["gaia_dr3", "route", "mjd", "exposure", "snr", "lines", "v", "e", "sky_offset_kms"]].to_string()); print(O.T.to_string())
