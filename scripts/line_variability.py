"""Hydrogen and helium line strengths of GALEX J085307.3-002204 (Gaia DR3 3074374016095070464, SDSS-V sdss_id 74320911) in every
available spectrum, their variability statistics, and the ZTF period.

Spectra:
- SDSS-V DR20 mwmVisit (Astra 0.8.1) visits. For visits with in_stack = True the XCSAO velocity shift is removed
  (wavelength = grid x (1 + v_xcsao/c)); visits with in_stack = False are listed but not used. Pixels with ivar <= 0 or
  non-finite flux get zero weight.
- SDSS and DESI DR1 spectra within 3 arcsec from SPARCL (NOIRLab Astro Data Lab), native wavelengths.
Equivalent widths (pseudo-EW, Angstrom; absorption positive): for each line a weighted linear continuum through the sidebands
(3 rounds of 3-sigma clipping), EW = sum (1 - f/c) dlambda over the line window; the error combines the flux ivar and the
continuum-fit covariance (inflated by the sideband chi2_r when > 1). Windows and sidebands (vacuum, Angstrom) are in LINES; the
Na D region 5886-5900 A is masked inside the He I 5876 window.
Calibrated errors (SDSS-V visits): the formal error times a factor per line, interpolated in log(error), from
data/line_variability_calibration.csv (group H: visit residuals of all SnowWhite class DA objects with at least 2 visits in the
screen described on docs/line_variability.md). Indices: Z_H = sum of r/sigma over H-alpha, H-beta, H-gamma divided by sqrt(N)
(r = residual about the line's weighted mean), Z_He the same over He I 4471, 4922, 5876, 6678, Z_A = (Z_He - Z_H)/sqrt(2); each
index is divided by its robust width in the S/N bin of the visit (data/line_variability_index_calibration.csv, group H).
p_ZH, p_ZHe, p_ZA: chi2 of the index across visits (sum of squares, dof n - 1). Ratio R = EW(He 4471 + He 5876)/EW(H-beta +
H-alpha) per visit and the chi2 of R about its weighted mean. r_HHe: correlation of Z_H and Z_He over the visits.
ZTF (data/ztf_3074374016095070464.csv, IRSA light-curve service, 1.5 arcsec): catflags = 0, magerr < 0.2; flux normalised by
its median per filter, field, CCD and quadrant; 5-sigma clip; times mjd + exptime/2 (UTC, exposure start) converted to BJD_TDB.
Lomb-Scargle over 0.1-50 c/d (g, r and joint), Baluev false-alarm probability; 200 shuffles of the fluxes within each
filter/field/CCD/quadrant group (seed 1); spectral-window power at the peak; per-season power at the peak; sinusoid fit at the
refined frequency (uncertainty: half-range of the fine grid with chi2 <= chi2_min + chi2_r). Neighbours: ZTF g light curves of
the other ZTF objects within 90 arcsec with at least 100 points (data/ztf_3074374016095070464_neighbours_g.csv), Lomb-Scargle
power at the target frequency. Phases of the spectra on this ephemeris: phase 0 = flux maximum; SDSS-V visits from tai_beg to
tai_end; SDSS and DESI at their dateobs (phase uncertainty propagated from the frequency uncertainty).
Usage: python line_variability.py  -> ../tables/line_variability_3074374016095070464.csv, ../tables/line_variability.csv"""
import os, json, time, numpy as np, pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from scipy import stats
from sdssv import fetch, C

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(R, "data")
GAIA, SID = "3074374016095070464", "74320911"
SRC = pd.read_csv(os.path.join(D, "line_variability_sources.csv"), dtype={"gaia_dr3": str, "sdss_id": str}).set_index("gaia_dr3").loc[GAIA]
C0 = SkyCoord(SRC.ra_deg * u.deg, SRC.dec_deg * u.deg)
PALOMAR = EarthLocation.from_geodetic(lon=-116.8597 * u.deg, lat=33.3563 * u.deg, height=1712 * u.m)
APO = EarthLocation.from_geodetic(lon=-105.8203 * u.deg, lat=32.7803 * u.deg, height=2788 * u.m)
KPNO = EarthLocation.from_geodetic(lon=-111.5997 * u.deg, lat=31.9634 * u.deg, height=2120 * u.m)
LINES = {  # name: (line window, sidebands, masks inside the window)
    "Ha": ((6535, 6595), [(6470, 6510), (6612, 6648)], []),
    "Hb": ((4840, 4885), [(4780, 4828), (4960, 4995)], []),
    "Hg": ((4322, 4360), [(4280, 4305), (4410, 4440)], []),
    "He4471": ((4455, 4490), [(4410, 4440), (4500, 4540)], []),
    "He4922": ((4908, 4938), [(4890, 4904), (4960, 4995)], []),
    "He5876": ((5856, 5900), [(5800, 5845), (5910, 5960)], [(5886, 5900)]),
    "He6678": ((6665, 6692), [(6645, 6660), (6700, 6730)], []),
}
HL, HEL = ["Ha", "Hb", "Hg"], ["He4471", "He4922", "He5876", "He6678"]


def ew_line(w, f, iv, win, sides, masks):
    s = np.zeros(len(w), bool)
    for a, b in sides:
        s |= (w > a) & (w < b)
    s &= iv > 0
    if s.sum() < 8:
        return np.nan, np.nan
    x = w[s] - win[0]; y = f[s]; wt = iv[s]
    for _ in range(3):
        X = np.vstack([np.ones_like(x), x]).T; WX = X * wt[:, None]; cov = np.linalg.inv(X.T @ WX)
        p = cov @ (WX.T @ y); r = (y - X @ p) * np.sqrt(wt); k = np.abs(r) < 3 * max(1.0, np.std(r))
        if k.all() or k.sum() < 8:
            break
        x, y, wt = x[k], y[k], wt[k]
    cov = cov * max(1.0, np.sum(r[k] ** 2) / max(1, k.sum() - 2))
    m = (w >= win[0]) & (w <= win[1])
    for a, b in masks:
        m &= ~((w > a) & (w < b))
    m &= iv > 0
    if m.sum() < 0.6 * ((w >= win[0]) & (w <= win[1])).sum():
        return np.nan, np.nan
    xm = w[m] - win[0]; cm = p[0] + p[1] * xm
    if np.any(cm <= 0):
        return np.nan, np.nan
    dl = np.gradient(w)[m]; ew = np.sum((1 - f[m] / cm) * dl)
    g = np.array([np.sum(f[m] / cm ** 2 * dl), np.sum(f[m] / cm ** 2 * dl * xm)])
    return ew, np.sqrt(np.sum(dl ** 2 / (iv[m] * cm ** 2)) + g @ cov @ g)


def bjd(t, loc):
    t = Time(t, location=loc); return (t.tdb + t.light_travel_time(C0)).jd


def sdssv_rows():
    out = []
    with fits.open(fetch(SID, "visit")) as h:
        for i in (1, 2):
            if h[i].data is None or len(h[i].data) == 0:
                continue
            hd = h[i].header; wg = 10 ** (hd["CRVAL"] + hd["CDELT"] * np.arange(hd["NPIXELS"]))
            for r in h[i].data:
                ins = bool(r["in_stack"]); v = float(r["xcsao_v_rad"]); use = ins and np.isfinite(v)
                w = wg * (1 + v / C) if use else wg; f = np.array(r["flux"], float); iv = np.array(r["ivar"], float)
                iv[~np.isfinite(f) | ~np.isfinite(iv)] = 0; f[iv <= 0] = 0
                tb, te = (Time(float(r[k]) / 86400.0, format="mjd", scale="tai") for k in ("tai_beg", "tai_end"))
                d = dict(dataset="SDSS-V DR20 visit", identifier=f"mwmVisit {SID} MJD {int(r['mjd'])}", mjd=int(r["mjd"]), fieldid=int(r["fieldid"]),
                         bjd_start=bjd(tb, APO), bjd_end=bjd(te, APO), snr=round(float(r["snr"]), 1), in_stack=ins, used=use, _w=w, _f=f, _iv=iv)
                for n, (win, sides, masks) in LINES.items():
                    d[f"ew_{n}"], d[f"e_{n}"] = ew_line(w, f, iv, win, sides, masks)
                out.append(d)
    return out


def sparcl_rows(r_as=3.0):
    from sparcl.client import SparclClient
    for k in range(5):  # the client's connection check times out often; retry
        try:
            cl = SparclClient(); break
        except Exception:
            if k == 4:
                raise
            time.sleep(10)
    d = r_as / 3600; ra, dec = SRC.ra_deg, SRC.dec_deg
    f = cl.find(outfields=["sparcl_id"], constraints={"ra": [ra - d / np.cos(np.radians(dec)), ra + d / np.cos(np.radians(dec))], "dec": [dec - d, dec + d]}, limit=50)
    ids = [x["sparcl_id"] for x in f.records]; out = []
    if not ids:
        return out
    for x in cl.retrieve(uuid_list=ids, include=["flux", "ivar", "wavelength", "data_release", "specid", "dateobs", "specprimary"]).records:
        dt = x["dateobs"]; dt = json.loads(dt) if isinstance(dt, str) and dt.startswith("[") else dt; dt = str(dt[0] if isinstance(dt, (list, tuple)) else dt)
        t = Time(dt[:19].replace(" ", "T"), scale="utc"); w = np.array(x["wavelength"], float); fl = np.array(x["flux"], float); iv = np.array(x["ivar"], float)
        iv[~np.isfinite(fl) | ~np.isfinite(iv)] = 0; fl[iv <= 0] = 0
        rel = str(x["data_release"]); loc = KPNO if rel.startswith("DESI") else APO
        row = dict(dataset=rel, identifier=f"specid {x['specid']}" + ("" if x.get("specprimary", True) else " (not primary)"), mjd=round(t.mjd, 3),
                   bjd_start=bjd(t, loc), bjd_end=np.nan, snr=np.nan, in_stack=np.nan, used=np.nan, _w=w, _f=fl, _iv=iv)
        for n, (win, sides, masks) in LINES.items():
            row[f"ew_{n}"], row[f"e_{n}"] = ew_line(w, fl, iv, win, sides, masks)
        out.append(row)
    return out


def calibrate(V):
    """Calibrated errors and the variability indices of the used SDSS-V visits (group H calibration)."""
    Cal = pd.read_csv(os.path.join(D, "line_variability_calibration.csv")); IDX = pd.read_csv(os.path.join(D, "line_variability_index_calibration.csv"))
    for n in LINES:
        c = Cal[(Cal.grp == "H") & (Cal.line == n)].sort_values("k"); e0 = V[f"e_{n}"].values
        V[f"e_cal_{n}"] = e0 * np.interp(np.log(e0), c.log_sigma.values, c.factor.values)
        e, s = V[f"ew_{n}"].values, V[f"e_cal_{n}"].values; w = 1 / s ** 2; m = np.sum(w * e) / np.sum(w); V[f"zc_{n}"] = (e - m) / s
    sb = np.concatenate([IDX[(IDX["index"] == "ZA") & (IDX.grp == "H")].sort_values("snrbin").snr_lo.values, [np.inf]]); sb[0] = -np.inf
    V["snrbin"] = np.clip(np.searchsorted(sb, V.snr.values) - 1, 0, 5)
    CH = V[[f"zc_{n}" for n in HL]].sum(1) / np.sqrt(len(HL)); CHe = V[[f"zc_{n}" for n in HEL]].sum(1) / np.sqrt(len(HEL)); CA = (CHe - CH) / np.sqrt(2)
    out = {}
    for k, x in (("ZH", CH), ("ZHe", CHe), ("ZA", CA)):
        f = V.snrbin.map(IDX[(IDX["index"] == k) & (IDX.grp == "H")].set_index("snrbin").robust_std).values; x = x.values / f; V[k] = x
        out[f"chi2_{k}"] = float(np.sum(x ** 2)); out[f"p_{k}"] = float(stats.chi2.sf(np.sum(x ** 2), len(x) - 1))
    out["r_HHe"] = float(np.corrcoef(V.ZH, V.ZHe)[0, 1])
    he = V.ew_He4471 + V.ew_He5876; hh = V.ew_Hb + V.ew_Ha; she = np.hypot(V.e_cal_He4471, V.e_cal_He5876); shh = np.hypot(V.e_cal_Hb, V.e_cal_Ha)
    Rr = he / hh; sR = np.abs(Rr) * np.hypot(she / he, shh / hh); w = 1 / sR ** 2; mu = np.sum(w * Rr) / np.sum(w); c2 = float(np.sum(w * (Rr - mu) ** 2))
    out.update(R_min=float(Rr.min()), R_max=float(Rr.max()), chi2_R=c2, p_R=float(stats.chi2.sf(c2, len(V) - 1)))
    a, b = V.iloc[-2], V.iloc[-1]  # the two consecutive nights on the same field
    out["Hb_change_last_two_sigma"] = float((b.ew_Hb - a.ew_Hb) / np.hypot(a.e_cal_Hb, b.e_cal_Hb)); out["last_two_mjd"] = f"{int(a.mjd)}/{int(b.mjd)}"
    out["last_two_same_field"] = bool(a.fieldid == b.fieldid)
    return out


def ztf_series(path, maxerr=0.2):
    z = pd.read_csv(path); T, F, E, B, K = [], [], [], [], []
    for band in ("zg", "zr"):
        d = z[(z.filtercode == band) & (z.catflags == 0) & (z.magerr < maxerr)]
        f = 10 ** (-0.4 * d.mag.values); ef = f * d.magerr.values / 1.0857
        key = (band + "_" + d.field.astype(str) + "_" + d.ccdid.astype(str) + "_" + d.qid.astype(str)).values
        for k in np.unique(key):
            m = key == k; s = np.median(f[m]); f[m] /= s; ef[m] /= s
        ok = np.abs(f - 1) < 5 * 1.4826 * np.median(np.abs(f - 1))
        t = Time(d.mjd.values[ok] + d.exptime.values[ok] / 2 / 86400.0, format="mjd", scale="utc")
        T += list(bjd(t, PALOMAR)); F += list(f[ok]); E += list(ef[ok]); B += [band] * int(ok.sum()); K += list(key[ok])
    return tuple(map(np.array, (T, F, E, B, K)))


def refine(T, F, E, f0, df):
    fine = np.linspace(f0 - df, f0 + df, 4001); w = 1 / E ** 2; chi = []
    for f in fine:
        X = np.vstack([np.ones(len(T)), np.sin(2 * np.pi * f * T), np.cos(2 * np.pi * f * T)]).T
        p = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * F)); chi.append(np.sum(w * (F - X @ p) ** 2))
    chi = np.array(chi); i = int(np.argmin(chi)); fb = fine[i]; chir = chi[i] / (len(T) - 3); ok = fine[chi <= chi[i] + chir]
    X = np.vstack([np.ones(len(T)), np.sin(2 * np.pi * fb * T), np.cos(2 * np.pi * fb * T)]).T; cov = np.linalg.inv(X.T @ (X * w[:, None])) * max(chir, 1)
    p = cov @ (X.T @ (w * F)) / max(chir, 1); amp = np.hypot(p[1], p[2]); e_amp = np.sqrt(cov[1, 1] * (p[1] / amp) ** 2 + cov[2, 2] * (p[2] / amp) ** 2)
    e_f = max((ok.max() - ok.min()) / 2, (fine[1] - fine[0]))
    t0 = (0.25 - np.arctan2(p[2], p[1]) / (2 * np.pi)) / fb; tmax = t0 + np.ceil((2458000.0 - t0) * fb) / fb  # first maximum after BJD 2458000
    return fb, e_f, amp, e_amp, tmax, chir


def ztf_period():
    T, F, E, B, K = ztf_series(os.path.join(D, f"ztf_{GAIA}.csv")); fr = np.linspace(0.1, 50, 500000); out = {}
    for lab, m in (("g", B == "zg"), ("r", B == "zr"), ("joint", np.ones(len(T), bool))):
        L = LombScargle(T[m], F[m], E[m]); p = L.power(fr, method="fast"); i = int(np.argmax(p))
        out[f"ztf_{lab}_n"] = int(m.sum()); out[f"ztf_{lab}_peak_cd"] = float(fr[i]); out[f"ztf_{lab}_power"] = float(p[i])
        out[f"ztf_{lab}_fap"] = float(L.false_alarm_probability(p[i], minimum_frequency=0.1, maximum_frequency=50, method="baluev"))
    f0 = out["ztf_joint_peak_cd"]
    for lab in ("g", "r"):
        m = B == ("zg" if lab == "g" else "zr"); L = LombScargle(T[m], F[m], E[m]); pw = float(L.power(f0))
        out[f"ztf_{lab}_power_at_joint"] = pw
        out[f"ztf_{lab}_fap_at_joint"] = float(L.false_alarm_probability(pw, minimum_frequency=0.1, maximum_frequency=50, method="baluev"))
    W = LombScargle(T, np.ones(len(T)), fit_mean=False, center_data=False).power(np.array([f0])); out["window_power_at_peak"] = float(W[0])
    rng = np.random.default_rng(1); pobs = out["ztf_joint_power"]; mx = []
    for s in range(200):
        Fs = F.copy()
        for k in np.unique(K):
            m = np.where(K == k)[0]; Fs[m] = Fs[rng.permutation(m)]
        mx.append(LombScargle(T, Fs, E).power(fr, method="fast").max())
    out["shuffles_n"] = 200; out["shuffles_ge_observed"] = int(np.sum(np.array(mx) >= pobs)); out["shuffles_max_power"] = float(np.max(mx))
    yr = np.floor(Time(T, format="jd").decimalyear).astype(int)
    out["season_power_at_peak"] = "; ".join(f"{y}: {LombScargle(T[yr == y], F[yr == y], E[yr == y]).power(f0):.3f} (n = {int((yr == y).sum())})" for y in np.unique(yr) if (yr == y).sum() >= 40)
    fb, e_f, amp, e_amp, tmax, chir = refine(T, F - 1, E, f0, 0.002)
    out.update(frequency_cd=fb, e_frequency_cd=e_f, period_h=24 / fb, e_period_h=24 * e_f / fb ** 2, semi_amplitude_frac=amp, e_semi_amplitude_frac=e_amp, t_max_bjd=tmax, chi2_r=chir,
               bjd_first=float(T.min()), bjd_last=float(T.max()))
    nb = pd.read_csv(os.path.join(D, f"ztf_{GAIA}_neighbours_g.csv")); tgt = pd.read_csv(os.path.join(D, f"ztf_{GAIA}.csv")).oid.unique(); pw = []
    for oid, d in nb[(nb.catflags == 0) & (nb.filtercode == "zg")].groupby("oid"):
        if oid in tgt or len(d) < 100:
            continue
        f = 10 ** (-0.4 * d.mag.values); t = bjd(Time(d.mjd.values + d.exptime.values / 2 / 86400.0, format="mjd", scale="utc"), PALOMAR)
        pw.append(float(LombScargle(t, f, f * d.magerr.values / 1.0857).power(fb)))
    out["neighbours_n"] = len(pw); out["neighbours_power_range"] = f"{min(pw):.3f}-{max(pw):.3f}" if pw else ""
    return out


def phase(t, fb, tmax):
    return (t - tmax) * fb % 1.0


if __name__ == "__main__":
    S = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in sdssv_rows() + sparcl_rows()])
    V = S[S.used == True].sort_values("mjd").reset_index(drop=True)
    st = calibrate(V); zt = ztf_period(); fb, tmax, e_f = zt["frequency_cd"], zt["t_max_bjd"], zt["e_frequency_cd"]
    S["phase_start"] = phase(S.bjd_start, fb, tmax); S["phase_end"] = phase(S.bjd_end, fb, tmax)
    S["e_phase"] = np.abs(S.bjd_start - (zt["bjd_first"] + zt["bjd_last"]) / 2) * e_f
    for n in LINES:
        S[f"e_cal_{n}"] = S.identifier.map(dict(zip(V.identifier, V[f"e_cal_{n}"])))
    cols = ["dataset", "identifier", "mjd", "snr", "in_stack", "used", "bjd_start", "bjd_end", "phase_start", "phase_end", "e_phase"]
    for n in LINES:
        cols += [f"ew_{n}", f"e_{n}", f"e_cal_{n}"]
    S = S[cols].copy()
    for c in S.columns:
        if c.startswith(("ew_", "e_", "phase")):
            S[c] = S[c].astype(float).round(3 if "phase" in c else 2)
    S["bjd_start"] = S.bjd_start.round(5); S["bjd_end"] = S.bjd_end.round(5)
    S.to_csv(os.path.join(R, "tables", f"line_variability_{GAIA}.csv"), index=False)
    row = dict(gaia_dr3=GAIA, name=SRC["name"], sdss_id=SID, **{k: SRC[k] for k in SRC.index if k not in ("name", "sdss_id")}, n_visits_used=len(V), **st, **zt)
    pd.DataFrame([row]).to_csv(os.path.join(R, "tables", "line_variability.csv"), index=False)
    pd.set_option("display.width", 250); print(S.drop(columns=[c for c in S.columns if c.startswith("e_cal")]).to_string()); print(json.dumps({k: v for k, v in row.items()}, indent=1, default=str))
