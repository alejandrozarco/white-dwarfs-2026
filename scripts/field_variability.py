"""Zeeman splitting of H-alpha and H-beta in each SDSS-V visit of two magnetic white dwarfs, the change of the field between visits,
and their ZTF periods (data/field_variability_sources.csv).

Spectra: SDSS-V DR20 mwmVisit (Astra 0.8.1) visits with in_stack = True, XCSAO velocity shift removed (sdssv.visits). Each line
is fitted per visit with the three-component model of zeeman_split.py (linear continuum, broad concentric Gaussian, three
Gaussian components with free centres, widths and depths), seeded at the half-separation of a fit to the coadd of the same
visits. B = half-separation of the outer components / (4.67e-13 lambda0^2): 20.13 A/MG at H-alpha, 11.04 A/MG at H-beta.
Per visit: S = (B_Ha + B_Hb)/2 and D = (B_Ha - B_Hb)/2, with error e = sqrt(e_Ha^2 + e_Hb^2)/2. The formal errors are multiplied
by infl = max(1, rms of (D - mean D)/e with n/(n - 1)), i.e. by the scatter of the H-alpha/H-beta disagreement, which measures the
fit systematics; chi2 of S about its weighted mean (dof n - 1) gives p. dB_Ha and dB_Hb: change of each line's field between
the visits with the highest and the lowest S.
ZTF (data/ztf_<gaia_dr3>.csv, IRSA light-curve service, 1.5 arcsec): catflags = 0, magerr < 0.2; flux normalised by its median per
filter, field, CCD and quadrant; 5-sigma clip; mjd + exptime/2 (UTC) converted to BJD_TDB; Lomb-Scargle per filter over
0.1-50 c/d with the Baluev false-alarm probability of the highest peak and the next two peaks (more than 0.01 c/d apart).
Usage: python field_variability.py  -> ../tables/field_variability_visits.csv, ../tables/field_variability.csv"""
import os, numpy as np, pandas as pd
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from scipy import stats
from sdssv import visits, coadd
from zeeman_split import fit, L0, K, WIN

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(R, "data")
PALOMAR = EarthLocation.from_geodetic(lon=-116.8597 * u.deg, lat=33.3563 * u.deg, height=1712 * u.m)


def field(lam, f, iv, k, dl0):
    m = (lam > WIN[k][0]) & (lam < WIN[k][1]); p, Cv = fit(lam[m], f[m], iv[m], L0[k], dl0)
    return (p[10] + p[11]) / 2 / K[k], 0.5 * np.sqrt(Cv[10, 10] + Cv[11, 11] + 2 * Cv[10, 11]) / K[k]


def per_visit(sdss_id):
    vs = [x for x in visits(sdss_id) if x["in_stack"]]; lam, f, iv = coadd(vs); B0 = {}
    for k in ("Ha", "Hb"):
        B0[k] = field(lam, f, iv, k, 6.0 * K[k])[0]
    rows = []
    for v in vs:
        r = dict(mjd=v["mjd"], snr=round(v["snr"], 1), xcsao_v_kms=round(v["xcsao_v"], 1))
        for k in ("Ha", "Hb"):
            b, e = field(v["wave"], v["flux"], v["ivar"], k, B0[k] * K[k]); r[f"B_{k}_MG"] = b; r[f"e_B_{k}_MG"] = e
        rows.append(r)
    return pd.DataFrame(rows), B0


def variability(d):
    S = (d.B_Ha_MG + d.B_Hb_MG) / 2; Dd = (d.B_Ha_MG - d.B_Hb_MG) / 2; e = np.hypot(d.e_B_Ha_MG, d.e_B_Hb_MG) / 2; n = len(d)
    infl = max(1.0, float(np.sqrt(np.mean(((Dd - Dd.mean()) / e) ** 2) * n / (n - 1))))
    w = 1 / (e * infl) ** 2; mu = float(np.sum(w * S) / w.sum()); c2 = float(np.sum(w * (S - mu) ** 2)); i, j = int(np.argmax(S.values)), int(np.argmin(S.values))
    return dict(n_visits=n, B_mean_MG=mu, B_min_MG=float(S.min()), B_max_MG=float(S.max()), dB_range_MG=float(S.max() - S.min()),
                dB_Ha_MG=float(d.B_Ha_MG.values[i] - d.B_Ha_MG.values[j]), dB_Hb_MG=float(d.B_Hb_MG.values[i] - d.B_Hb_MG.values[j]),
                error_inflation=infl, chi2=c2, p=float(stats.chi2.sf(c2, n - 1)), mjd_high=int(d.mjd.values[i]), mjd_low=int(d.mjd.values[j])), S, e * infl


def ztf_ls(gaia, c0):
    z = pd.read_csv(os.path.join(D, f"ztf_{gaia}.csv")); out = {}
    for band in ("zg", "zr"):
        d = z[(z.filtercode == band) & (z.catflags == 0) & (z.magerr < 0.2)]
        f = 10 ** (-0.4 * d.mag.values); ef = f * d.magerr.values / 1.0857
        key = (d.field.astype(str) + "_" + d.ccdid.astype(str) + "_" + d.qid.astype(str)).values
        for k in np.unique(key):
            m = key == k; s = np.median(f[m]); f[m] /= s; ef[m] /= s
        ok = np.abs(f - 1) < 5 * 1.4826 * np.median(np.abs(f - 1))
        t = Time(d.mjd.values[ok] + d.exptime.values[ok] / 2 / 86400.0, format="mjd", scale="utc", location=PALOMAR); t = (t.tdb + t.light_travel_time(c0)).jd
        fr = np.linspace(0.1, 50, 400000); L = LombScargle(t, f[ok], ef[ok]); p = L.power(fr, method="fast"); top = []
        for j in np.argsort(p)[::-1]:
            if all(abs(fr[j] - x) > 0.01 for x, _ in top):
                top.append((fr[j], p[j]))
            if len(top) == 3:
                break
        b = band[1]; amp = float(np.sqrt(2) * np.std(L.model(t, top[0][0]) - np.mean(f[ok])))
        out.update({f"ztf_{b}_n": int(ok.sum()), f"ztf_{b}_period_h": 24 / top[0][0], f"ztf_{b}_power": float(top[0][1]),
                    f"ztf_{b}_fap": float(L.false_alarm_probability(top[0][1], minimum_frequency=0.1, maximum_frequency=50, method="baluev")),
                    f"ztf_{b}_next_peaks_h": "; ".join(f"{24 / a:.4f}" for a, _ in top[1:]), f"ztf_{b}_semi_amplitude_frac": amp})
    return out


if __name__ == "__main__":
    SRC = pd.read_csv(os.path.join(D, "field_variability_sources.csv"), dtype={"gaia_dr3": str, "sdss_id": str})
    VR, OR = [], []
    for s in SRC.itertuples():
        d, B0 = per_visit(s.sdss_id); st, S, eS = variability(d)
        d.insert(0, "gaia_dr3", s.gaia_dr3); d.insert(1, "name", s.name); d["B_mean_MG"] = S.values; d["e_B_mean_MG"] = eS.values; VR.append(d)
        OR.append(dict(gaia_dr3=s.gaia_dr3, name=s.name, sdss_id=s.sdss_id, **{k: getattr(s, k) for k in SRC.columns if k not in ("gaia_dr3", "name", "sdss_id")},
                       B_coadd_Ha_MG=B0["Ha"], B_coadd_Hb_MG=B0["Hb"], **st, **ztf_ls(s.gaia_dr3, SkyCoord(s.ra_deg * u.deg, s.dec_deg * u.deg))))
    V = pd.concat(VR).round(3); O = pd.DataFrame(OR)
    V.to_csv(os.path.join(R, "tables", "field_variability_visits.csv"), index=False); O.to_csv(os.path.join(R, "tables", "field_variability.csv"), index=False)
    pd.set_option("display.width", 250); print(V.to_string()); print(O.T.to_string())
