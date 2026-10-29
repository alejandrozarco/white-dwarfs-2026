"""Checks of the Ca II triplet emission of the five white dwarfs in gas_disc_white_dwarfs.csv (WD 0856+048, WD J1959+2208,
GALEX J0039-0356, SDSS J2054+1610, WDJ1448+3225): error estimates and continuum choice per spectrum, variability, the 8664 A
line profile, metal lines in the X-shooter UVB spectra of WD 0856+048, the BOSS 2010 and SDSS 2003 spectra of WD 0856+048, and
infrared photometry against white-dwarf model photometry. Spectra are retrieved and normalised as in gas_disc_epochs.py.
Usage: python gas_disc_checks.py [epochs|amplitude|xshooter|profiles|infrared]   (no argument: all)
cont_8600 is the continuum at 8600 A and line_flux = ew_A x cont_8600, in the flux units of each spectrum (1e-17 erg/s/cm2/A for
SDSS, BOSS, SDSS-V and DESI; erg/s/cm2/A for X-shooter); median_flux_4500_5000 is the median flux over 4500-5000 A (same units).
Writes ../tables/gas_disc_epoch_checks.csv, gas_disc_variability.csv, gas_disc_578709631539357440_amplitude_tests.csv,
gas_disc_578709631539357440_xshooter_metals.csv, gas_disc_line_profiles.csv, gas_disc_infrared.csv."""
import sys, os, io, numpy as np, pandas as pd, requests
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit, least_squares
from scipy.stats import chi2 as chi2dist
import gas_disc_epochs as G
from sdssv import CACHE, C

TAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables")
W, CAT = G.W, G.CAT
OBJ = ["578709631539357440", "1827014701883095680", "2527617665632689024", "1764314497240770176", "1283510882895711872"]
LINEFREE = ((8336, 8388), (8742, 8794), (8778, 8830))  # each about the width of one +-900 km/s line window
J0856 = "578709631539357440"


def objects():
    t = pd.read_csv(os.path.join(TAB, "gas_disc_white_dwarfs.csv"), dtype=str).set_index("gaia_dr3")
    return t.loc[OBJ]


def norm(w, f, iv, deg=2, rng=(8330, 8830), excl=1100):
    """gas_disc_epochs.norm with the continuum order, fitting range and line exclusion as parameters (defaults = gas_disc_epochs)."""
    w = np.asarray(w, float); f = np.asarray(f, float); iv = np.asarray(iv, float)
    sel = (w > rng[0]) & (w < rng[1]) & (iv > 0) & np.isfinite(f); use = sel.copy()
    for l in list(CAT) + [8448.7]:
        use &= np.abs(w / l - 1) * C > excl
    for _ in range(3):
        p = np.polyfit(w[use], f[use], deg, w=np.sqrt(iv[use])); P = np.polyval(p, w); use &= np.abs((f - P) * np.sqrt(iv)) < 3
    n = np.interp(W, w[sel], (f / P)[sel], left=np.nan, right=np.nan); v = np.interp(W, w[sel], (iv * P ** 2)[sel], left=0, right=0)
    v *= np.median(np.diff(W)) / np.median(np.diff(w[sel]))
    return n, np.where(np.isfinite(n), v, 0)


def continuum_8600(w, f, iv):
    """Value at 8600 A of the default (gas_disc_epochs) continuum fit, in the flux units of the spectrum."""
    w = np.asarray(w, float); f = np.asarray(f, float); iv = np.asarray(iv, float)
    use = (w > 8330) & (w < 8830) & (iv > 0) & np.isfinite(f)
    for l in list(CAT) + [8448.7]:
        use &= np.abs(w / l - 1) * C > 1100
    for _ in range(3):
        p = np.polyfit(w[use], f[use], 2, w=np.sqrt(iv[use])); use &= np.abs((f - np.polyval(p, w)) * np.sqrt(iv)) < 3
    return float(np.polyval(p, 8600))


def ew(n, v):
    dl = np.gradient(W); win = np.zeros(len(W), bool)
    for l in CAT:
        win |= np.abs(W / l - 1) * C < 900
    ok = win & (v > 0) & np.isfinite(n)
    return float(np.sum(((n - 1) * dl)[ok])), float(np.sqrt(np.sum((dl ** 2 / v)[ok])))


def ew_linefree_err(n, v):
    """sqrt(3 <EW_k^2>) over the three line-free windows: the scatter of a three-line equivalent width from residual structure."""
    dl = np.gradient(W); vals = []
    for lo, hi in LINEFREE:
        m = (W > lo) & (W < hi) & (v > 0) & np.isfinite(n)
        if m.sum() > 30:
            vals.append(float(np.sum(((n - 1) * dl)[m])))
    return float(np.sqrt(3 * np.mean(np.square(vals)))) if vals else np.nan


def template(n):
    return np.where(G.WIN, gaussian_filter1d(np.nan_to_num(n - 1), 1), 0.0)


FIT = np.abs(W - 8590) < 240  # the fitting range of amp_rel_ref in gas_disc_epochs.measure


def amp(n, v, T):
    ok = FIT & (v > 0) & np.isfinite(n)
    return float(np.sum(v[ok] * T[ok] * (n[ok] - 1)) / np.sum(v[ok] * T[ok] ** 2)), float(1 / np.sqrt(np.sum(v[ok] * T[ok] ** 2)))


def spectra(g, r):
    sid = r.sdss_id if isinstance(r.sdss_id, str) and r.sdss_id not in ("", "nan") and g != "1283510882895711872" else "-"
    co, per, sp, xs = G.all_spectra(g, sid, float(r.ra_deg), float(r.dec_deg))
    ref = co if co is not None else max(sp + xs, key=G.snr)
    return ([co] if co is not None else []) + per + sp + xs, ref


def epochs():
    rows = []
    for g, r in objects().iterrows():
        allsp, ref = spectra(g, r)
        T = template(G.norm(ref["w"], ref["f"], ref["iv"])[0])
        for s in allsp:
            n, v = G.norm(s["w"], s["f"], s["iv"]); e, ee = ew(n, v); A, eA = amp(n, v, T); c86 = continuum_8600(s["w"], s["f"], s["iv"])
            rows.append(dict(gaia_dr3=g, name=r["name"], dataset=s["dataset"], identifier=s["identifier"], date_utc=s["date"], mjd=s["mjd"],
                             ew_A=round(e, 2), ew_err_A=round(ee, 2), ew_err_linefree_A=round(ew_linefree_err(n, v), 2),
                             ew_A_cont1=round(ew(*norm(s["w"], s["f"], s["iv"], 1))[0], 2), ew_A_cont3=round(ew(*norm(s["w"], s["f"], s["iv"], 3))[0], 2),
                             amp_rel_ref=round(A, 3), amp_err=round(eA, 3), cont_8600=float(f"{c86:.4g}"), line_flux=float(f"{e * c86:.4g}"),
                             median_flux_4500_5000=float(f"{np.nanmedian(np.asarray(s['f'], float)[(np.asarray(s['w']) > 4500) & (np.asarray(s['w']) < 5000)]):.4g}") if np.asarray(s['w']).min() < 4500 else np.nan))
    t = pd.DataFrame(rows); t.to_csv(os.path.join(TAB, "gas_disc_epoch_checks.csv"), index=False); print(t.to_string())
    out = []
    for g, d in t.groupby("gaia_dr3", sort=False):
        d = d[(d.dataset != "SDSS-V coadd") & (d.ew_A > 3 * d.ew_err_A)]
        o = dict(gaia_dr3=g, name=d.name.iloc[0], n_spectra=len(d), ew_min_A=d.ew_A.min(), ew_max_A=d.ew_A.max())
        for lab, err in (("formal", d.ew_err_A), ("formal_linefree", np.hypot(d.ew_err_A, d.ew_err_linefree_A))):
            w = 1 / err ** 2; m = np.sum(w * d.ew_A) / np.sum(w); c2 = float(np.sum(w * (d.ew_A - m) ** 2))
            o.update({f"chi2_{lab}": round(c2, 1), f"p_{lab}": float(f"{chi2dist.sf(c2, len(d) - 1):.2g}")})
        o["dof"] = len(d) - 1; out.append(o)
    v = pd.DataFrame(out); v.to_csv(os.path.join(TAB, "gas_disc_variability.csv"), index=False); print(v.to_string())


def amplitude():
    """WD 0856+048: scale of the SDSS-V coadd profile in the BOSS 2010 and SDSS 2003 spectra under 24 continuum choices (the
    reference profile re-derived with each), an empirical null from circular shifts of the residuals, and injection-recovery."""
    r = objects().loc[J0856]; allsp, co = spectra(J0856, r); rows = []
    tests = [s for s in allsp if s["dataset"] in ("BOSS-DR17", "SDSS-DR17")]
    for s in tests:
        lab = f"{s['dataset']} {s['date']}"
        for deg in (1, 2, 3, 4):
            for rng in ((8330, 8830), (8250, 8950), (8380, 8780)):
                for ex in (1100, 1400):
                    T = template(norm(co["w"], co["f"], co["iv"], deg, rng, ex)[0])
                    rows.append(dict(spectrum=lab, test="continuum", setting=f"order {deg}, {rng[0]}-{rng[1]} A, lines masked +-{ex} km/s",
                                     value=round(amp(*norm(s["w"], s["f"], s["iv"], deg, rng, ex), T)[0], 3)))
        T = template(G.norm(co["w"], co["f"], co["iv"])[0]); n, v = G.norm(s["w"], s["f"], s["iv"]); A, eA = amp(n, v, T)
        rows.append(dict(spectrum=lab, test="default", setting="amp_rel_ref +- formal error", value=f"{A:.3f} +- {eA:.3f}"))
        reg = np.where((W > 8330) & (W < 8830) & np.isfinite(n))[0]; res, vv = (n - 1)[reg], v[reg]; null = []
        for k in range(40, len(reg) - 40):
            rr = np.full(len(W), np.nan); vx = np.zeros(len(W)); rr[reg] = np.roll(res, k); vx[reg] = np.roll(vv, k); null.append(amp(rr + 1, vx, T)[0])
        null = np.array(null)
        for nm, val in (("shifts", len(null)), ("mean", null.mean()), ("standard deviation", null.std()), ("95th percentile", np.percentile(null, 95)), ("maximum", null.max())):
            rows.append(dict(spectrum=lab, test="circular-shift null", setting=nm, value=round(float(val), 3)))
        prof = np.interp(s["w"], W, T, left=0, right=0)
        for a in (0.5, 1.0):
            n2, v2 = G.norm(s["w"], s["f"] * (1 + a * prof), s["iv"] / (1 + a * prof) ** 2)
            rows.append(dict(spectrum=lab, test="injection", setting=f"coadd profile x {a} added before normalisation; recovered amplitude", value=round(amp(n2, v2, T)[0], 3)))
    t = pd.DataFrame(rows); t.to_csv(os.path.join(TAB, f"gas_disc_{J0856}_amplitude_tests.csv"), index=False); print(t.to_string())


def xshooter():
    """WD 0856+048 X-shooter UVB: Mg II 4481 and Ca II K absorption. EW with a linear continuum from sidebands (error scaled by
    the sideband scatter), the same with three other window/sideband choices, and a Gaussian fit for the velocity."""
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    r = objects().loc[J0856]; ra, dec = float(r.ra_deg), float(r.dec_deg)
    uvb = G.eso_xshooter(ra, dec, arm="UVB"); site = EarthLocation.of_site("paranal"); sc = SkyCoord(ra * u.deg, dec * u.deg)
    lines = {"Mg II 4481": G.air_to_vac(np.array([4481.2]))[0], "Ca II K": G.air_to_vac(np.array([3933.663]))[0]}
    choices = ((2.0, (3, 9)), (1.5, (3, 8)), (2.5, (4, 10)), (2.0, (4, 14)))
    rows = []
    for s in uvb:
        w, f, iv = s["w"], s["f"], s["iv"]; ok = (iv > 0) & np.isfinite(f); dl = np.gradient(w)
        bc = sc.radial_velocity_correction("barycentric", obstime=Time(s["mjd"], format="mjd"), location=site).to(u.km / u.s).value
        for nm, l0 in lines.items():
            vals = []
            for half, sb in choices:
                d = w - l0; side = ok & (np.abs(d) > sb[0]) & (np.abs(d) < sb[1]); win = ok & (np.abs(d) < half)
                p = np.polyfit(d[side], f[side], 1, w=np.sqrt(iv[side])); c = np.polyval(p, d)
                scale = max(1.0, float(np.sqrt(np.mean(((f - c) * np.sqrt(iv))[side] ** 2))))
                vals.append((float(np.sum(((1 - f / c) * dl)[win])), float(np.sqrt(np.sum((dl ** 2 / iv / c ** 2)[win & ok])) * scale)))
            m = ok & (np.abs(w - l0) < 9); x = (w[m] / l0 - 1) * C
            g = lambda x, a, vel, sg, c0, c1: (c0 + c1 * x / 300) * (1 - a * np.exp(-0.5 * ((x - vel) / sg) ** 2))
            pp, cv = curve_fit(g, x, f[m], p0=[0.05, 30, 25, np.median(f[m]), 0], sigma=1 / np.sqrt(iv[m]), maxfev=20000); ep = np.sqrt(np.diag(cv))
            others = [x_[0] for x_ in vals[1:]]
            rows.append(dict(line=nm, vacuum_wavelength_A=round(l0, 2), identifier=s["identifier"], mjd=s["mjd"], ew_A=round(vals[0][0], 3), ew_err_A=round(vals[0][1], 3),
                             ew_other_windows_A=f"{min(others):.3f} to {max(others):.3f}", gauss_depth=round(pp[0], 3), e_gauss_depth=round(ep[0], 3),
                             v_topocentric_kms=round(pp[1], 1), e_v_kms=round(ep[1], 1), barycentric_correction_kms=round(bc, 2), v_barycentric_kms=round(pp[1] + bc, 1)))
    t = pd.DataFrame(rows); t.to_csv(os.path.join(TAB, f"gas_disc_{J0856}_xshooter_metals.csv"), index=False); print(t.to_string())


def profiles():
    """Profile of the isolated Ca II 8664.52 line (normalised as in gas_disc_epochs.py; 40 km/s steps, 1-step smoothing): peak
    velocities and separation and the central minimum relative to the lower peak (given where that minimum is below 0.8, i.e.
    double-peaked), the FWHM and the velocities where the profile crosses 25% of its peak."""
    pick = {J0856: ("SDSS-V coadd", "ESO XSHOOTER VIS"), "1827014701883095680": ("SDSS-V coadd", "SDSS-V visit"), "2527617665632689024": ("SDSS-V coadd",),
            "1764314497240770176": ("SDSS-V coadd",), "1283510882895711872": ("DESI-DR1",)}
    vel = np.arange(-1300, 1301, 40.0); rows = []
    for g, r in objects().iterrows():
        allsp, _ = spectra(g, r)
        for s in allsp:
            if s["dataset"] not in pick[g]:
                continue
            n, v = G.norm(s["w"], s["f"], s["iv"]); vv = (W / CAT[2] - 1) * C; ok = np.isfinite(n) & (v > 0) & (np.abs(vv) < 1500)
            P = gaussian_filter1d(np.interp(vel, vv[ok], n[ok] - 1), 1.0); pk = P.max()
            above = lambda fr: (vel[np.where(P > fr * pk)[0][0]], vel[np.where(P > fr * pk)[0][-1]])
            b, rd = vel < -50, vel > 50; vb, vr = vel[b][np.argmax(P[b])], vel[rd][np.argmax(P[rd])]; low = min(P[b].max(), P[rd].max())
            trough = P[(vel > -100) & (vel < 100)].min() / low if low > 0 else np.nan
            h50, h25 = above(0.5), above(0.25); dp = bool(trough < 0.8)
            rows.append(dict(gaia_dr3=g, name=r["name"], dataset=s["dataset"], identifier=s["identifier"], date_utc=s["date"], peak_over_continuum=round(pk, 2),
                             double_peaked=dp, v_blue_peak_kms=int(vb) if dp else "", v_red_peak_kms=int(vr) if dp else "", peak_separation_kms=int(vr - vb) if dp else "",
                             central_minimum_over_lower_peak=round(float(trough), 2) if dp else "", fwhm_kms=int(h50[1] - h50[0]), v25_lo_kms=int(h25[0]), v25_hi_kms=int(h25[1])))
    t = pd.DataFrame(rows)
    for c in ("v_blue_peak_kms", "v_red_peak_kms", "peak_separation_kms"):
        t[c] = pd.to_numeric(t[c]).astype("Int64")
    t.to_csv(os.path.join(TAB, "gas_disc_line_profiles.csv"), index=False); print(t.to_string())


ZP = dict(J=1594.0, H=1024.0, Ks=666.7, W1=309.54, W2=171.787)  # Jy, Vega (2MASS, WISE)
MODEL = {"DA": "https://www.astro.umontreal.ca/~bergeron/CoolingModels/Tables/Table_DA",
         "DB": "https://www.astro.umontreal.ca/~bergeron/CoolingModels/Tables/Table_DB"}
ATM = {J0856: "DA", "1827014701883095680": "DB", "2527617665632689024": "DA", "1764314497240770176": "DA", "1283510882895711872": "DB"}


def model_colour(atm):
    p = os.path.join(CACHE, f"Table_{atm}")
    if not os.path.exists(p):
        os.makedirs(CACHE, exist_ok=True); open(p, "w").write(requests.get(MODEL[atm], timeout=120).text)
    L = open(p).read().splitlines(); hdr = L[1].replace("log g", "logg").split()
    M = np.array([[float(x) for x in l.split()] for l in L[2:] if len(l.split()) == len(hdr)]); ix = {k: hdr.index(k) for k in ("G3", "W1", "W2", "J", "H", "Ks")}

    def colour(teff, logg, band):
        lg = min(max(logg, 7.0), 9.0); lo = min(np.floor(lg * 2) / 2, 8.5); v = []
        for l in (lo, lo + 0.5):
            q = M[np.isclose(M[:, 1], l)]; q = q[np.argsort(q[:, 0])]; v.append(np.interp(teff, q[:, 0], q[:, ix["G3"]] - q[:, ix[band]]))
        return float(np.interp(lg, [lo, lo + 0.5], v))
    return colour


def irsa(table, cols, ra, dec, r):
    """Cone search of an IRSA catalogue (catwise_2020, fp_psc) through the IRSA TAP service; sorted by separation (arcsec)."""
    q = f"SELECT ra, dec, {cols} FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{r / 3600}))=1"
    t = pd.read_csv(io.StringIO(requests.get("https://irsa.ipac.caltech.edu/TAP/sync", params=dict(QUERY=q, FORMAT="csv"), timeout=300).text))
    t["sep"] = np.hypot((t.ra - ra) * np.cos(np.radians(dec)), t.dec - dec) * 3600
    return t.sort_values("sep")


def ps1(ra, dec, r):
    """Pan-STARRS1 DR2 mean PSF magnitudes (MAST catalogs API) of the nearest source within r arcsec."""
    u_ = ("https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?" + f"ra={ra}&dec={dec}&radius={r / 3600}&nDetections.gt=1&columns="
          + "[raMean,decMean,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag]")
    t = pd.read_csv(io.StringIO(requests.get(u_, timeout=300).text)); t["sep"] = np.hypot((t.raMean - ra) * np.cos(np.radians(dec)), t.decMean - dec) * 3600
    return t.sort_values("sep").iloc[0]


def legacy(ra, dec, rad=8.0):
    d = rad / 3600
    q = ("SELECT ra, dec, type, flux_z, flux_w1, flux_ivar_w1, flux_w2, flux_ivar_w2, fracflux_w1, fracflux_w2 FROM ls_dr10.tractor WHERE "
         f"ra BETWEEN {ra - d / np.cos(np.radians(dec))} AND {ra + d / np.cos(np.radians(dec))} AND dec BETWEEN {dec - d} AND {dec + d}")
    t = pd.read_csv(io.StringIO(requests.get("https://datalab.noirlab.edu/tap/sync", params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q), timeout=300).text))
    if "ra" not in t.columns:
        return t
    t["sep"] = np.hypot((t.ra - ra) * np.cos(np.radians(dec)), t.dec - dec) * 3600
    return t.sort_values("sep")


def bb(lam_um, T):
    nu = 2.99792458e14 / lam_um; return nu ** 3 / np.expm1(6.62607015e-34 * nu / (1.380649e-23 * T)) / 1e40


def infrared():
    """LS DR10 Tractor W1/W2 (forced, deblended) and CatWISE2020 / 2MASS (IRSA) against Montreal model photometry (pure-H for DA, pure-He
    for DB/DBA) at the tabulated Teff and log g, scaled to Gaia G. Model uncertainty taken as 10%. For WD J1959+2208 (blended with
    a Gaia source 0.97 arcsec away) the Pan-STARRS1 + 2MASS + CatWISE fluxes are also fitted with two blackbodies (8% errors)."""
    rows = []
    for g, r in objects().iterrows():
        ra, dec, Gm, teff, logg = float(r.ra_deg), float(r.dec_deg), float(r.G), float(r.teff_K), float(r.logg); col = model_colour(ATM[g])
        mod = {b: ZP[b] * 1e6 * 10 ** (-0.4 * (Gm - col(teff, logg, b))) for b in ZP}
        o = dict(gaia_dr3=g, name=r["name"], model=f"Montreal {ATM[g]} {teff:.0f} K, log g {logg}", W1_model_uJy=round(mod["W1"], 1), W2_model_uJy=round(mod["W2"], 1))
        ls = legacy(ra, dec)
        if "sep" in ls.columns and len(ls) and ls.sep.iloc[0] < 1.5:
            m = ls.iloc[0]
            for b in ("W1", "W2"):
                f_, e_ = m[f"flux_{b.lower()}"] * 3.631, 3.631 / np.sqrt(m[f"flux_ivar_{b.lower()}"])
                o.update({f"ls_{b}_uJy": round(f_, 1), f"ls_e_{b}_uJy": round(e_, 1), f"ls_{b}_ratio": round(f_ / mod[b], 2),
                          f"ls_{b}_excess_sigma": round((f_ - mod[b]) / np.hypot(e_, 0.1 * mod[b]), 1)})
            o["ls_fracflux_w1"] = round(m.fracflux_w1, 2)
            o["ls_neighbours_8arcsec"] = "; ".join(f"{x.sep:.1f} arcsec {x.type} W1 {x.flux_w1 * 3.631:.1f} uJy" for x in ls.iloc[1:].itertuples()) or "none"
        else:
            o["ls_W1_uJy"] = "outside Legacy Surveys DR10" if not len(ls) or "sep" not in ls.columns else "no source within 1.5 arcsec"
        cw = irsa("catwise_2020", "w1mpro_pm, w1sigmpro_pm, w2mpro_pm, w2sigmpro_pm", ra, dec, 3)
        if len(cw):
            c0 = cw.iloc[0]; o["catwise_sep_arcsec"] = round(float(c0.sep), 2)
            for b in ("W1", "W2"):
                f_ = ZP[b] * 1e6 * 10 ** (-0.4 * c0[f"{b.lower()}mpro_pm"]); e_ = f_ * c0[f"{b.lower()}sigmpro_pm"] / 1.0857
                o.update({f"catwise_{b}_uJy": round(f_, 1), f"catwise_e_{b}_uJy": round(e_, 1), f"catwise_{b}_ratio": round(f_ / mod[b], 2)})
        else:
            o["catwise_sep_arcsec"] = "no CatWISE2020 source within 3 arcsec"
        tm = irsa("fp_psc", "j_m, h_m, k_m", ra, dec, 2)
        if len(tm):
            for b, k in (("J", "j_m"), ("H", "h_m"), ("Ks", "k_m")):
                f_ = ZP[b] * 1e6 * 10 ** (-0.4 * tm.iloc[0][k]); o[f"2mass_{b}_ratio"] = round(f_ / mod[b], 2)
        if g == "1827014701883095680":
            o.update(two_blackbodies(ra, dec, tm, cw))
        rows.append(o)
    t = pd.DataFrame(rows); t.to_csv(os.path.join(TAB, "gas_disc_infrared.csv"), index=False); print(t.T.to_string())


def two_blackbodies(ra, dec, tm, cw):
    lam = dict(g=0.481, r=0.617, i=0.752, z=0.866, y=0.962, J=1.235, H=1.662, Ks=2.159, W1=3.368, W2=4.618)
    ps = ps1(ra, dec, 1.5); obs = {b: 3631e6 * 10 ** (-0.4 * ps[f"{b}MeanPSFMag"]) for b in "grizy"}
    obs.update({b: ZP[b] * 1e6 * 10 ** (-0.4 * tm.iloc[0][k]) for b, k in (("J", "j_m"), ("H", "h_m"), ("Ks", "k_m"))})
    obs.update({b: ZP[b] * 1e6 * 10 ** (-0.4 * cw.iloc[0][f"{b.lower()}mpro_pm"]) for b in ("W1", "W2")})
    err = {b: max(0.08 * obs[b], 1.0) for b in obs}
    model = lambda p, l: np.exp(p[1]) * bb(l, p[0]) + np.exp(p[3]) * bb(l, p[2])
    res = lambda p: [(obs[b] - model(p, lam[b])) / err[b] for b in obs]
    best = None
    for T1 in (12000, 15400, 20000):
        for T2 in (2500, 3500, 4500, 6000):
            f = least_squares(res, [T1, np.log(obs["g"] / bb(0.481, T1)) * 0.9, T2, np.log(obs["Ks"] / bb(2.159, T2)) * 0.5], bounds=([10000, -90, 1800, -90], [25000, 90, 8000, 90]))
            c2 = float(np.sum(np.square(f.fun)))
            if best is None or c2 < best[0]:
                best = (c2, f.x)
    c2, p = best; cool = lambda l: np.exp(p[3]) * bb(l, p[2])
    one = least_squares(lambda q: [(obs[b] - np.exp(q[1]) * bb(lam[b], q[0])) / err[b] for b in obs], [15000, np.log(obs["g"] / bb(0.481, 15000))], bounds=([10000, -90], [25000, 90]))
    return dict(twobb_T_hot_K=int(p[0]), twobb_T_cool_K=int(p[2]), twobb_chi2=round(c2, 1), twobb_dof=len(obs) - 4,
                twobb_cool_fraction_W1=round(cool(3.368) / model(p, 3.368), 2), twobb_cool_AB_0p64um=round(-2.5 * np.log10(cool(0.64) / 3.631e9), 2),
                onebb_chi2=round(float(np.sum(np.square(one.fun))), 1))


if __name__ == "__main__":
    steps = dict(epochs=epochs, amplitude=amplitude, xshooter=xshooter, profiles=profiles, infrared=infrared)
    for k in (sys.argv[1:] or list(steps)):
        steps[k]()
