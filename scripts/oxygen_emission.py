"""O I 7774 and 8446 emission of Gaia DR3 6488421971834040576 (WDJ235008.79-601857.00) in its SDSS-V DR20 BOSS spectra, with line
widths and the 8446/7774 flux ratio of WD J0914+1914 (Gaia DR3 635879254021761024; Gaensicke et al. 2019) measured with the same code.
Vacuum wavelengths; the pipeline's heliocentric frame. O I reference wavelengths: gf-weighted vacuum centroids of the multiplets
(NIST air wavelengths and log gf, Morton 2000 air-to-vacuum), 7775.51 and 8448.80 A. Fluxes in 1e-17 erg/s/cm2.
Usage: python oxygen_emission.py [lines|frames|sed|variability]   (no argument: all, in this order)

lines        BOSS v6_2_1 spec-full files from the SDSS DR20 SAS (field 109553, MJD 60203 and 60206; WD J0914+1914: field 102026,
             MJD 60653). Spectra: each exposure (MJD_EXP_* extensions), each visit coadd (COADD), the inverse-variance stack of the four
             exposures, and the visit coadd of WD J0914+1914. Pixels with AND_MASK != 0 are not used. Per spectrum and line:
             - window EW, the metric of the screen: linear continuum through 7705-7762 / 7788-7845 A (8395-8438 / 8458-8492 A) with
               3-sigma clipping, EW over 7768-7784 A (8441.5-8455.5 A), error from IVAR plus the continuum-level error, s = -EW / error;
             - single-line fit over +-70 A around the line at +50 km/s: linear continuum + one Gaussian (centre within 20 A, sigma 1-20 A), errors
               scaled by sqrt(max(1, chi2_r)); EW = -flux / continuum; sigma deconvolved with the pipeline WDISP; Delta BIC against the
               continuum alone, against continuum + a free multiple of the pipeline sky model, and of two Gaussians (common sigma)
               against one (positive values favour the first model); largest single-pixel share of the window excess;
             - joint fit of both lines (as in frames, with the pipeline sky model of the spectrum).
             Over the four exposures: inverse-variance mean of the single-line centres and their ratio. On the four-exposure stack:
             EWs of other lines with a Gaussian of the joint-fit velocity and width fixed on a polynomial continuum (+-40 A; +-25 A and
             quadratic for [S II] 4068/4076 next to H-delta); emission limit = max(0, -EW) + 3 sigma.
frames       spCFrame-r2 files (red camera, all fibres, IVAR as delivered) of the four exposures (20585, 20701, 20702, 20703; fibre 447)
             and of five exposures in which fibre 447 was elsewhere: 20586 (field 109553, MJD 60203; a sky position 2 arcmin away), 20571
             (field 106186) and 20559 (106406) on MJD 60203, 20695 (106276) and 20690 (106146) on MJD 60206.
             - joint fit of fibre 447 on native pixels: per +-45 A window a linear continuum, a free multiple beta of the fibre's sky model
               (common to both windows), and two Gaussians with a common velocity and width (fluxes free); errors from the Jacobian;
               Delta chi2 against the same model without the Gaussians (4 parameters);
             - O I nightglow in fibre 447's sky model: Gaussian + constant over 7766-7786 A; for 8446 a Gaussian (centre 8447.5-8450.5 A,
               sigma 1.0-2.5 A) + linear term over 8438-8452 A;
             - mask bits other than NOPLUG within +-20 A of 7777 and 8451 A in fibre 447;
             - window z for every fibre: flux minus a clipped linear continuum (25-70 A from the centre) summed over +-15 A, divided by
               its error; centres 7777.0 and 8450.5 A; sky fibres, fibres 446 and 448, and all fibres with both windows above 2.5;
             - per-exposure chance rate: the same window z for every other science and sky fibre at window pairs (lambda, R lambda),
               R = 8448.80 / 7775.51, lambda = 7000-7950 A in 10 A steps (7560-7740 A and +-40 A around 7777 A excluded); fraction of
               pairs with both z at or above fibre 447's values. Expected number in the screen: the night-1 rate times the lowest night-2
               rate (the three night-2 exposures counted once) times the 50,886 SDSS-V objects screened.
sed          Gaia DR3 astrometry (ESA archive TAP) and SDSS-V SnowWhite class and parameters (CAS); Montreal pure-H synthetic photometry
             (Table_DA; Bedard et al. 2020) fitted to data/oxygen_emission_photometry.csv with E(B-V) = 0.015 (SFD, DES DR2 EBV_SFD98) and
             the Gaia DR3 parallax (+0.025 mas zero point) as a Gaussian prior: Teff / log g / parallax grid, interpolation in log Teff
             and log g, 1-sigma ranges from Delta chi2 <= 1. Band subsets: all; UV + Gaia; Gaia + DES + VHS; Gaia only. The same
             photometry at the SnowWhite and Gentile Fusillo et al. (2021) parameters. Gaia DR3 sources within 30 arcsec.
variability  NSC DR2 single-epoch DECam photometry (NOIRLab Data Lab TAP; flags 0, error < 0.2 mag) and SkyMapper DR4 per-image
             photometry (SkyMapper TAP; flags 0, nimaflags 0, error < 0.3 mag): per band N, rms, chi2_r and the 95 % upper limit on a
             sinusoid's semi-amplitude from the excess variance, A = sqrt(2 (max(0, var - <err^2>) + 1.645 sqrt(2/(N-1)) var)).
             ATLAS forced photometry (data/atlas_forced_photometry_6488421971834040576.txt; duJy > 0, err = 0, chi/N < 3, three 5-sigma
             MAD clips per filter, per-unit mean removed): Lomb-Scargle 0.05-50 c/d; a 4-harmonic series at the Gaia DR3 RR Lyrae period
             of Gaia DR3 6488422354087356160 (14.4 arcsec away) fitted and removed; top residual peak (+-0.01 c/d around its harmonics and
             their +-1, +-2 sidereal-day aliases excluded) against the 99th percentile of the top peak of 200 shuffled residual sets.
Writes ../tables/oxygen_emission_exposures.csv (one row per spectrum or frame) and ../tables/oxygen_emission.csv (this object, and WD
J0914+1914 as the comparison). Each step keeps its results in ../data/cache/oxygen_emission_steps.json, and both tables are rebuilt from
that file after every step, so a step can be rerun on its own. The frames step downloads nine files of about 92 MB."""
import os, sys, io, json, time, numpy as np, pandas as pd, requests
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit, least_squares
from scipy.stats import chi2 as chi2dist
from sdssv import CACHE, C, cas

TAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables")
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
STEPS = os.path.join(CACHE, "oxygen_emission_steps.json")
SAS = "https://data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1"
GAIA, NAME = "6488421971834040576", "WDJ235008.79-601857.00"
J0914, J0914_NAME = "635879254021761024", "WD J0914+1914"
SPEC = {GAIA: [(109553, 60203, "63050396476576493"), (109553, 60206, "63050396476576493")], J0914: [(102026, 60653, "63050394819193233")]}
FIBRE = 447
S1_FRAMES = (20585, 20701, 20702, 20703)
CMP_FRAMES = ((109553, 20586), (106186, 20571), (106406, 20559), (106276, 20695), (106146, 20690))
N_SDSSV = 50886  # SDSS-V DR20 objects in the screen
GF21 = {GAIA: (18767, 7.41), J0914: (20389, 8.01)}  # Gentile Fusillo et al. (2021) H-atmosphere Teff, log g (VizieR J/MNRAS/508/3877)
TIC = {GAIA: "2055094754"}
RRL_P = 0.71399638385608877  # Gaia DR3 vari_rrlyrae pf (d) of Gaia DR3 6488422354087356160
EBV, PLX_ZP = 0.015, 0.025


def air2vac(a):
    s2 = (1e4 / np.asarray(a, float)) ** 2
    return a * (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


OI = {"o7774": (air2vac(np.array([7771.944, 7774.166, 7775.388])), 10 ** np.array([0.369, 0.223, 0.002])),
      "o8446": (air2vac(np.array([8446.250, 8446.360, 8446.760])), 10 ** np.array([-0.463, 0.236, 0.014]))}
REF = {k: float(np.sum(l * g) / g.sum()) for k, (l, g) in OI.items()}
LINES = ("o7774", "o8446")
FIT = {k: (v * (1 + 50 / C) - 70, v * (1 + 50 / C) + 70) for k, v in REF.items()}  # +-70 A around the line at +50 km/s
WINDOW = {"o7774": (7776.2, (7768.0, 7784.0), [(7705.0, 7762.0), (7788.0, 7845.0)]), "o8446": (8448.7, (8441.5, 8455.5), [(8395.0, 8438.0), (8458.0, 8492.0)])}
FIXED_SKY = [(5574.9, 5582.9), (5888.6, 5894.6), (5894.6, 5900.6), (6298.0, 6306.0), (6361.5, 6369.5), (6870.0, 6895.0), (7596.0, 7705.0)]
ZWIN = {"o7774": 7777.0, "o8446": 8450.5}  # window-z centres in the frames step
OTHER = {"HeI5876": (5877.25, 40, 1), "HeI6678": (6679.99, 40, 2), "HeII4686": (4687.02, 40, 2), "[SII]4068": (4069.75, 25, 2),
         "[SII]4076": (4077.50, 25, 2), "SI9213": (air2vac(9212.86), 40, 1), "OI9263": (air2vac(9263.0), 40, 1), "[OI]6300": (6302.05, 40, 2),
         "[OI]6363": (6365.54, 40, 2), "CaII8498": (8500.35, 40, 1), "CaII8542": (8544.44, 40, 1), "CaII8662": (8664.52, 40, 1),
         "MgII4481": (air2vac(4481.13), 40, 2), "FeII5169": (air2vac(5169.03), 40, 1)}


def load_steps():
    return json.load(open(STEPS)) if os.path.exists(STEPS) else {}


def save_step(name, res):
    S = load_steps(); S[name] = res; os.makedirs(CACHE, exist_ok=True)
    json.dump(S, open(STEPS, "w"), indent=1, default=lambda x: None if x is None else float(x))


def download(url, name, minsize):
    path = os.path.join(CACHE, name); os.makedirs(CACHE, exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) > minsize:
        return path
    for k in range(3):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                if r.status_code == 200:
                    with open(path + ".part", "wb") as fh:
                        for chunk in r.iter_content(1 << 20):
                            fh.write(chunk)
                    if os.path.getsize(path + ".part") > minsize and open(path + ".part", "rb").read(6) == b"SIMPLE":
                        os.replace(path + ".part", path); return path
        except requests.RequestException:
            pass
        time.sleep(5)
    raise IOError(f"download failed: {url}")


def spec_full(field, mjd, catid):
    f = f"{field:06d}"; name = f"spec-{f}-{mjd}-{catid}.fits"
    return download(f"{SAS}/spectra/daily/full/{f[:3]}XXX/{f}/{mjd}/{name}", name, 100000)


def spcframe(field, expo):
    f = f"{field:06d}"; name = f"spCFrame-r2-{expo:08d}.fits"
    return download(f"{SAS}/fields/{f[:3]}XXX/{f}/{name}", name, 50_000_000)


def read_spec(path, field):
    """Exposures (MJD_EXP_*) and the visit coadd (COADD) of one spec-full file."""
    out = []
    with fits.open(path) as h:
        mjd = int(h[0].header["MJD"])
        for k in range(1, len(h)):
            nm = h[k].name
            if not (nm.startswith("MJD_EXP") or nm == "COADD"):
                continue
            d = h[k].data; hd = h[k].header; am = np.asarray(d["AND_MASK"]).astype(np.int64); iv = np.asarray(d["IVAR"], float); ex = nm != "COADD"
            tai = hd.get("TAI-BEG") if ex else None
            out.append(dict(kind="exposure" if ex else "visit coadd", mjd=mjd, field=field, exposure=int(hd["EXPOSURE"]) if ex else None,
                            date_utc=Time(tai / 86400, format="mjd", scale="tai").utc.isot[:19] if tai else None, exptime_s=hd.get("EXPTIME") if ex else None,
                            airmass=hd.get("AIRMASS") if ex else None, seeing_arcsec=hd.get("SEEING50") if ex else None, helio_rv_kms=hd.get("HELIO_RV") if ex else None,
                            w=10 ** np.asarray(d["LOGLAM"], float), f=np.asarray(d["FLUX"], float), iv=np.where(am == 0, iv, 0.0), am=am, iv0=iv,
                            sky=np.asarray(d["SKY"], float), wdisp=np.asarray(d["WDISP"], float)))
    return out


def stack(ss):
    """Inverse-variance stack on the wavelength grid of the first spectrum (flux, sky model and WDISP)."""
    w = ss[0]["w"]; A = lambda s, k: s[k] if len(s["w"]) == len(w) and np.allclose(s["w"], w) else np.interp(w, s["w"], s[k], left=0, right=0)
    IV = np.array([A(s, "iv") for s in ss]); W = IV.sum(0); nz = np.where(W > 0, W, 1)
    f = np.where(W > 0, (np.array([A(s, "f") for s in ss]) * IV).sum(0) / nz, 0); sky = np.where(W > 0, (np.array([A(s, "sky") for s in ss]) * IV).sum(0) / nz, 0)
    return dict(w=w, f=f, iv=W, iv0=W, am=np.zeros(len(w), np.int64), sky=sky, wdisp=np.mean([A(s, "wdisp") for s in ss], axis=0))


def window_ew(w, f, iv, line):
    lam, win, cont = WINDOW[line]; x = (w - lam) / 100
    iv = np.where(np.isfinite(f) & np.isfinite(iv) & (iv > 0), iv, 0.0); sky = np.zeros(len(w), bool); ic = np.zeros(len(w), bool)
    for lo, hi in FIXED_SKY:
        sky |= (w >= lo) & (w <= hi)
    for lo, hi in cont:
        ic |= (w >= lo) & (w <= hi)
    use = ic & (iv > 0) & ~sky
    for _ in range(3):
        p, cov = np.polyfit(x[use], f[use], 1, w=np.sqrt(iv[use]), cov="unscaled")
        new = use & (np.abs((f - np.polyval(p, x)) * np.sqrt(iv)) < 3)
        if new.sum() == use.sum() or new.sum() < 8:
            break
        use = new
    c = np.polyval(p, x); ec0 = float(np.sqrt(max(cov[-1, -1], 0)))
    iw = (w >= win[0]) & (w <= win[1]); good = iw & (iv > 0) & ~sky & (c > 0)
    if good.sum() < max(2, 0.5 * iw.sum()):
        return np.nan, np.nan
    dl = np.abs(np.gradient(w))[good]; cc = c[good]
    ew = float(np.sum((1 - f[good] / cc) * dl)); e = float(np.hypot(np.sqrt(np.sum((dl / cc) ** 2 / iv[good])), abs(ec0 / np.median(cc)) * dl.sum()))
    return ew, e


def cfit(fn, x, y, e, p0, bounds):
    p, c = curve_fit(fn, x, y, p0=p0, sigma=e, absolute_sigma=True, bounds=bounds, maxfev=40000)
    return p, c, float(np.sum(((y - fn(x, *p)) / e) ** 2))


def single(w, f, iv, line, wdisp=None, sky=None):
    """Linear continuum + one Gaussian over the FIT range, with model comparisons."""
    lam0 = REF[line]; lo, hi = FIT[line]; m = (w > lo) & (w < hi) & (iv > 0) & np.isfinite(f); x, y, e = w[m], f[m], 1 / np.sqrt(iv[m]); n = len(x)
    g = lambda xx, c0, c1, a, mu, s: c0 + c1 * (xx - lam0) / 100 + a * np.exp(-0.5 * ((xx - mu) / s) ** 2)
    p, c, chi1 = cfit(g, x, y, e, [np.median(y), 0, 0.3 * np.median(y), lam0, 6.0], ([-np.inf] * 3 + [lam0 - 20, 1.0], [np.inf] * 3 + [lam0 + 20, 20.0]))
    r = max(1.0, chi1 / (n - 5)); pe = np.sqrt(np.diag(c) * r); fl = p[2] * p[4] * np.sqrt(2 * np.pi)
    J = np.array([0, 0, p[4] * np.sqrt(2 * np.pi), 0, p[2] * np.sqrt(2 * np.pi)]); efl = float(np.sqrt(J @ c @ J * r))
    o = {"ew_gauss_A": -fl / p[0], "e_ew_gauss_A": efl / p[0], "flux_single": fl, "e_flux_single": efl, "centre_A": p[3], "e_centre_A": pe[3],
         "v_single_kms": (p[3] / lam0 - 1) * C, "e_v_single_kms": pe[3] / lam0 * C, "sigma_A": p[4], "e_sigma_A": pe[4]}
    if wdisp is not None:
        si = float(np.median(wdisp[m]) * 1e-4 * np.log(10) * lam0); sd = np.sqrt(max(p[4] ** 2 - si ** 2, 0))
        o.update({"sigma_inst_kms": si / lam0 * C, "sigma_deconv_kms": sd / lam0 * C, "e_sigma_kms": pe[4] / lam0 * C})
    A = np.vstack([np.ones(n), (x - lam0) / 100]).T / e[:, None]; q = np.linalg.lstsq(A, y / e, rcond=None)[0]
    o["dbic_line_vs_continuum"] = float(np.sum((y / e - A @ q) ** 2) - chi1 - 3 * np.log(n))
    if sky is not None:
        B = np.vstack([np.ones(n), (x - lam0) / 100, sky[m]]).T / e[:, None]; q = np.linalg.lstsq(B, y / e, rcond=None)[0]
        o["sky_scale_single"] = float(q[2]); o["dbic_line_vs_sky"] = float(np.sum((y / e - B @ q) ** 2) - chi1 - 2 * np.log(n))
    g2 = lambda xx, c0, c1, a1, a2, mu, d, s: c0 + c1 * (xx - lam0) / 100 + a1 * np.exp(-0.5 * ((xx - mu + d / 2) / s) ** 2) + a2 * np.exp(-0.5 * ((xx - mu - d / 2) / s) ** 2)
    best = None
    for d0 in (4.0, 8.0, 12.0, 16.0):
        try:
            p2, _, ch2 = cfit(g2, x, y, e, [p[0], p[1], p[2] / 2, p[2] / 2, p[3], d0, max(1.5, p[4] / 2)],
                              ([-np.inf, -np.inf, 0, 0, lam0 - 20, 1.0, 1.0], [np.inf] * 4 + [lam0 + 20, 40.0, 20.0]))
            best = (ch2, p2) if best is None or ch2 < best[0] else best
        except (RuntimeError, ValueError):
            pass
    if best is not None:
        o["dbic_two_vs_one"] = float(chi1 - best[0] - 2 * np.log(n))
    wl, wh = WINDOW[line][1]; iw = (w >= wl) & (w <= wh) & (iv > 0)
    cm = (((w > lo) & (w < wl - 6)) | ((w > wh + 6) & (w < hi))) & (iv > 0); pc = np.polyfit(w[cm], f[cm], 1, w=np.sqrt(iv[cm])); ex = (f - np.polyval(pc, w))[iw]
    o["max_pixel_share"] = float(ex.max() / ex.sum()) if ex.sum() > 0 else np.nan
    return o, p


def joint(w, f, iv, sky=None):
    """Both O I lines: per +-45 A window a linear continuum, beta x sky model (common), two Gaussians with common velocity and width."""
    ok = (iv > 0) & np.isfinite(f) & (np.isfinite(sky) if sky is not None else True); L = [REF[k] for k in LINES]
    ms = [(np.abs(w - l) < 45) & ok for l in L]; x = [w[m] for m in ms]; y = [f[m] for m in ms]; e = [1 / np.sqrt(iv[m]) for m in ms]
    s = [sky[m] if sky is not None else np.zeros(m.sum()) for m in ms]

    def resid(p, lines=True):
        out = []
        for j, l in enumerate(L):
            mm = p[1 + 2 * j] + p[2 + 2 * j] * (x[j] - l) + p[0] * s[j]
            if lines:
                sg = l * p[8] / C; mm = mm + p[5 + j] / (np.sqrt(2 * np.pi) * sg) * np.exp(-0.5 * ((x[j] - l * (1 + p[7] / C)) / sg) ** 2)
            out.append((y[j] - mm) / e[j])
        return np.concatenate(out)
    lo = np.array([-1, -1e4, -1, -1e4, -1, -300, -300, -1500, 30.]); hi = np.array([1, 1e4, 1, 1e4, 1, 300, 300, 1500, 1000.])
    p0 = np.clip([0, np.median(y[0]), 0, np.median(y[1]), 0, 30, 30, 50, 250], lo + 1e-3, hi - 1e-3)
    rg = least_squares(resid, p0, bounds=(lo, hi)); rn = least_squares(lambda p: resid(p, False), p0, bounds=(lo, hi))
    pe = np.sqrt(np.diag(np.linalg.pinv(rg.jac.T @ rg.jac))); d = float((rn.fun ** 2).sum() - (rg.fun ** 2).sum())
    return {"joint_flux_o7774": rg.x[5], "joint_e_flux_o7774": pe[5], "joint_flux_o8446": rg.x[6], "joint_e_flux_o8446": pe[6], "joint_v_kms": rg.x[7],
            "joint_e_v_kms": pe[7], "joint_sigma_v_kms": rg.x[8], "joint_e_sigma_v_kms": pe[8], **({"joint_sky_beta": rg.x[0], "joint_e_sky_beta": pe[0]} if sky is not None else {}),
            "joint_dchi2": d, "joint_p": float(chi2dist.sf(d, 4)), "joint_n_pixels": len(rg.fun)}


def fixed_line(w, f, iv, lam, v, sig_kms, half, order):
    mu = lam * (1 + v / C); s = max(sig_kms / C * lam, 1.4); m = (np.abs(w - mu) < half) & (iv > 0)
    x, y, e = w[m], f[m], 1 / np.sqrt(iv[m])
    A = np.vstack([((x - mu) / half) ** k for k in range(order + 1)] + [np.exp(-0.5 * ((x - mu) / s) ** 2)]).T
    cov = np.linalg.inv((A / e[:, None] ** 2).T @ A); p = cov @ ((A / e[:, None] ** 2).T @ y); r = max(1, np.sum(((y - A @ p) / e) ** 2) / (len(x) - order - 2))
    return float(-p[-1] * s * np.sqrt(2 * np.pi) / p[0]), float(np.sqrt(cov[-1, -1] * r) * s * np.sqrt(2 * np.pi) / p[0])


def measure(s):
    o = {}; w, f, iv = s["w"], s["f"], s["iv"]; m = (np.abs(w - 7776.2) < 40) & (iv > 0); o["snr_7774"] = float(np.median(f[m] * np.sqrt(iv[m])))
    for line in LINES:
        ew, e = window_ew(w, f, iv, line); o.update({f"{line}_ew_window_A": ew, f"{line}_e_ew_window_A": e, f"{line}_s_window": -ew / e})
        g, _ = single(w, f, iv, line, s["wdisp"], s.get("sky")); o.update({f"{line}_{k}": v for k, v in g.items()})
        o[f"{line}_n_masked"] = int(((w >= WINDOW[line][1][0]) & (w <= WINDOW[line][1][1]) & ((s["am"] != 0) | (s["iv0"] <= 0))).sum())
    o.update(joint(w, f, iv, s.get("sky")))
    a, ea, b, eb = o["joint_flux_o8446"], o["joint_e_flux_o8446"], o["joint_flux_o7774"], o["joint_e_flux_o7774"]
    o["joint_ratio_8446_7774"] = a / b; o["joint_e_ratio_8446_7774"] = a / b * np.hypot(ea / a, eb / b)
    return o


def spectra():
    """(gaia, identifier, dataset, spectrum) for every spectrum of this topic."""
    out = []
    for g, specs in SPEC.items():
        ex = []
        for field, mjd, cat in specs:
            for s in read_spec(spec_full(field, mjd, cat), field):
                ident = f"exposure {s['exposure']}" if s["kind"] == "exposure" else f"visit {mjd}"
                out.append((g, ident, f"SDSS-V DR20 BOSS {s['kind']} (spec-full)", s)); ex += [s] if s["kind"] == "exposure" else []
        if g == GAIA:
            s = stack(ex); s.update(kind="stack", field=specs[0][0], mjd=None)
            out.append((g, f"exposures {', '.join(str(x['exposure']) for x in ex)}", "SDSS-V DR20 BOSS exposures stacked (spec-full)", s))
    return out


def wmean(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float); w = 1 / e ** 2; mu = float(np.sum(v * w) / w.sum())
    return mu, float(1 / np.sqrt(w.sum())), float(np.sum(((v - mu) / e) ** 2))


def lines():
    res = {"spectra": {}, "objects": {GAIA: {}, J0914: {}}}
    for g, ident, dataset, s in spectra():
        if g == J0914 and s["kind"] != "visit coadd":
            continue
        o = measure(s); meta = {k: s.get(k) for k in ("field", "mjd", "exposure", "date_utc", "exptime_s", "airmass", "seeing_arcsec", "helio_rv_kms")}
        if g == GAIA:
            res["spectra"][ident] = dict(gaia_dr3=g, name=NAME, dataset=dataset, identifier=ident, fibre=FIBRE, **meta, **o)
        if g == J0914 or s["kind"] == "stack":
            ob = res["objects"][g]; ob["lines_spectrum"] = f"{dataset}: {ident}"
            ob.update({k: o[k] for k in o if k.startswith(("joint_", "snr")) or k.endswith(("sigma_deconv_kms", "e_sigma_kms", "sigma_inst_kms", "ew_gauss_A", "ew_window_A", "_s_window"))})
            if g == GAIA:
                for name, (lam, half, order) in OTHER.items():
                    ew, e = fixed_line(s["w"], s["f"], s["iv"], lam, o["joint_v_kms"], o["joint_sigma_v_kms"], half, order)
                    ob.update({f"ew_{name}_A": ew, f"e_ew_{name}_A": e, f"ul3_{name}_A": max(0.0, -ew) + 3 * e})
                k = ["CaII8498", "CaII8542", "CaII8662"]; ew = sum(ob[f"ew_{x}_A"] for x in k); e = float(np.sqrt(sum(ob[f"e_ew_{x}_A"] ** 2 for x in k)))
                ob.update({"ew_CaII_triplet_A": ew, "e_ew_CaII_triplet_A": e, "ul3_CaII_triplet_A": max(0.0, -ew) + 3 * e})
        print(f"{g} {ident}: window EW 7774 {o['o7774_ew_window_A']:+.2f} +- {o['o7774_e_ew_window_A']:.2f}, 8446 {o['o8446_ew_window_A']:+.2f} +- {o['o8446_e_ew_window_A']:.2f}; "
              f"single centres {o['o7774_centre_A']:.2f} +- {o['o7774_e_centre_A']:.2f}, {o['o8446_centre_A']:.2f} +- {o['o8446_e_centre_A']:.2f}; joint F {o['joint_flux_o7774']:.1f} "
              f"+- {o['joint_e_flux_o7774']:.1f}, {o['joint_flux_o8446']:.1f} +- {o['joint_e_flux_o8446']:.1f}, v {o['joint_v_kms']:+.0f} +- {o['joint_e_v_kms']:.0f}, sigma_v "
              f"{o['joint_sigma_v_kms']:.0f} +- {o['joint_e_sigma_v_kms']:.0f}, dchi2 {o['joint_dchi2']:.1f}", flush=True)
    e = [r for r in res["spectra"].values() if r["dataset"].endswith("exposure (spec-full)")]; ob = res["objects"][GAIA]
    for line in LINES:
        ob[f"{line}_exposures_ew_window_chi2"] = wmean([r[f"{line}_ew_window_A"] for r in e], [r[f"{line}_e_ew_window_A"] for r in e])[2]
        ob[f"{line}_exposures_centre_A"], ob[f"{line}_exposures_e_centre_A"], ob[f"{line}_exposures_centre_chi2"] = wmean([r[f"{line}_centre_A"] for r in e], [r[f"{line}_e_centre_A"] for r in e])
    a, ea, b, eb = ob["o8446_exposures_centre_A"], ob["o8446_exposures_e_centre_A"], ob["o7774_exposures_centre_A"], ob["o7774_exposures_e_centre_A"]
    ob.update({"centre_ratio": a / b, "e_centre_ratio": a / b * np.hypot(ea / a, eb / b), "centre_ratio_oi": REF["o8446"] / REF["o7774"],
               "centre_v_o7774_kms": (b / REF["o7774"] - 1) * C, "centre_v_o8446_kms": (a / REF["o8446"] - 1) * C})
    print({k: round(v, 5) if isinstance(v, float) else v for k, v in ob.items() if "exposures" in k or "centre" in k})
    save_step("lines", res)


def wz(wl, f, iv, c, half=15.0, inner=25.0, outer=70.0):
    """Window excess over +-half A above a clipped linear continuum from inner-outer A on both sides: (excess, error)."""
    m = (np.abs(wl - c) > inner) & (np.abs(wl - c) < outer) & (iv > 0)
    if m.sum() < 10:
        return np.nan, np.nan
    x, y, w = wl[m] - c, f[m], iv[m]; keep = np.ones(len(x), bool)
    for _ in range(3):
        p = np.polyfit(x[keep], y[keep], 1, w=np.sqrt(w[keep])); r = (y - np.polyval(p, x)) * np.sqrt(w)
        nk = np.abs(r) < 3 * max(1.4826 * np.median(np.abs(r[keep])), 1)
        if (nk == keep).all():
            break
        keep = nk
    mw = (np.abs(wl - c) <= half) & (iv > 0)
    if mw.sum() < 5:
        return np.nan, np.nan
    dl = np.abs(np.gradient(wl))[mw]
    return float(((f[mw] - np.polyval(p, wl[mw] - c)) * dl).sum()), float(np.sqrt((dl ** 2 / iv[mw]).sum()))


def nightglow(w, sky):
    o = {}; m = (w > 7766) & (w < 7786)
    g = lambda x, a, mu, s, c: a / (np.sqrt(2 * np.pi) * s) * np.exp(-0.5 * ((x - mu) / s) ** 2) + c
    p, _ = curve_fit(g, w[m], sky[m], p0=[50, 7776.0, 1.6, np.percentile(sky[m], 10)], maxfev=20000)
    o.update(nightglow_flux_o7774=p[0], nightglow_v_o7774_kms=(p[1] / REF["o7774"] - 1) * C, nightglow_sigma_o7774_A=abs(p[2]))
    m = (w > 8438) & (w < 8452); g2 = lambda x, a, mu, s, c, d: g(x, a, mu, s, c) + d * (x - 8449)
    p, _ = curve_fit(g2, w[m], sky[m], p0=[15, 8449, 1.5, 2, 0], bounds=([-100, 8447.5, 1.0, -50, -5], [300, 8450.5, 2.5, 50, 5]))
    o.update(nightglow_flux_o8446=p[0], nightglow_v_o8446_kms=(p[1] / REF["o8446"] - 1) * C, nightglow_sigma_o8446_A=p[2])
    return o


NULL_LAMS = [l for l in np.arange(7000, 7960, 10.0) if not (7560 < l < 7740) and abs(l - 7777) > 40]


def frames():
    res = {"spectra": {}, "objects": {GAIA: {}}}; R = REF["o8446"] / REF["o7774"]
    for field, expo in [(109553, e) for e in S1_FRAMES] + list(CMP_FRAMES):
        own = expo in S1_FRAMES
        with fits.open(spcframe(field, expo), memmap=False) as h:
            hd = h[0].header; F = h[0].data.astype(float); IV = h[1].data.astype(float); MK = h[2].data; WL = h[3].data.astype(float); SK = h[6].data.astype(float); pm = h[5].data
        WL = 10 ** WL if np.nanmedian(WL) < 10 else WL; IV = IV * np.isfinite(F); F = np.nan_to_num(F)
        ot = np.array([str(x).strip().lower() for x in pm["OBJTYPE"]]); fib = np.asarray(pm["FIBERID"]); i = int(np.where(fib == FIBRE)[0][0])
        o = dict(gaia_dr3=GAIA, name=NAME, dataset="SDSS-V DR20 BOSS spCFrame-r2", identifier=f"exposure {expo}" if own else f"exposure {expo} (fibre {FIBRE} elsewhere)",
                 field=field, mjd=int(hd["MJD"]), exposure=expo, design_id=hd.get("DESIGNID"), helio_rv_kms=hd.get("HELIO_RV"), fibre=FIBRE,
                 fibre_objtype=ot[i], fibre_catalogid=str(pm["CATALOGID"][i]).strip())
        if not own:
            o["date_utc"] = Time(hd["TAI-BEG"] / 86400, format="mjd", scale="tai").utc.isot[:19]; o["airmass"] = hd.get("AIRMASS"); o["seeing_arcsec"] = hd.get("SEEING50"); o["exptime_s"] = hd.get("EXPTIME")
        o.update({f"frame_{k}": v for k, v in joint(WL[i], F[i], IV[i], SK[i]).items()})
        o.update(nightglow(WL[i], SK[i]))
        z = {}
        for line, c in ZWIN.items():
            m = np.abs(WL[i] - c) < 20; o[f"{line}_mask_pixels_fibre"] = int(((MK[i][m] & ~1) != 0).sum())
            z[line] = np.array([wz(WL[j], F[j], IV[j], c) for j in range(len(fib))]); zz = z[line][:, 0] / z[line][:, 1]
            o.update({f"{line}_window_flux_frame": z[line][i, 0], f"{line}_e_window_flux_frame": z[line][i, 1], f"{line}_z_frame": zz[i]})
            sky = (ot == "sky") & np.isfinite(zz); o.update({f"{line}_sky_fibres_n": int(sky.sum()), f"{line}_sky_fibres_z_median": float(np.median(zz[sky])),
                                                            f"{line}_sky_fibres_z_max": float(zz[sky].max())})
            for nb in (FIBRE - 1, FIBRE + 1):
                o[f"{line}_z_fibre_{nb}"] = float(zz[fib == nb][0])
            z[line] = zz
        both = (z["o7774"] > 2.5) & (z["o8446"] > 2.5)
        o["sky_fibres_both_z_gt_2p5"] = int((both & (ot == "sky")).sum()); o["other_fibres_both_z_gt_2p5"] = int((both & (fib != FIBRE)).sum())
        o["fibre_446_objtype"], o["fibre_448_objtype"] = ot[fib == FIBRE - 1][0], ot[fib == FIBRE + 1][0]
        if own:  # per-exposure chance rate from all other science and sky fibres of the frame
            A, B = [], []
            for j in np.where(np.isin(ot, ["sky", "science"]) & (fib != FIBRE))[0]:
                sl = (WL[j] > 6900) & (WL[j] < 8800)
                for l in NULL_LAMS:
                    a = wz(WL[j][sl], F[j][sl], IV[j][sl], l); b = wz(WL[j][sl], F[j][sl], IV[j][sl], l * R); A.append(a[0] / a[1]); B.append(b[0] / b[1])
            A, B = np.array(A), np.array(B); ok = np.isfinite(A) & np.isfinite(B); A, B = A[ok], B[ok]; hit = (A >= o["o7774_z_frame"]) & (B >= o["o8446_z_frame"])
            o.update(null_pairs=len(A), null_pairs_both=int(hit.sum()), null_rate=float(hit.mean()))
        res["spectra"][o["identifier"]] = o
        print(f"{o['identifier']} ({o['fibre_objtype']}): joint F {o['frame_joint_flux_o7774']:.1f} +- {o['frame_joint_e_flux_o7774']:.1f}, {o['frame_joint_flux_o8446']:.1f} +- "
              f"{o['frame_joint_e_flux_o8446']:.1f}, v {o['frame_joint_v_kms']:+.0f} +- {o['frame_joint_e_v_kms']:.0f}, sigma_v {o['frame_joint_sigma_v_kms']:.0f}, beta "
              f"{o['frame_joint_sky_beta']:+.3f}, dchi2 {o['frame_joint_dchi2']:.1f}; nightglow {o['nightglow_flux_o7774']:.1f} at {o['nightglow_v_o7774_kms']:+.0f} km/s, "
              f"{o['nightglow_flux_o8446']:.1f}; z {o['o7774_z_frame']:+.2f} / {o['o8446_z_frame']:+.2f}; sky fibres both > 2.5: {o['sky_fibres_both_z_gt_2p5']}"
              + (f"; chance rate {o['null_rate']:.1e} ({o['null_pairs_both']} of {o['null_pairs']})" if own else ""), flush=True)
    S = [res["spectra"][f"exposure {e}"] for e in S1_FRAMES]; ob = res["objects"][GAIA]
    for k in ("flux_o7774", "flux_o8446", "v_kms", "sigma_v_kms"):
        ob[f"frames_{k}"], ob[f"frames_e_{k}"], ob[f"frames_chi2_{k}"] = wmean([r[f"frame_joint_{k}"] for r in S], [r[f"frame_joint_e_{k}"] for r in S])
    a, ea, b, eb = ob["frames_flux_o8446"], ob["frames_e_flux_o8446"], ob["frames_flux_o7774"], ob["frames_e_flux_o7774"]
    ob.update(frames_ratio_8446_7774=a / b, frames_e_ratio_8446_7774=a / b * np.hypot(ea / a, eb / b))
    ob["frames_beta_max_abs"] = float(max(abs(r["frame_joint_sky_beta"]) for r in S))
    n1 = res["spectra"][f"exposure {S1_FRAMES[0]}"]
    ob.update(night1_line_over_nightglow_o7774=n1["frame_joint_flux_o7774"] / n1["nightglow_flux_o7774"], night1_line_over_nightglow_o8446=n1["frame_joint_flux_o8446"] / n1["nightglow_flux_o8446"])
    v = [r["nightglow_v_o7774_kms"] for r in S]; ob["nightglow_v_o7774_kms"] = float(np.mean(v)); ob["nightglow_e_v_o7774_kms"] = float(np.std(v, ddof=1) / np.sqrt(len(v)))
    cmp_sky = [r for r in res["spectra"].values() if "elsewhere" in r["identifier"] and r["fibre_objtype"] == "sky"]
    for line in LINES:
        m, e, _ = wmean([r[f"frame_joint_flux_{line}"] for r in cmp_sky], [r[f"frame_joint_e_flux_{line}"] for r in cmp_sky])
        ob[f"fibre_blank_sky_flux_{line}"], ob[f"fibre_blank_sky_e_flux_{line}"] = m, e
        ob[f"frames_flux_minus_blank_{line}"] = ob[f"frames_flux_{line}"] - m; ob[f"frames_e_flux_minus_blank_{line}"] = float(np.hypot(ob[f"frames_e_flux_{line}"], e))
    ob["fibre_blank_sky_frames"] = ", ".join(str(r["exposure"]) for r in cmp_sky)
    rates = {e: res["spectra"][f"exposure {e}"]["null_rate"] for e in S1_FRAMES}; night2 = min(rates[e] for e in S1_FRAMES[1:])
    ob.update(null_rate_night1=rates[S1_FRAMES[0]], null_rate_night2_min=night2, null_expected_sdssv=rates[S1_FRAMES[0]] * night2 * N_SDSSV)
    print({k: round(v, 4) if isinstance(v, float) else v for k, v in ob.items()})
    save_step("frames", res)


def tap(url, q, post=False, tries=4):
    kw = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)
    for k in range(tries):  # the ESA archive answers 503 when overloaded
        r = requests.post(url, data=kw, timeout=600) if post else requests.get(url, params=kw, timeout=600)
        if r.status_code not in (502, 503, 504) or k == tries - 1:
            break
        time.sleep(30)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


def sed():
    GA = "https://gea.esac.esa.int/tap-server/tap/sync"
    gaia = tap(GA, "SELECT source_id, ra, dec, pmra, pmdec, parallax, parallax_error, ruwe, phot_g_mean_mag, bp_rp, phot_variable_flag, has_epoch_photometry "
                   f"FROM gaiadr3.gaia_source WHERE source_id IN ({GAIA}, {J0914})", post=True).set_index("source_id")
    me = gaia.loc[int(GAIA)]
    nb = tap(GA, "SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, phot_variable_flag FROM gaiadr3.gaia_source WHERE "
                 f"1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {me.ra}, {me.dec}, {30 / 3600}))", post=True)
    nb["sep"] = np.hypot((nb.ra - me.ra) * np.cos(np.radians(me.dec)), nb.dec - me.dec) * 3600; nb = nb[nb.source_id != int(GAIA)].sort_values("sep")
    sw = cas("SELECT sdss_id, gaia_dr3_source_id, classification, p_da, p_cv, p_da_ms, teff, e_teff, logg, e_logg FROM snow_white_boss_star "
             f"WHERE gaia_dr3_source_id IN ({GAIA}, {J0914})").set_index("gaia_dr3_source_id")
    res = {"objects": {}}
    for g in (GAIA, J0914):
        x = gaia.loc[int(g)]; s = sw.loc[int(g)] if int(g) in sw.index else pd.Series(dtype=object, index=sw.columns)
        s = s.iloc[0] if isinstance(s, pd.DataFrame) else s
        res["objects"][g] = dict(ra_deg=x.ra, dec_deg=x.dec, G=x.phot_g_mean_mag, bp_rp=x.bp_rp, parallax_mas=x.parallax, e_parallax_mas=x.parallax_error,
                                 M_G=x.phot_g_mean_mag + 5 * np.log10(x.parallax / 100), pmra=x.pmra, pmdec=x.pmdec, ruwe=x.ruwe,
                                 gaia_phot_variable_flag=x.phot_variable_flag, gaia_epoch_photometry=str(x.has_epoch_photometry).lower() in ("true", "1"),
                                 sdss_id=str(s.sdss_id), snowwhite_class=s.classification, snowwhite_p_da=s.p_da, snowwhite_p_cv=s.p_cv, snowwhite_p_da_ms=s.p_da_ms,
                                 snowwhite_teff_K=s.teff, snowwhite_e_teff_K=s.e_teff, snowwhite_logg=s.logg, snowwhite_e_logg=s.e_logg,
                                 gf21_teff_K=GF21[g][0], gf21_logg=GF21[g][1], tic=TIC.get(g, ""))
    ob = res["objects"][GAIA]
    ob["gaia_neighbours_30arcsec"] = "; ".join(f"Gaia DR3 {int(r.source_id)} at {r.sep:.1f} arcsec, G {r.phot_g_mean_mag:.2f}, BP-RP {r.bp_rp:.2f}"
                                               + (f", {r.phot_variable_flag}" if r.phot_variable_flag == "VARIABLE" else "") for r in nb.itertuples())
    P = pd.read_csv(os.path.join(DATA, "oxygen_emission_photometry.csv")); p = os.path.join(CACHE, "Table_DA")
    if not os.path.exists(p):
        os.makedirs(CACHE, exist_ok=True); open(p, "w").write(requests.get("https://www.astro.umontreal.ca/~bergeron/CoolingModels/Tables/Table_DA", timeout=120).text)
    L = open(p).read().splitlines(); nh = len(L[1].split()) - 1  # 'log g' is two tokens
    A = np.array([[float(v) for v in l.split()] for l in L[2:] if len(l.split()) == nh])
    Tq = 10 ** np.arange(np.log10(10000), np.log10(45000), 0.004); Gq = np.arange(7.0, 9.0001, 0.02); gs = np.unique(A[:, 1]); cols = list(P.table_da_column) + [2]
    blk = np.zeros((len(Tq), len(gs), len(cols)))
    for j, lg in enumerate(gs):
        B = A[A[:, 1] == lg]; B = B[np.argsort(B[:, 0])]
        for k, c in enumerate(cols):
            blk[:, j, k] = np.interp(np.log10(Tq), np.log10(B[:, 0]), B[:, c])
    M = np.zeros((len(Tq), len(Gq), len(cols)))
    for i, lg in enumerate(Gq):
        j = int(np.clip(np.searchsorted(gs, lg) - 1, 0, len(gs) - 2)); a = (lg - gs[j]) / (gs[j + 1] - gs[j]); M[:, i] = (1 - a) * blk[:, j] + a * blk[:, j + 1]
    mass = M[:, :, -1]; M = M[:, :, :-1]; plx, eplx = me.parallax + PLX_ZP, me.parallax_error; Pg = np.arange(0.9, 2.6, 0.01); mu = 5 * np.log10(100 / Pg)
    obs = (P.mag - P.R * EBV).values; sig = np.hypot(P.e_mag_used, P.sys_mag).values

    def chi(bands, T=None, lg=None):
        k = [list(P.band).index(b) for b in bands]
        mod = M[:, :, k] if T is None else M[[int(np.argmin(np.abs(Tq - T)))]][:, [int(np.argmin(np.abs(Gq - lg)))]][:, :, k]
        model = mod[:, :, None, :] + mu[None, None, :, None]
        return np.sum(((obs[k] - model) / sig[k]) ** 2, axis=-1) + ((Pg - plx) / eplx)[None, None, :] ** 2, model, k
    for tag, bands in (("sed", list(P.band)), ("sed_uv_gaia", ["FUV", "NUV", "G", "BP", "RP"]), ("sed_gaia_des_vhs", ["G", "BP", "RP", "g", "r", "i", "z", "Y", "J"]),
                       ("sed_gaia", ["G", "BP", "RP"])):
        X, model, k = chi(bands); i, j, q = np.unravel_index(np.argmin(X), X.shape); c0 = X[i, j, q]; ok = X.min(axis=2) - c0 <= 1
        ob.update({f"{tag}_teff_K": Tq[i], f"{tag}_teff_lo_K": Tq[ok.any(1)].min(), f"{tag}_teff_hi_K": Tq[ok.any(1)].max(), f"{tag}_logg": Gq[j],
                   f"{tag}_logg_lo": Gq[ok.any(0)].min(), f"{tag}_logg_hi": Gq[ok.any(0)].max(), f"{tag}_mass": mass[i, j], f"{tag}_mass_lo": mass[ok].min(),
                   f"{tag}_mass_hi": mass[ok].max(), f"{tag}_parallax_mas": Pg[q], f"{tag}_chi2": c0, f"{tag}_dof": len(bands) - 2})
        if tag == "sed":
            ob.update({f"sed_resid_{b}": obs[kk] - model[i, j, q, n] for n, (b, kk) in enumerate(zip(bands, k))})
        print(f"{tag}: Teff {Tq[i]:.0f} [{Tq[ok.any(1)].min():.0f}-{Tq[ok.any(1)].max():.0f}] log g {Gq[j]:.2f} [{Gq[ok.any(0)].min():.2f}-{Gq[ok.any(0)].max():.2f}] "
              f"M {mass[i, j]:.3f} plx {Pg[q]:.2f} chi2 {c0:.1f} / {len(bands) - 2}")
    for tag, (T, lg) in (("snowwhite", (ob["snowwhite_teff_K"], ob["snowwhite_logg"])), ("gf21", GF21[GAIA])):
        X, model, k = chi(list(P.band), T, lg); q = int(np.argmin(X[0, 0]))
        ob.update({f"sed_at_{tag}_chi2": X[0, 0, q], f"sed_at_{tag}_parallax_mas": Pg[q], f"sed_at_{tag}_resid_FUV": obs[0] - model[0, 0, q, 0],
                   f"sed_at_{tag}_resid_NUV": obs[1] - model[0, 0, q, 1]})
        print(f"photometry at {tag} Teff {T:.0f} log g {lg:.2f}: chi2 {X[0, 0, q]:.1f}, parallax {Pg[q]:.2f} mas")
    save_step("sed", res)


def stats(m, e):
    m, e = np.asarray(m, float), np.asarray(e, float); w = 1 / e ** 2; mu = np.sum(w * m) / w.sum(); n = len(m); var = np.var(m, ddof=1)
    return dict(n=n, rms=float(np.sqrt(var)), chi2r=float(np.sum(((m - mu) / e) ** 2) / (n - 1)),
                A95=float(np.sqrt(2 * (max(0.0, var - np.mean(e ** 2)) + 1.645 * np.sqrt(2 / (n - 1)) * var))))


def variability():
    ra, dec = 357.53671148528, -60.31592222605; d = 2 / 3600; dr = d / np.cos(np.radians(dec)); res = {"objects": {GAIA: {}}}; ob = res["objects"][GAIA]
    box = lambda a, b: f"{a} BETWEEN {ra - dr} AND {ra + dr} AND {b} BETWEEN {dec - d} AND {dec + d}"
    DL = "https://datalab.noirlab.edu/tap/sync"; o = tap(DL, f"SELECT id, ra, dec FROM nsc_dr2.object WHERE {box('ra', 'dec')}")
    n = tap(DL, f"SELECT mjd, filter, mag_auto, magerr_auto, flags FROM nsc_dr2.meas WHERE objectid='{o.id.iloc[0]}'")
    n = n[(n["flags"] == 0) & (n.magerr_auto < 0.2)]; ob["nsc_dr2_object"] = str(o.id.iloc[0]); ob["nsc_dr2_mjd_span"] = f"{n.mjd.min():.0f}-{n.mjd.max():.0f}"
    for b in "grizY":
        x = n[n["filter"] == b]
        if len(x) > 2:
            ob.update({f"des_{b}_{k}": v for k, v in stats(x.mag_auto, x.magerr_auto).items()})
    SM = "https://api.skymapper.nci.org.au/public/tap/sync"; o = tap(SM, f"SELECT object_id, raj2000, dej2000 FROM dr4.master WHERE {box('raj2000', 'dej2000')}")
    s = tap(SM, f"SELECT image_id, filter, mag_psf, e_mag_psf, flags, nimaflags FROM dr4.photometry WHERE object_id={int(o.object_id.iloc[0])}")
    s = s[(s["flags"] == 0) & (s.nimaflags == 0) & (s.e_mag_psf < 0.3)]; ob["skymapper_dr4_object"] = str(int(o.object_id.iloc[0]))
    for b in "uvgri":
        x = s[s["filter"] == b]
        if len(x) > 2:
            ob.update({f"skymapper_{b}_{k}": v for k, v in stats(x.mag_psf, x.e_mag_psf).items()})
    P = pd.read_csv(os.path.join(DATA, "oxygen_emission_photometry.csv")).set_index("band"); fl = lambda b: 10 ** (-0.4 * (P.mag[b] - 23.9))
    TOT = {"o": (fl("r") + fl("i")) / 2, "c": (fl("g") + fl("r")) / 2}  # total flux in uJy from DES
    L = [l for l in open(os.path.join(DATA, f"atlas_forced_photometry_{GAIA}.txt")).read().splitlines() if l.strip()]
    t = pd.DataFrame([l.split() for l in L[1:]], columns=L[0].lstrip("#").split())
    for c in ("MJD", "uJy", "duJy", "chi/N", "err"):
        t[c] = t[c].astype(float)
    t = t[(t.duJy > 0) & (t["chi/N"] < 3) & (t.err == 0)].copy(); t["unit"] = t.Obs.str[:3]
    fr = np.linspace(1 / 20, 50, 100000); keep = np.ones(len(fr), bool); rng = np.random.default_rng(7)
    for h in range(1, 40):
        for a in (-2, -1, 0, 1, 2):
            keep &= np.abs(fr - (h / RRL_P + a * 1.0027)) > 0.01
    amp = lambda x, y, W, f0: float(np.hypot(*np.linalg.lstsq(np.vstack([np.ones(len(x)), np.sin(2 * np.pi * f0 * x), np.cos(2 * np.pi * f0 * x)]).T * np.sqrt(W)[:, None],
                                                            y * np.sqrt(W), rcond=None)[0][1:]))
    for b in ("o", "c"):
        x = t[t.F == b]
        for _ in range(3):
            med = np.median(x.uJy); mad = 1.4826 * np.median(np.abs(x.uJy - med)); x = x[np.abs(x.uJy - med) < 5 * max(mad, 1)]
        y = x.uJy.values.copy(); e = x.duJy.values; mj = x.MJD.values; W = 1 / e ** 2
        for u in x.unit.unique():
            m = (x.unit == u).values; y[m] -= np.average(y[m], weights=W[m])
        f0 = fr[int(np.argmax(LombScargle(mj, y, e).power(fr)))]
        X = np.vstack([np.ones(len(mj))] + [fn(2 * np.pi * h * mj / RRL_P) for h in range(1, 5) for fn in (np.sin, np.cos)]).T
        p = np.linalg.lstsq(X * np.sqrt(W)[:, None], y * np.sqrt(W), rcond=None)[0]; r = y - X @ p
        pw = LombScargle(mj, r, e).power(fr); pw[~keep] = 0; f1 = fr[int(np.argmax(pw))]; A1 = amp(mj, r, W, f1); mx = []
        for _ in range(200):
            rs = rng.permutation(r); ps = LombScargle(mj, rs, e).power(fr[::5]); mx.append(amp(mj, rs, W, fr[::5][int(np.argmax(ps))]))
        ob.update({f"atlas_{b}_n": len(mj), f"atlas_{b}_mjd_span": f"{mj.min():.0f}-{mj.max():.0f}", f"atlas_{b}_median_error_uJy": float(np.median(e)),
                   f"atlas_{b}_total_flux_uJy": TOT[b], f"atlas_{b}_top_period_d": 1 / f0, f"atlas_{b}_top_semi_amplitude_uJy": amp(mj, y, W, f0),
                   f"atlas_{b}_rrl_semi_amplitude_uJy": float(np.hypot(p[1], p[2])), f"atlas_{b}_chi2r_before_rrl": float(np.sum(W * y ** 2) / (len(y) - 1)),
                   f"atlas_{b}_chi2r_after_rrl": float(np.sum(W * r ** 2) / (len(y) - 9)), f"atlas_{b}_residual_top_period_d": 1 / f1,
                   f"atlas_{b}_residual_top_semi_amplitude_uJy": A1, f"atlas_{b}_shuffled_p99_semi_amplitude_uJy": float(np.quantile(mx, 0.99))})
        print(f"ATLAS {b}: n {len(mj)}; top P {1 / f0:.5f} d ({amp(mj, y, W, f0):.1f} uJy); after the RR Lyrae fit chi2_r {ob[f'atlas_{b}_chi2r_before_rrl']:.2f} -> "
              f"{ob[f'atlas_{b}_chi2r_after_rrl']:.2f}; residual top P {1 / f1:.4f} d, {A1:.1f} uJy; shuffled 99 % {np.quantile(mx, 0.99):.1f} uJy; total {TOT[b]:.0f} uJy", flush=True)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in ob.items() if k.startswith(("des_", "skymapper_"))})
    save_step("variability", res)


DIG = (("helio_rv", 2), ("_kms", 1), ("dbic", 1), ("dchi2", 1), ("percentile", 1), ("_uJy", 1), ("chi2", 2), ("centre_ratio", 5), ("ratio", 4), ("resid", 3), ("A95", 4), ("rms", 4), ("beta", 3),
       ("scale", 3), ("share", 3), ("_K", 0), ("_z_", 2), ("_s_", 2), ("snr", 2), ("flux", 2), ("centre", 3), ("_A", 3), ("period_d", 5), ("M_G", 2))


def rnd(col, v):
    if not isinstance(v, float) or not np.isfinite(v) or col in ("ra_deg", "dec_deg") or col.startswith(("null_", "joint_p", "frame_joint_p")):
        return v
    for key, nd in DIG:
        if key in col:
            return round(v, nd)
    return round(v, 4)


def write_tables():
    S = load_steps(); rows = {}
    for step in ("frames", "lines"):  # for the four exposures the spec-full (lines) values of shared keys are kept
        for ident, d in S.get(step, {}).get("spectra", {}).items():
            rows.setdefault(ident, {}).update({k: v for k, v in d.items() if v is not None or k not in rows.get(ident, {})})
    for ident in set(S.get("frames", {}).get("spectra", {})) & set(S.get("lines", {}).get("spectra", {})):
        rows[ident]["dataset"] = "SDSS-V DR20 BOSS exposure (spec-full and spCFrame-r2)"
    if rows:
        t = pd.DataFrame(list(rows.values())); own = {f"exposure {e}" for e in S1_FRAMES}
        t["sort"] = t.identifier.map(lambda s: 0 if s in own else 1 if s.startswith("visit ") else 2 if s.startswith("exposures") else 3)
        t = t.sort_values(["sort", "mjd", "exposure"]).drop(columns="sort"); t = t[[c for c in t.columns if t[c].notna().any()]]
        for c in t.columns:
            t[c] = [rnd(c, v) for v in t[c]]
        for c in ("field", "mjd", "exposure", "design_id", "null_pairs", "null_pairs_both"):
            if c in t.columns:
                t[c] = t[c].astype("Int64")
        first = ["gaia_dr3", "name", "dataset", "identifier", "field", "mjd", "exposure", "design_id", "date_utc", "exptime_s", "airmass", "seeing_arcsec",
                 "helio_rv_kms", "fibre", "fibre_objtype", "fibre_catalogid"]
        t = t[[c for c in first if c in t.columns] + [c for c in t.columns if c not in first]]
        t.to_csv(os.path.join(TAB, "oxygen_emission_exposures.csv"), index=False)
    obj = {GAIA: dict(gaia_dr3=GAIA, name=NAME, role="object"), J0914: dict(gaia_dr3=J0914, name=J0914_NAME, role="comparison (line widths and flux ratio)")}
    for step in ("sed", "lines", "frames", "variability"):
        for g, d in S.get(step, {}).get("objects", {}).items():
            obj[g].update(d)
    t = pd.DataFrame(list(obj.values()))
    for c in t.columns:
        t[c] = [rnd(c, v) for v in t[c]]
        if c.endswith(("_n", "_dof", "_pixels", "_K")) and pd.api.types.is_numeric_dtype(t[c]):
            t[c] = t[c].round().astype("Int64")
    t.to_csv(os.path.join(TAB, "oxygen_emission.csv"), index=False)


if __name__ == "__main__":
    steps = dict(lines=lines, frames=frames, sed=sed, variability=variability)
    for k in (sys.argv[1:] or list(steps)):
        steps[k](); write_tables(); print(k, "done", flush=True)
