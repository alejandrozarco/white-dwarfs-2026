# Hydrogen and helium line strengths that change between visits

A screen of the SDSS-V DR20 visit spectra for white dwarfs whose hydrogen and helium lines change in opposite directions between visits, the signature of a "double-faced" white dwarf (a surface with hydrogen-rich and helium-rich regions that rotate into view). One candidate is listed: GALEX J085307.3-002204 (Gaia DR3 3074374016095070464). Methods are below and in [METHODS.md](../METHODS.md#methods).

## GALEX J085307.3-002204 (candidate)

Gaia DR3 3074374016095070464, SDSS J085307.37-002203.9, WDJ085307.38-002203.89; G = 18.80, BP-RP = -0.27, M_G = 11.01. SDSS-V SnowWhite class DA (36.9 kK); Gentile Fusillo et al. (2021) H-atmosphere fit 21.6 kK, 0.75 Msun. No other Gaia DR3 source within 5 arcsec.

Equivalent widths in every spectrum (`tables/line_variability_3074374016095070464.csv`; Angstrom, absorption positive; SDSS-V errors calibrated as described below, BOSS and DESI errors formal):

| spectrum | MJD | S/N | ZTF phase | EW H-alpha | EW H-beta | EW H-gamma | EW He I 4471 | EW He I 4922 | EW He I 5876 |
|---|---|---|---|---|---|---|---|---|---|
| SDSS-V visit (field 15164) | 59200 | 17.8 | 0.27-0.60 | 3.57 ± 0.76 | 6.52 ± 0.29 | 5.26 ± 0.30 | 2.24 ± 0.25 | 0.89 ± 0.26 | 1.38 ± 0.36 |
| SDSS-V visit (field 15166) | 59250 | 17.2 | 0.47-0.83 | 5.72 ± 0.85 | 8.16 ± 0.28 | 7.36 ± 0.30 | 2.06 ± 0.26 | 0.59 ± 0.26 | 1.33 ± 0.37 |
| SDSS-V visit (field 15166) | 59251 | 14.2 | 0.70-0.88 | 8.65 ± 0.99 | 11.20 ± 0.33 | 10.14 ± 0.36 | 1.63 ± 0.35 | -0.15 ± 0.36 | 1.34 ± 0.45 |
| BOSS DR17 (SPARCL) | 55509.48 | - | 0.72 ± 0.09 | 6.76 ± 0.65 | 9.67 ± 0.29 | 8.29 ± 0.29 | 1.39 ± 0.28 | 0.08 ± 0.32 | 1.45 ± 0.38 |
| DESI DR1 (SPARCL) | 59589.31 | - | 0.72 | 8.13 ± 1.17 | 10.47 ± 0.72 | 9.06 ± 0.69 | 0.05 ± 0.69 | 0.34 ± 0.79 | 2.28 ± 0.86 |

- Across the three SDSS-V visits EW(H-beta) rises from 6.52 to 11.20 A while EW(He I 4471) falls from 2.24 to 1.63 A. The H-beta change between MJD 59250 and 59251 (consecutive nights, same field) is 7.0 sigma.
- Variability statistics of the three visits (`tables/line_variability.csv`): hydrogen index p_ZH = 2.9e-38, helium index p_ZHe = 0.037, anticorrelation index p_ZA = 4.8e-30; ratio R = EW(He 4471 + 5876)/EW(H-beta + H-alpha) from 0.15 to 0.36, p_R = 1.4e-3; correlation of the hydrogen and helium indices -0.997.
- Other DA white dwarfs on the same fields show mean index offsets below 0.3 sigma on all three nights (17 and 33 plate-mates); the continuum colour change of the MJD 59200 visit is shared by its plate-mates (median blue/red flux ratio 0.72). These two checks are from the screen and are not recomputed by the scripts.
- Existing classifications: SIMBAD WD* with spectral type DAB (from the SDSS DR14 white-dwarf catalogue, Kepler et al. 2019); MWDD DA; DESI DR1 class DA (Swan et al. 2026; three classifiers). No ADS full-text match for seven aliases, no VSX entry, no period in Oliveira da Rosa et al. (2024) (checked 2026-10-18).

**ZTF period.** The ZTF g + r light curve (720 points, 2018-2025) has its highest Lomb-Scargle peak at f = 4.076139 ± 0.000022 c/d, P = 5.88793 ± 0.00003 h, semi-amplitude 2.96 ± 0.32%, Baluev false-alarm probability 1.2e-12.
- The peak is the highest in ZTF g alone (false-alarm probability 1.1e-5). In ZTF r its false-alarm probability is 0.018 (the highest r peak is at 2.0734 c/d, 0.007).
- None of 200 shuffles of the fluxes within each filter, field, CCD and quadrant reaches the observed power (highest shuffled power 0.052, observed 0.108). The spectral window has power 0.0005 at the peak.
- Power at the peak in each season from 2018 to 2025: 0.13, 0.08, 0.21, 0.08, 0.06, 0.11, 0.04, 0.21.
- The eight other ZTF objects within 90 arcsec with at least 100 g-band points have power 0.003-0.032 at this frequency.
- TESS: no 2-minute data (T = 19.07).

**Phases.** On this ephemeris (phase 0 = flux maximum), EW(H-beta) increases and EW(He I 4471) decreases with the phase of the three SDSS-V visits. Each visit spans 0.2-0.4 cycles, and three points do not establish a modulation. The BOSS and DESI phases carry only the formal frequency error and assume no cycle-count error over 4,000 days.

**Decisive test.** Time-resolved spectroscopy covering at least one 5.9-h period (a 4-8 m telescope at G = 18.8), with H-beta and He I 4471 measured against phase.

<img src="../figures/line_variability/3074374016095070464.png" width="760">

## Screen

- **Sample.** All 26,757 SDSS-V DR20 objects with at least two BOSS visits in the mwmVisit files, 79,437 visits. Visits with in_stack = False (13,593) are excluded; for the others the XCSAO shift is removed. The Astra mwmStar coadds are not used.
- **Visit S/N cut.** Per telescope, log S/N = poly2(G) + log(exposure time); visits more than a factor 2 above or 3 below the prediction are dropped (1,335 of 65,844).
- **Equivalent widths.** H-alpha, H-beta, H-gamma, He I 4471, 4922, 5876 and 6678, as in `scripts/line_variability.py`.
- **Error calibration.** One factor per line and atmosphere group, from the robust width of single-visit residuals in 12 bins of formal error. The null samples are SnowWhite class DA (hydrogen group) and DB (helium group). The factors are 1.4-1.8 for DA and 1.0-2.4 for DB (about 300 DB objects). The indices are rescaled in 6 S/N bins (`data/line_variability_calibration.csv`, `data/line_variability_index_calibration.csv`).
- **Analysed.** 21,162 objects with at least two usable visits (6,514 with at least three): 18,092 in the hydrogen group, 300 helium, 116 mixed (DAB, DBA, DA/DB), 277 DC, 2,377 other.
- **Null rate** at p_ZA < 1e-4 (`tables/line_variability_screen.csv`): 172 of 18,092 hydrogen-group objects (0.95%), 2 of 300 helium (0.67%), 0 of 116 mixed and 5 of 277 DC with at least two visits; 1.7% (hydrogen) and 1.9% (helium) with at least three visits.
- **Injection-recovery** (`tables/line_variability_injection.csv`). A surface fraction a = a_max (1 + cos phi)/2 with a random phase per visit is replaced by the template of the other atmosphere (DA sdss_id 75787962, 18.6 kK; DB sdss_id 54394901) in 600 random DAs and all 179 DBs with at least three visits. Recovery at p_ZA < 1e-4 with at least three visits:

  | targets | a_max (median swing between visits) | S/N < 8 | S/N 8-15 | S/N > 15 |
  |---|---|---|---|---|
  | DA | 0.15 (0.10-0.12) | 0.10 | 0.37 | 0.77 |
  | DA | 0.3 (0.21-0.24) | 0.36 | 0.82 | 0.93 |
  | DA | 0.6 (0.41-0.49) | 0.72 | 0.96 | 0.98 |
  | DB | 0.3 (0.21-0.26) | 0.57 | 0.72 | 0.97 |
  | DB | 0.6 (0.42-0.52) | 1.00 | 0.93 | 1.00 |

  The same flag on the uninjected targets: 1-3%.
- **Positive controls.** None of the known double-faced or spectroscopically variable white dwarfs can be recovered from the SDSS-V visits. J0910+2105 and GD 323 have one visit each. GALEX J0718+3731 has three visits and no variability, but its 11-min rotation is averaged within each visit. Janus (ZTF J2033+3229), SDSS J0847+4842, LB 8915, Feige 7, KUV 02196+2816 and the candidates WDJ2131+0255 and WDJ0917+0010 (Moss et al. 2025; Munday et al. 2026) have no SDSS-V visit spectra.
- **Negative controls.**
  - Ten unresolved DA+DB pairs with at least two visits all have p_ZA ≥ 0.03 and a constant ratio, p_R ≥ 0.12 (`tables/line_variability_dadb_controls.csv`).
  - 89 pairs of visits of one object on the same night (`tables/line_variability_sameday_pairs.csv`): robust width of the change in Z_A 0.85, none beyond 3 sigma.
  - 33,137 pairs of consecutive visits on different nights: robust width 0.97; 1.3% beyond 3 sigma and 0.53% at p < 1e-4, against 0.27% and 0.01% for a normal distribution.
- **Vetting.** 44 flagged objects were inspected by eye: every flag of the anticorrelation, helium-appearance and mixed-group tests, the ranked species-appearance list, and the 12 strongest single-species hydrogen flags. Apart from GALEX J085307.3-002204 the flags were traced to:
  - one visit with lines displaced by about 70 A (failed removal of the XCSAO shift);
  - one featureless, weak-line or other-source visit;
  - main-sequence stars classified DA:;
  - diffuse nebular emission at low Galactic latitude;
  - dilution of all lines in one visit with a constant line ratio;
  - sky or cosmic-ray spikes, broken continua and low-S/N noise.

  140 further single-species hydrogen flags in the white-dwarf locus were not inspected.
- **Limits.** Visits last 15-120 min, so rotation periods below about 2 h are averaged out. Visits with in_stack = False (17% of all visits) and BOSS per-exposure spectra were not used.

The screen tables (`line_variability_injection.csv`, `line_variability_screen.csv`, `line_variability_dadb_controls.csv`, `line_variability_sameday_pairs.csv`) are exported from the full run and are not recomputed by the scripts; the per-object tables are (`scripts/line_variability.py`).
