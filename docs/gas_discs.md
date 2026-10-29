# Ca II triplet emission (gaseous debris discs)

Tables:
- `gas_disc_white_dwarfs.csv`: six white dwarfs with Ca II triplet emission (8500, 8544, 8665 Å) and no existing emission classification: four from the SDSS-V screen (two with double-peaked and two with single-peaked profiles) and two from the DESI DR1 screen (one tentative) in SIMBAD, MWDD, the DESI DR1 white-dwarf catalogues (Swan et al. 2026; Amorim et al. 2026), VizieR, the gaseous-disc lists of Saker et al. (2025) and Ma et al. (2025), or ADS full text;
- `gas_disc_epochs_<gaia_dr3>.csv`: the Ca II emission equivalent width, profile scale and Gaussian centroid and width in every available spectrum;
- `desi_gas_disc_screen.csv`: DESI DR1 screen output for the two DESI stars and the known Ca II emitters in DESI DR1;
- `gas_disc_screen.csv`: SDSS-V screen output for the four SDSS-V stars, the known gaseous-disc white dwarfs with SDSS-V spectra and the DESI EDR emission candidates of Ma et al. (2025) with SDSS-V spectra;
- `gas_disc_epoch_checks.csv` (`gas_disc_checks.py`): for every spectrum of the five non-tentative stars, the equivalent width with its formal error, an error from line-free windows, the equivalent width with continuum orders 1 and 3, the profile scale, the continuum at 8600 Å, the line flux and the median flux at 4500-5000 Å;
- `gas_disc_variability.csv`: $`\chi^2`$ of the equivalent widths of the dated spectra of each star, with formal errors and with formal and line-free errors combined;
- `gas_disc_line_profiles.csv`: peak velocities, central minimum, FWHM and 25% crossing velocities of the Ca II 8664.52 Å line;
- `gas_disc_infrared.csv`: Legacy Surveys DR10, CatWISE2020 and 2MASS fluxes against Montreal model photometry, and a two-blackbody fit for WD J1959+2208;
- `gas_disc_578709631539357440_amplitude_tests.csv` and `gas_disc_578709631539357440_xshooter_metals.csv`: checks of the BOSS 2010 and SDSS 2003 spectra and the X-shooter UVB metal lines of WD 0856+048.

Methods are in [METHODS.md](../METHODS.md#methods).

## WD 0856+048 (Gaia DR3 578709631539357440)
G = 18.34, parallax 3.82 ± 0.19 mas; catalogued DA (SIMBAD, MWDD, SnowWhite, DESI DR1). The Ca II triplet shows double-peaked emission in all seven SDSS-V visits (2021-03-19 to 2023-11-22; z per visit 9 to 22), in the DESI spectrum of 2022-03-07 and in both ESO X-shooter spectra of 2025-04-23. The emission equivalent width (three lines) is 17.6 ± 0.5 Å in the SDSS-V coadd (visits 10.5-21.4 Å), 18.7 ± 1.8 Å (DESI 2022), 31.2 ± 0.8 and 32.9 ± 0.4 Å (X-shooter 2025) and 3.1 ± 2.3 Å in the BOSS spectrum of 2010-12-05. The SDSS spectrum of 2003-01-31 gives 15.0 ± 6.8 Å. In X-shooter the profile peaks are at about −310 and +360 km/s, and the UVB arm shows Mg II 4481 and Ca II K absorption.

The emission strength varies between epochs. With errors from line-free windows added to the formal errors, the ten dated spectra with an equivalent width above three times its error give $`\chi^2`$ = 120.6 for 9 degrees of freedom; the X-shooter values exceed the SDSS-V coadd under every continuum choice tested (30.3-33.8 Å against 14.6-17.6 Å; `gas_disc_epoch_checks.csv`). The scale of the SDSS-V coadd profile (`amp_rel_ref`) is 0.016 ± 0.093 in the BOSS spectrum of 2010-12-05, 0.63 ± 0.26 in the SDSS spectrum of 2003-01-31, 0.78-1.18 in the SDSS-V visits, 0.91 ± 0.07 in DESI 2022 and 1.69 and 1.78 in X-shooter 2025. No emission is detected in the BOSS 2010 spectrum. Checks of that spectrum (`gas_disc_578709631539357440_amplitude_tests.csv`):
- the scale is −0.13 to 0.26 for 24 continuum choices (order 1-4, three fitting ranges, two line masks; the reference profile re-derived with each);
- applying the fit to its residuals circularly shifted by 40 to 379 pixels gives a mean of 0.06, a standard deviation of 0.08 and a maximum of 0.26 (340 shifts);
- the coadd profile added to the spectrum at scale 0.5 and 1.0 before normalisation is recovered at 0.52 and 1.02;
- its continuum at 8600 Å, 3.46 × 10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹, is within the range of the SDSS-V visits (2.49-3.86). The BOSS 2010 and SDSS-V spectra come from the BOSS spectrographs on the 2.5-m Sloan telescope.

The SDSS 2003 spectrum (S/N about 2.5 per pixel near 8600 Å) gives a scale of −0.51 to 0.78 over the same continuum choices, and its circular-shift scatter is 0.25: it does not distinguish between no emission and emission at the SDSS-V level. There are no spectra between 2010-12-05 and 2021-03-19.

Line profile (`gas_disc_line_profiles.csv`, Ca II 8664.52 Å): in the second X-shooter spectrum the peaks are at −300 and +340 km/s, the central minimum is 0.57 of the lower peak, the FWHM is 800 km/s and the profile crosses 25% of its peak at −420 and +460 km/s; in the SDSS-V coadd the peaks are at −300 and +180 km/s.

X-shooter UVB (resolving power about 5,450; `gas_disc_578709631539357440_xshooter_metals.csv`): Mg II 4481 absorption with equivalent widths of 0.131 ± 0.025 and 0.227 ± 0.014 Å in the two spectra (0.112-0.146 and 0.207-0.231 Å with three other continuum windows), and Ca II K absorption with Gaussian depths of 0.15 ± 0.02 and 0.17 ± 0.02. The Ca II K equivalent width depends on the continuum window (0.005-0.127 Å). Barycentric velocities are 13.7 ± 3.1 and 26.2 ± 1.3 km/s (Mg II) and 23.9 ± 4.6 and 26.5 ± 2.3 km/s (Ca II K).

Infrared (`gas_disc_infrared.csv`): Legacy Surveys DR10 forced photometry gives W1 = 19.1 ± 1.8 μJy and W2 = 11.1 ± 4.2 μJy, against 8.3 and 4.5 μJy for the Montreal pure-H model at 18,181 K and log g 8.15 scaled to Gaia G: 2.3 and 2.5 times the model, 5.4σ and 1.6σ with a 10% model uncertainty. A galaxy (Legacy Surveys type REX, W1 13.5 μJy) lies 3.3″ away (fracflux_w1 = 0.66). The CatWISE2020 source, 1.6″ from the white dwarf, has W1 = 35.7 ± 2.2 μJy.

Left: the Ca II triplet in the BOSS, SDSS-V, DESI and X-shooter spectra (3 Å bins). Right: equivalent width against date.

<img src="../figures/gas_discs/578709631539357440_epochs.png" width="800">

Left: the BOSS spectrum of 2010-12-05 with the SDSS-V coadd profile at scale 1. Middle: the scale of that profile in every spectrum. Right: the scale in the BOSS 2010 and SDSS 2003 spectra for 24 continuum choices.

<img src="../figures/gas_discs/578709631539357440_amplitude.png" width="800">

X-shooter 2025-04-23: the three Ca II lines in velocity (left), Mg II 4481 and Ca II K (right).

<img src="../figures/gas_discs/578709631539357440_xshooter.png" width="800">

ZTF (69 r-band points, 2019-2026; `ztf_lightcurve.py`): seasonal medians within 5% of each other (errors 1-3%); no significant period between 0.02 and 300 c/d.

## WD J1959+2208 (Gaia DR3 1827014701883095680)
G = 16.78, parallax 5.73 ± 0.07 mas; catalogued DB (SIMBAD, MWDD) and DBA/DB (SnowWhite). Both SDSS-V visits (2024-08-10 and 2024-08-11) show double-peaked Ca II emission, with peaks at about −260 and +430 km/s in the coadd. The equivalent width is 17.8 ± 0.5 Å and 26.6 ± 0.5 Å in the two visits. No other spectrum was found in SPARCL (SDSS, BOSS, DESI) or the ESO archive.

Gaia DR3 lists a second source 0.97″ away (Gaia DR3 1827014701900880256, G = 18.44, no BP−RP colour, parallax 0.20 ± 0.25 mas), inside the 2″ BOSS fibre. Its light is included in both SDSS-V spectra, so the Ca II emission is not yet attributed to the white dwarf alone; a spectrum that resolves the two sources would settle it.

Between the two visits the Ca II line flux (equivalent width times the continuum at 8600 Å) changes from 519 to 691 × 10⁻¹⁷ erg s⁻¹ cm⁻², while the continuum changes by −11% at 8600 Å and by −4% at 4500-5000 Å (`gas_disc_epoch_checks.csv`). The 2MASS and CatWISE2020 fluxes, which include the neighbour, are 2.5 (J) to 7.1 (Ks) times the Montreal pure-He model of the white dwarf (15,417 K, log g 7.57, scaled to Gaia G). A fit of two blackbodies to the Pan-STARRS1, 2MASS and CatWISE2020 fluxes ($`\chi^2`$ = 27.1 for 6 degrees of freedom, against 287.9 for one blackbody; the hot component reaches the 25,000 K limit of the fit) gives a 3,537 K component that supplies 90% of the W1 flux and has an AB magnitude of 18.6 at 0.64 μm, close to the G = 18.44 of the neighbour (`gas_disc_infrared.csv`).

<img src="../figures/gas_discs/1827014701883095680_sdssv.png" width="800">

ZTF (about 4,000 g, r and i points, 2018-2026; `ztf_lightcurve.py`): rms 1.4-2.8% against errors of 1.3-1.9%. Three nights of high-cadence data (2018-08-13, 2018-08-14 and a 6.5-hour 30-s sequence on 2024-07-05) have no 10-minute bin more than 2.5% below the median. The highest peaks of the combined periodogram are at 1, 2 and 3 c/d (daily sampling); no peak is shared by all data sets.

## GALEX J0039−0356 (Gaia DR3 2527617665632689024)
G = 18.91, parallax 2.54 ± 0.25 mas; catalogued DA (SIMBAD, MWDD, SnowWhite; SnowWhite Teff 22,800 K, log g 8.15). All four SDSS-V visits (2022-11-22 to 2024-08-08) show single-peaked Ca II triplet emission. The equivalent width is 34.4 ± 1.1 Å in the coadd and 23-49 Å in the visits; the Gaussian FWHM is 340 km/s (not corrected for the instrumental resolution) and the centroid moves between −17 and −70 km/s. H-alpha and H-beta show no emission. No other spectrum was found in SPARCL or the ESO archive.

Legacy Surveys DR10 forced WISE photometry (`wise_excess.py`) gives W1 = 7.65 ± 0.54 and W2 = 9.03 ± 1.23 nanomaggies, 7.3 and 16.3 times the Rayleigh-Jeans extrapolation of the z-band flux (excess 12σ in W1 and 7σ in W2). Two galaxies (Legacy Surveys type REX) lie 6.4″ and 7.1″ away, with W1 fluxes of 3.4 and 9.7 nanomaggies; the Tractor fit separates them (fracflux_w1 = 0.47). ZTF (about 840 points, 2018-2026) shows no dips and no significant period.

Against the Montreal pure-H model at 22,826 K and log g 8.15 scaled to Gaia G (W1 4.2, W2 2.2 μJy), the Legacy Surveys fluxes (27.8 ± 1.9 and 32.8 ± 4.5 μJy) are 6.7 and 14.7 times the model (`gas_disc_infrared.csv`). With errors from line-free windows added, the four visit equivalent widths give $`\chi^2`$ = 6.2 for 3 degrees of freedom (p = 0.10), and the 2023-10-08 value is 42.1-66.6 Å for continuum orders 1-3 (`gas_disc_epoch_checks.csv`): differences between the visits are not established.

<img src="../figures/gas_discs/2527617665632689024_sdssv.png" width="850">

## SDSS J2054+1610 (Gaia DR3 1764314497240770176)
G = 18.41, parallax 2.11 ± 0.16 mas; listed as a photometric white-dwarf candidate (SIMBAD WD?), with no earlier spectrum. The SDSS-V spectrum (two visits, 2024-06-06 and 2024-06-07) is a DA (SnowWhite Teff 27,500 K, log g 8.06) with single-peaked emission in the Ca II triplet (equivalent width 24 ± 3 Å, FWHM 265 km/s), O I 8446 and O I 7774, and no H-alpha or H-beta emission. The star is outside the Legacy Surveys DR10 footprint and not detected in CatWISE2020. The equivalent widths of the two visits agree ($`\chi^2`$ = 0.1 for 1 degree of freedom). ZTF (about 1,900 points) shows no significant period; its faint outliers are mostly at airmass 1.9-2.4 or in shallow images.

<img src="../figures/gas_discs/1764314497240770176_sdssv.png" width="850">

## WDJ1448+3225 (Gaia DR3 1283510882895711872; DESI DR1 screen)
G = 19.38, parallax 1.96 ± 0.26 mas; catalogued DBA (SIMBAD, MWDD, DESI DR1; DESI Teff 19,800 K, log g 8.32). Not in the SDSS-V SnowWhite table. The Ca II triplet shows double-peaked emission in the BOSS spectrum of 2010-03-18 and in both DESI spectra (2021-05-20 and 2022-06-05), with equivalent widths of 19.2 ± 3.8, 37.4 ± 3.7 and 23.5 ± 2.5 Å. H-alpha shows no emission. With errors from line-free windows added, the three equivalent widths give $`\chi^2`$ = 6.2 for 2 degrees of freedom (p = 0.046). Legacy Surveys DR10 forced photometry gives W1 = −0.6 ± 1.7 μJy and W2 = −5.4 ± 3.5 μJy, against 2.8 and 1.5 μJy for the Montreal pure-He model at 19,802 K and log g 8.32; the star is not detected, and no other Legacy Surveys source lies within 8″ (`gas_disc_infrared.csv`). An earlier version of this page described the WISE photometry as dominated by a neighbour; the large fracflux_w1 (10) follows from the near-zero flux.

<img src="../figures/gas_discs/1283510882895711872_epochs.png" width="600">

## WDJ1611+4017 (Gaia DR3 1379988076130545536; DESI DR1 screen; tentative)
G = 17.75, parallax 3.20 ± 0.08 mas; catalogued DA (DESI Teff 16,600 K, log g 7.55). The DESI spectrum of 2021-05-22 shows weak Ca II triplet emission (equivalent width 5.9 ± 0.4 Å; screen z_cat 16.8); the BOSS spectrum of 2012-05-25 gives 3.5 ± 1.1 Å. The single SDSS-V visit (2023-06-24) has S/N 2.6 in this region. The Legacy Surveys WISE fluxes exceed the Rayleigh-Jeans extrapolation (4× in W1, 6× in W2) but the source is blended (fracflux_w1 = 2.1).

<img src="../figures/gas_discs/1379988076130545536_epochs.png" width="600">

## DESI DR1 screen
`desi_gas_disc_screen.py` applies the same method to the 44,301 usable DESI DR1 spectra of the white dwarfs in the Amorim et al. (2026) class file (retrieved with SPARCL), with neighbours from the colours and absolute magnitudes in that file. The 99.9th percentiles of z_cat are 14.3 (DA-type), 12.5 (other white-dwarf classes), 38 (classes with an M-dwarf companion) and 104 (CVs). The known emitters give z_cat 39.7 (WD J1829+4537), 14.7 (WD J2307−0002), 12.3 (WD 0856+048) and 9.2 (WD J0857−2245). Of 46 spectra with z_cat > 8 and consistent lines, 23 have a DESI DR1 class of WD+MS, CV or STAR, and four DA with narrow Ca II and H-alpha emission are 1-2 magnitudes brighter than single white dwarfs of their colour.

## SDSS-V screen
`gas_disc_screen.py` screens the SDSS-V visit spectra of all 50,960 objects with a SnowWhite classification. The CaT region of each spectrum is divided by the median of its 100 nearest neighbours in the Gaia colour-magnitude diagram, and the residual is matched-filtered with single- and double-peaked emission templates. On the white-dwarf locus (parallax/error > 3, M_G > 8.5), the 99th and 99.9th percentiles of z_cat are 8.8 and 19.1 for DA-type spectra (30,090) and 14.7 and 27.6 for other white-dwarf classes (2,702). MS and CV classes, whose companions or accretion discs show Ca II emission, have higher values (99.9th percentiles 38 and 78).

Validation (`gas_disc_screen.csv`): of the eight known gaseous-disc white dwarfs with Ca II emission that have SDSS-V spectra, six have z_cat between 11.6 and 81.9 (WD 0842+572, SDSS J0738+1835, WD J1930−5028, WD J0529−3401, WD J1829+4537 and WD J2307−0002). SDSS J0234−0406 (z_cat 2.4) and WD 1622+587 (2.0; S/N 4.7) are not recovered. WD J0914+1914, whose disc shows H, O and S but no Ca II emission, has z_cat 3.9. The seven DESI EDR emission candidates of Ma et al. (2025) with SDSS-V spectra have z_cat 0.6 to 3.3.
