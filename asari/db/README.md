# Databases for testing purpose

For more comprehensive databases, please retrieve the full copies to use with Asari.

- core_blood_metabolites_20250212.json, examples from CSM
- gcms_corelib_20260930.msp, 501 example spectra of high-resolution predicted GC-MS library
- lowres_test_gcms.msp, 1000 example low-resolution spectra from MoNA
- v3_authlib2024_hilicpos.json, LC-MS authentic standards from Li lab
- calibration_mz.json, reference m/z values for future m/z calibration, as {"pos": [...], "neg": [...]}.
  Combined from (1) [M+H]+ / [M-H]- of compounds in KEGG global metabolic map (map01100, formulas with C,H,N,O,P,S,F,Cl,Br,I,Se only; m/z 50-2000),
  and (2) top 2000 CSM r1 KDE consensus mass registries (by KDE_segment_density) per mode.
  Values within 2 ppm are merged (theoretical KEGG value kept if present), rounded to 6 decimals.
