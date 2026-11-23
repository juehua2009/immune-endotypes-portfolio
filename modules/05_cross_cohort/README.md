# Module 5 — Cross-cohort replication of a fibroblast circuit in IBD (paper module)

**Working title:** A myeloid–fibroblast cytokine circuit in inflammatory bowel disease:
cross-cohort single-cell replication and cross-species validation.

**Framing (decided Nov 22, 2026):** a practice-scale short paper, not a grand atlas.
The large integrations are published (scIBD, Nat Comput Sci 2023; gut cell atlas,
Nature 2024) and the IAF expansion story is published at Nature level (GLIS3/IAF,
2025). This module claims only what module 1 newly suggested and those papers did
not quantify: the circuit read out at three points, replicated cohort by cohort.

**The circuit (from module 1, Smillie UC):**
1. OSM (inflammatory monocytes/DCs) -> OSMR (stromal/endothelial): interaction absent
   in healthy colon, switched on in UC (LIANA strength 0.00 -> 1.06 -> 1.16).
2. IL11 -> IL11RA fibroblast loop: detectable only in inflamed tissue.
3. IL1B -> SIGIRR (the inhibitory IL-1 receptor): strongest in healthy colon and
   weakening with inflammation — a loss-of-braking gradient (20.49 -> 10.45 -> 5.88).

**Datasets:**
| Cohort | Content | Source | Status |
|---|---|---|---|
| Smillie et al., Cell 2019 | Human UC colon, 365k cells | SCP259 / scverse h5ad | DONE (module 1) |
| Martin et al., Cell 2019 | Human Crohn's ileum (GIMATS / anti-TNF resistance) | GEO GSE134809 | VERIFIED (RAW.tar 579M; 22 ileal samples = 11 involved + 11 uninvolved, raw MTX counts; 9 PBMC samples in the series are excluded — different compartment, later batch; no healthy controls, contrast = involved vs uninvolved) |
| Kinchen et al., Cell 2018 | Human healthy + UC colon; mouse healthy + DSS colon | GEO GSE114374 | VERIFIED (4 matrices, see scripts/) |
| Jasso et al., PLoS Biol 2022 | Mouse chronic DSS, stroma only (IL-11+ fibroblasts) | GEO GSE172261 / SCP1711 | optional second mouse arm |

Kinchen is the keystone: one series carries human UC AND mouse DSS from the same
lab and pipeline, so the cross-species check needs no extra cohort.

**Cohort roles (set after Kinchen notebook 01, Nov 22, 2026):** the Kinchen human
matrices are stromal/mesenchymal-enriched — fibroblasts and endothelium dominate,
immune and epithelial cells are sparse (MS4A1/MUC2/CHGA absent from the matrix),
and the values are pre-normalized (no raw counts, so no DESeq2). Roles:
Smillie = full circuit (has every compartment); Martin = full-circuit replication
in Crohn's (immune-rich lamina propria, raw counts in the RAW.tar); Kinchen =
fibroblast-side readouts (OSMR/SIGIRR/CHI3L1 within fibroblasts, by condition)
plus the mouse DSS cross-species arm. Kinchen notebook 01 first pass: OSMR and
CHI3L1 sit in the expected fibroblast clusters and SIGIRR is broadly detectable,
but IL11/IL13RA2 do not resolve at whole-data resolution 0.5 — the inflammatory
subset is diluted inside large fibroblast clusters, so notebook 02 subclusters
the fibroblast compartment before any replication verdict.

**Analysis plan:**
1. Per-cohort QC with the module 1 template (same gates, same readouts).
2. Human integration: Smillie + Martin + Kinchen-human, scVI with dataset as batch;
   label harmonization to the Smillie Major levels first, fine labels by transfer.
3. Circuit replication, per cohort separately (each cohort = one replication):
   inflammatory fibroblast abundance, OSM/IL11 program scores, pseudobulk DE of the
   watch genes, LIANA watch pairs incl. IL1B-SIGIRR.
4. Cross-species: score the human fibroblast program in Kinchen mouse DSS vs mouse
   healthy (ortholog-mapped gene symbols); optional Jasso stroma as second arm.
5. One summary figure per readout, forest-style across cohorts.

Analyses of public datasets; all data credited to original authors.

**Results so far (Nov 22, 2026):**
- Kinchen nb01: stromal-enriched data; OSMR/CHI3L1 in the expected fibroblast
  clusters; SIGIRR broadly detectable; IL11/IL13RA2 unresolved at res 0.5.
- Kinchen nb02 (fibroblast subclustering, 17 subclusters): REPLICATION CONFIRMED.
  Subcluster 6 is the canonical inflammatory fibroblast — the only subcluster
  co-expressing IL11, IL13RA2 and IL24, with CHI3L1/MMP3/IL6/CXCL1 high — and it
  is essentially UC-specific (~14% of UC fibroblasts vs ~1% of healthy).
  Healthy fibroblasts are dominated by ADAMDEC1/WNT2B homeostatic subclusters.
  Mirrors Smillie (inflammatory fibroblasts ~0 in healthy, expanded in disease).

- Martin nb03 (Nov 22, 2026): loaded clean — 22 ileal samples, 15 clusters, full
  lamina propria (epithelial 5/6, T 0/1/8/9, B 2, plasma 3/7/12/13, myeloid 4,
  fibroblast 10, endothelial 11). The circuit's two ends sit in the right
  compartments: myeloid cluster 4 is the OSM+/IL1B+ source (IL1B strongest in the
  dataset); fibroblast cluster 10 carries OSMR/CHI3L1/IL1R1/MMP3 with only IL11/
  IL13RA2 specks at res 0.5 — same dilution pattern as Kinchen, so notebook 04
  subclusters Martin fibroblasts and checks involved-enrichment, plus the
  cluster x condition composition (GIMATS components) across the whole data.

- Martin nb04 (Nov 22, 2026): fibroblast subclustering (20 subclusters) finds
  subcluster 7 = the canonical inflammatory fibroblast with the FULL signature
  (IL11 + IL13RA2 + IL24 + CHI3L1 + MMP3 + CXCL1 + IL1R1 + PDPN) — the richest
  IAF readout of the three cohorts. CHI3L1-high subclusters 6/15/3 occupy an
  almost purely Involved fibroblast territory on the UMAP. Whole-data
  composition: myeloid cluster 4 and B cluster 2 skew to Involved (~1.75x and
  ~2.4x by bar share), consistent with GIMATS; total fibroblast share is flat —
  the change is within-fibroblast state composition, as in Smillie and Kinchen.
  State-level replication now confirmed in all three human cohorts; remaining
  module work: Martin LIANA (nb05), human integration, Kinchen mouse arm, and
  the cross-cohort summary figure with exact enrichment numbers.

- Martin nb05 LIANA (Nov 22, 2026): OSM-OSMR REPLICATED as a switch — Involved
  1.43 (Myeloid -> Fibroblast), Uninvolved 0.00; OSM-LIFR and TNF-TNFRSF1A also
  Involved-only. IL11-IL11RA absent at compartment resolution (IL11 is confined
  to fibroblast subcluster 7, diluted by pooling — resolution artifact, tested
  at subcluster resolution in nb05b). IL1B-SIGIRR present in both states at
  similar strength (3.02 vs 2.66, different targets) — the Smillie brake
  gradient does NOT replicate in the involved-vs-uninvolved contrast, which has
  no healthy arm; the SIGIRR claim rests on Smillie + Kinchen, reported as is.

- Martin nb05b (fibroblast-split LIANA, Nov 22, 2026): the OSM switch sharpens —
  Involved best pair is Myeloid -> Fib_CHI3L1 at 2.39 (vs 1.43 pooled), absent in
  Uninvolved. The label counts are themselves a headline: Fib_IAF = 148 Involved
  vs 1 Uninvolved; Fib_CHI3L1 = 445 vs 2 — inflammatory fibroblast states make
  up 44% of involved fibroblasts vs 0.2% of uninvolved, a near-binary expansion.
  IL11-IL11RA stays below LIANA detection (the IAF group is only 148 cells);
  the loop is claimed at expression level for Martin, communication level for
  Smillie.

- Kinchen mouse nb06 (Nov 22, 2026): CROSS-SPECIES REPLICATION — the four
  highest-program fibroblast clusters (1, 2, 6, 7; program scores 0.55-1.05)
  are entirely DSS-specific (0.000 share in healthy) and together make up ~74%
  of DSS fibroblasts; Il11 peaks exactly in clusters 1/2; Osm sits only in the
  small immune cluster 13, Osmr across fibroblasts. Caveats recorded: Chi3l1 is
  absent from the mouse matrices (authors' filtering; score ran on 8/9 genes),
  and condition is confounded with file in Kinchen's deposit (healthy and DSS
  are separate matrices) — stated in methods. Analysis core of the paper is now
  complete; remaining: summary figure notebook (exact numbers, one uniform IAF
  rule across cohorts) and the manuscript draft.

- Cross-cohort summary (nb07, Nov 22, 2026) — FINAL NUMBERS (inflammatory
  fibroblast share of fibroblasts, one uniform rule):
  Smillie UC: Healthy 0.24% / Non-inflamed 7.08% / Inflamed 12.90%.
  Kinchen UC (subcluster 6): Healthy 0.03% / UC 13.94%.
  Martin CD: IAF 0.06% -> 11.06%; CHI3L1-high 0.12% -> 33.26%;
  combined 0.19% -> 44.32% (Uninvolved -> Involved).
  Kinchen mouse: Healthy 0.00% / DSS 73.82%.
  Saved: results/cross_cohort_summary.csv + figures/cross_cohort_iaf_share.png.
  ANALYSIS COMPLETE. Remaining: manuscript draft only.

**Response-validation design (added Nov 22, 2026):** single-cell anti-TNF
response data are genuinely scarce (Martin is effectively the only cohort).
The field-standard cross-validation route is pretreatment bulk biopsy cohorts
with mucosal-healing response labels — verified on GEO Nov 22, 2026:
GSE12251 (UC, infliximab, response at 8 wk), GSE16879 (UC + CD, infliximab,
4-6 wk), GSE73661 (UC, infliximab subset, 4-6 wk), GSE23597 (UC, infliximab).
Plan (notebook 09): score the IAF program (9-gene set + the DE-derived
signature from notebook 08A) in each pretreatment biopsy and test
responder vs non-responder separation per cohort. This is the route West 2017
(OSM) and Martin 2019 (GIMATS bulk validation) used. Note: generic bulk
non-response signatures are occupied territory (e.g. a 2026 68-gene consensus
analysis of the same cohorts; an IGFBP5 fibroblast multi-omics paper) — our
claim stays cell-state-specific: the signature of one defined fibroblast
state, discovered at single-cell resolution.

- Notebook 08 results (Nov 22, 2026): (A) IAF molecular identity — pseudobulk
  DESeq2 IAF vs other fibroblasts (Martin, Involved): 25 genes up at LFC>1
  (MMP1 +7.2, SAA1 +8.0, IL24 +5.6, MASP1, WFDC2, CXCL5, CSF3, MMP3, CEMIP,
  TWIST2 — matrix proteases + acute phase + neutrophil chemotaxis) and a
  STRONGER down side (MGP, NDUFA4L2, PPP1R14A, SPARCL1, IGFBP7, GPX3 —
  homeostatic fibroblast identity lost). The 25-gene up set is the transfer
  signature for the bulk response validation (notebook 09). (B) Per-patient
  burden is subgroup-concentrated, not graded: Martin Involved median 0.000 /
  mean 0.276 / max 0.837 (n=11); Smillie Inflamed median 0.051 / max 0.535;
  Kinchen libraries are pooled (n=2) — footnote only. (C) OSM+ myeloid fraction
  is dilution-confounded (Martin 9.9% vs 6.2%; Smillie reversed) — NOT used;
  ligand side rests on inflammatory-monocyte expansion + receptor switch.
  (E) IAF burden does NOT correlate with GIMATS components within Martin
  patients (plasma r=0.19, myeloid 0.03, T -0.18; n=11) — the state is an
  independent stromal axis, not a GIMATS passenger. (D) decoupler PROGENy/
  DoRothEA pending an API fix for the installed version.
  (D) Drivers (decoupler 2.2 PROGENy/DoRothEA, Martin fibroblasts, Involved):
  TNFa +4.82, NFkB +4.04, JAK-STAT +3.10 are the top pathway activities;
  TGFb is DOWN (-0.75) - an inflammatory, not classic-fibrotic, state.
  TF activities agree: RELA +3.03, NFKB1 +2.81, FOSL1/JUN (AP-1), IRF1,
  STAT2 up; SRF -1.44, ELK1, MEF2C down (homeostatic program switched off);
  TWIST2 is up both as a gene (DE) and as TF activity. Discussion point:
  co-driving by JAK-STAT (gp130/OSM family) alongside TNF offers a mechanism
  for anti-TNF insufficiency in high-burden patients. Notebook 08 COMPLETE.
  Remaining: notebook 09 (bulk response validation) + manuscript draft.
