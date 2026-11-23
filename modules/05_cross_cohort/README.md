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
