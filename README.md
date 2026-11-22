# Immune Endotypes Across Inflammatory Diseases: A Single-Cell and Spatial Analysis Portfolio

**Core question.** Using public single-cell and spatial datasets, identify the
cell states associated with disease activity and treatment response in
inflammatory diseases — and connect those states to clinical drug targets.
This is the endotype / pharmacodynamic-biomarker / patient-selection problem,
worked at cellular resolution.

Each disease module runs the **same standard pipeline** and ends with the
**same deliverable**: a translational summary figure that maps the cell states
found to the clinical targets they implicate.

## Modules (ordered by the disease areas of translational immunology teams)

| # | Module | Data (public) | Key questions |
|---|--------|---------------|----------------|
| 1 | **IBD (gastroenterology)** — core module | Ulcerative colitis & Crohn's disease single-cell atlases (Smillie et al., *Cell* 2019, UC; Martin et al., *Cell* 2019, CD) | Which cell subsets expand in inflamed tissue? Is there a cell module associated with anti-TNF non-response (inflammatory fibroblasts, OSM-pathway monocytes)? |
| 2 | **Skin: atopic dermatitis vs psoriasis** | Human skin single-cell atlas (Reynolds et al., *Science* 2021) + disease datasets | Which immune axis dominates each disease (Type 2 in AD vs IL-23/Type 17 in psoriasis)? In which cells are drug targets expressed — IL23A (risankizumab), JAK1 (upadacitinib), TNF (adalimumab)? A target-engagement / PD hypothesis at cell resolution. |
| 3 | **Spatial transcriptomics** | Public 10x Visium (IBD or skin) | Deconvolution against the module 1/2 single-cell references (cell2location / RCTD); spatial neighborhoods of inflammation (Squidpy); do the endotype cell states co-localize in tissue? |
| 4 | **RA synovium or IPF lung** (optional) | AMP RA synovial atlas / IPF single-cell atlases | Does the same pipeline recover the published cell-composition endotypes in a fourth disease area? |

## Standard pipeline (every module)

1. **QC, normalization, integration** — scVI (Python) / Harmony (R)
2. **Annotation** — automated (CellTypist / Azimuth) + manual curation against canonical immune markers
3. **Differential analysis** — pseudobulk per sample, DESeq2 (the pharma-standard design: the sample, not the cell, is the replicate)
4. **Cell–cell communication** — CellChat / LIANA
5. **Endotype & target mapping** — cell-state modules scored per sample; drug-target expression mapped per cell type; one translational summary figure per module

Module 1 is implemented in **both** Seurat (R) and Scanpy (Python);
later modules are Python-first.

## Key figures

Each module folder carries 3–5 figures with a one-line conclusion in its own
README. *(Populated as modules are completed.)*

## Reproducibility

- `envs/environment.yml` — Python environment (Scanpy stack)
- `envs/r_packages.R` — R environment (Seurat stack)
- Raw data is **not** committed; each module README gives the source and
  accession of its dataset and the download steps.

---

*Analyses of public datasets; all data credited to original authors.*
