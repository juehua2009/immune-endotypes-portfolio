# A myeloid–fibroblast cytokine circuit in inflammatory bowel disease: cross-cohort single-cell replication and cross-species validation

Juehua Yu

*Draft — Methods and Results only (Introduction and Discussion to follow). Analyses of public datasets; all data credited to original authors.*

## Methods

### Data sources

Three published human single-cell RNA-sequencing cohorts and one mouse arm were analysed, each cohort processed independently and treated as one replication unit; datasets were never merged into a single integrated object for the primary analyses.

(i) **Smillie et al. (2019)**, ulcerative colitis (UC) colon atlas: the processed count matrix (365,492 cells × 18,172 genes) was obtained from the scverse example-data server. The matrix carries the original authors' annotations (51 cell subsets with hierarchical major labels) and sample metadata. The cohort comprises 12 healthy subjects and 18 UC subjects, the UC subjects contributing paired inflamed and non-inflamed biopsies. (ii) **Martin et al. (2019)**, Crohn's disease (CD) ileum (GEO GSE134809): per-sample count matrices (matrix market format) for the 22 ileal tissue samples (11 involved, 11 uninvolved) were assembled from the series RAW archive; the 9 peripheral-blood samples in the series were excluded (different compartment, separate batch). The series has no healthy arm; the contrast is involved versus uninvolved ileum. (iii) **Kinchen et al. (2018)** (GEO GSE114374): the deposited human healthy-colon and UC expression matrices were used as provided. These matrices are stromal/mesenchymal-enriched (fibroblast and endothelial clusters dominate; MS4A1, MUC2 and CHGA are absent) and contain pre-normalized expression values rather than raw counts; the cohort therefore contributes abundance, program-score and communication readouts but no count-based differential expression. (iv) **Mouse arm**: the mouse healthy and dextran sulfate sodium (DSS) colitis matrices from the same series (GSE114374) were analysed for cross-species validation (3,760 healthy and 3,491 DSS cells before quality control).

### Preprocessing and quality control

Processing followed one template for every cohort. For count data (Smillie, Martin), raw counts were retained in a separate layer, cells with fewer than 200 detected genes were removed, libraries were normalized to 10,000 counts per cell and log-transformed. For the Smillie atlas this left 351,952 cells. Highly variable genes (3,000, with condition as batch key) were used for principal-component analysis, neighbourhood graph construction, UMAP embedding and Leiden clustering (resolution 0.5). For the Smillie atlas, a scVI model (3,000 highly variable genes selected with the seurat_v3 flavour, sample as batch key, 30 latent dimensions) provided a batch-corrected embedding as a clustering cross-check (28 clusters, concordant with the 29 Leiden clusters from the PCA embedding and with the authors' labels). The pre-normalized Kinchen matrices were used as provided after the same cell-level quality gates. Martin samples were assembled per sample, quality-controlled with the same gates, and clustered jointly (15 clusters), with all 22 samples represented across clusters.

### Annotation

For Smillie, the authors' published subset labels were used directly. For Martin and Kinchen, clusters were annotated from canonical markers (EPCAM, COL1A1/DCN, PECAM1, PTPRC, CD3D, MS4A1, JCHAIN, LYZ2, S100A8). In Martin, this resolved epithelial (clusters 5, 6), T (0, 1, 8, 9), B (2), plasma (3, 7, 12, 13), myeloid (4), fibroblast (10) and endothelial (11) compartments.

### Fibroblast subclustering and the uniform inflammatory-fibroblast rule

Because the inflammatory fibroblast state is diluted at whole-data resolution, the fibroblast compartment of each cohort was subset and re-clustered: Smillie by the authors' stromal labels; Kinchen human by whole-data clusters in which more than half of cells expressed COL1A1, re-embedded and re-clustered (Leiden resolution 1.0, 17 subclusters); Martin cluster 10, re-clustered (20 subclusters); Kinchen mouse by the same COL1A1 rule. One uniform rule then defined the inflammatory fibroblast (IAF) state in every cohort: within fibroblasts, the label or subcluster with the highest mean inflammatory-program score (below). For the mouse arm, all fibroblast clusters with mean program score above 0.5 were combined. State abundance is reported throughout as the share of fibroblasts (Smillie: share of the stromal compartment, which is fibroblast-dominated) per condition.

### Program scoring

The inflammatory fibroblast program is the nine-gene set IL11, IL13RA2, IL24, CHI3L1, MMP3, MMP10, IL6, CXCL1, CXCL5, scored per cell with scanpy's `score_genes`. The myeloid ligand-side program (OSM, OSMR, IL1B, IL6, TNF, S100A8, S100A9) was scored the same way. For the mouse arm, mouse orthologue symbols were used (Il11, Il13ra2, Il24, Mmp3, Mmp10, Il6, Cxcl1, Cxcl5); CHI3L1/Chi3l1 is absent from the deposited mouse matrices, so the mouse score uses eight of the nine genes. Human and mouse matrices were never merged; cross-species comparison is at program level only.

### Pseudobulk differential expression

Count-based differential expression used DESeq2 (via pydeseq2) on pseudobulk profiles: counts summed per subject × condition × cell label (Smillie) or per sample × state (Martin), groups with fewer than 20 cells dropped, genes retained at ≥10 counts in ≥3 samples, and a label tested only when at least 3 samples were available per arm. Smillie contrasts were inflamed versus healthy per label and per major compartment. The IAF signature was derived in Martin: IAF (subcluster 7) versus all other fibroblasts within involved samples, using only patients contributing both states (paired design). Kinchen data were excluded from differential expression (no raw counts).

### Cell–cell communication

Communication was inferred per condition with LIANA (rank-aggregate method, consensus resource, minimum expression proportion 0.1, 100 permutations), on label-stratified subsamples capped at 1,000 cells per condition × label. A pre-specified watch list (OSM–OSMR, OSM–LIFR, IL11–IL11RA, IL1B–SIGIRR, IL1B–IL1R1, TNF–TNFRSF1A, CCL2–CCR2, IL6–IL6R) was extracted with complex-aware matching (ligand exact; receptor matched as substring, because the resource names receptor complexes, e.g. OSMR as IL6ST_OSMR). Strength is reported as the inverse magnitude rank (1/magnitude_rank; higher = stronger), best source–target pair per condition. In Martin, communication was computed at compartment level and again after splitting fibroblasts into IAF (subcluster 7), CHI3L1-high (subclusters 6, 15, 3) and other fibroblasts.

### Pathway and transcription-factor activity

PROGENy pathway activities (top 500 responsive genes per pathway) and DoRothEA transcription-factor activities (confidence levels A and B) were computed with decoupler (version 2.2.0) using the ULM method (minimum 3 targets) on log-normalized expression in Martin fibroblasts, and compared as group-mean differences between IAF and other fibroblasts within involved samples.

### Per-patient burden and module correlations

Per-patient IAF burden is the share of a patient's fibroblasts assigned to the inflammatory state(s), computed per cohort from the objects above (Smillie by subject; Martin by sample; Kinchen human by library — the deposited Kinchen libraries pool donors, n = 2 per condition, and are reported for completeness only). In Martin, per-patient composition shares of plasma cells, myeloid cells and T cells were computed from whole-data clusters, and Pearson correlations with IAF burden (share of all cells) were calculated across the 11 involved samples.

### Bulk-cohort response validation

Four public pretreatment bulk-biopsy cohorts with infliximab response labels were scored with the IAF signatures: GSE12251 (UC; week-8 histologic mucosal-healing response, all biopsies pretreatment), GSE16879 (UC and CD; biopsies before the first infliximab infusion, response at 4–6 weeks), GSE23597 (UC; week-0 biopsies from the non-placebo arms, week-8 response) and GSE73661 (UC, infliximab arm; week-0 biopsies, response derived from the same individual's week-4/6 Mayo endoscopic subscore, ≤1 responder, ≥2 non-responder). Series matrices were obtained from GEO; probe-to-gene mapping was extracted from the GEO platform tables (GPL570 for GSE12251/GSE16879/GSE23597; GPL6244 for GSE73661, on which TWIST2 has no mapped transcript, so SIG25 uses 24 genes there). Probes were averaged per gene, each gene was z-scored across the analysis subset, and the signature score is the mean z-score over available signature genes. Responders and non-responders were compared by Mann-Whitney U (one-sided, non-responders higher); discrimination is reported as the area under the ROC curve with non-response as the positive class. GSE16879 was analysed at biopsy level and stratified by disease.

### Software

Analyses were performed in Python 3.11 with scanpy 1.11, anndata 0.11, scvi-tools 1.4.2, pydeseq2 0.4, liana 1.10.0 and decoupler 2.2.0 (numpy 1.26, scipy 1.12). All code and result tables are available in the accompanying repository.

## Results

### One rule recovers the same fibroblast state in three human cohorts

Applying the uniform rule — the fibroblast label with the highest inflammatory-program score — recovered a coherent inflammatory fibroblast state in each cohort without manual selection. In Smillie, the rule selected the original authors' own "Inflammatory Fibroblasts" label, validating the rule against independent annotation. In Kinchen human colon, fibroblast subcluster 6 was the only subcluster co-expressing IL11, IL13RA2 and IL24, with CHI3L1, MMP3, IL6 and CXCL1 high, while healthy fibroblasts were dominated by ADAMDEC1/WNT2B homeostatic subclusters. In Martin Crohn's ileum, subcluster 7 carried the full signature (IL11, IL13RA2, IL24, CHI3L1, MMP3, CXCL1, IL1R1, PDPN), and three further subclusters (6, 15, 3) formed a CHI3L1-high inflammatory territory; the two groups are reported separately and combined below.

### The state expands in disease in every cohort, by state switching rather than compartment growth

IAF share of fibroblasts was near zero in every baseline and rose sharply in disease (Table 1, Fig. 1). In Smillie UC: 0.24% in healthy colon, 7.08% in non-inflamed UC and 12.90% in inflamed UC. In Kinchen: 0.03% in healthy versus 13.94% in UC. In Martin: IAF (subcluster 7) 0.06% in uninvolved versus 11.06% in involved ileum; CHI3L1-high fibroblasts 0.12% versus 33.26%; combined 0.19% versus 44.32%. In Martin the total fibroblast share of all cells was unchanged between conditions — the expansion is a change of fibroblast state composition, not of compartment size.

Transcriptional analysis in Smillie corroborated the abundance readout at compartment level. Myeloid pseudobulk (inflamed versus healthy) yielded 1,550 differentially expressed genes, including S100A8 (log2FC +4.60, padj 0.003), S100A9 (+2.00, padj 0.004) and CHI3L1 (+3.79, padj 0.035). Stromal pseudobulk yielded 3,568 differentially expressed genes, with IL13RA2 (+6.33, padj 0.0016), CHI3L1 (+5.15, padj 0.0006) and the OSM receptor OSMR (+0.65, padj 0.013) significantly up in inflamed tissue. The inflammatory fibroblast and inflammatory monocyte labels themselves could not be tested in the healthy contrast because both populations are essentially absent from healthy colon (fewer than three healthy pseudobulk groups) — their disease association is carried by their appearance, not by within-label expression change. Consistently, the direct inflamed-versus-non-inflamed pseudobulk contrast within these labels returned no differentially expressed genes: between UC states, the program is already switched on, and what changes is how many cells carry it. Inflammatory monocytes showed the same behaviour as the fibroblast state — absent in healthy colon, expanded in UC, with their OSM/IL1B program score already high in non-inflamed UC mucosa.

### IAF burden is concentrated in a patient subgroup

Per-patient burden was not a graded function of disease status (Fig. 2). In Martin involved ileum (n = 11), median burden was 0.000 while the mean was 0.276 and the maximum 0.837: half of the patients carried essentially no IAF, while in the highest-burden patient more than four-fifths of fibroblasts were in the inflammatory states. In Smillie inflamed UC (n = 18), median burden was 0.051 with a maximum of 0.535, against a healthy maximum of 0.010. The state therefore behaves as a patient-level variable with a subgroup distribution — an endotype-like pattern — rather than as a uniform readout of inflammation.

### The myeloid–fibroblast circuit replicates at communication level, with one non-replication

In Smillie, LIANA recovered the circuit as a switch: OSM–OSMR strength (myeloid sources to stromal and endothelial targets) was 0.00 in healthy colon, 1.06 in non-inflamed and 1.16 in inflamed UC; the IL11–IL11RA fibroblast loop was detectable only in inflamed tissue (2.24); IL1B–SIGIRR, the inhibitory IL-1 receptor interaction, was strongest in healthy colon (20.49) and weakened with inflammation (10.45, then 5.88), a loss-of-braking gradient; TNF–TNFRSF1A was constitutive (2.4–2.7 in all conditions). In Martin, the OSM switch replicated with the correct anatomy: OSM–OSMR was absent in uninvolved ileum and present in involved ileum (1.43, myeloid → fibroblast), sharpening to 2.39 for the CHI3L1-high fibroblast target when fibroblasts were split; OSM–LIFR (myeloid → endothelium, 1.00) and TNF–TNFRSF1A (1.64) were likewise involved-only. Two readouts did not replicate and are reported as such. IL11–IL11RA remained below LIANA detection in Martin at both resolutions (the IAF group contains only 148 cells), so the IL11 loop is supported at communication level in Smillie and at expression level in Martin. The IL1B–SIGIRR brake gradient did not replicate in Martin (3.02 involved versus 2.66 uninvolved, different dominant targets); the Martin contrast lacks a healthy arm, and the gradient claim therefore rests on Smillie.

### Molecular identity: an NF-κB-driven inflammatory program with loss of homeostatic identity

Pseudobulk differential expression between IAF and other fibroblasts within involved Martin samples (paired by patient) defined the state's molecular identity (Fig. 3). Twenty-five genes were up-regulated at log2FC > 1 (FDR < 0.05), headed by SAA1 (+7.99), MMP1 (+7.19), IL24 (+5.58), MASP1 (+5.09) and WFDC2 (+4.96), and including CXCL5, CSF3, MMP3, CEMIP and TWIST2 — a matrix-protease, acute-phase and neutrophil-chemoattractant program. The down-regulated side was stronger still: homeostatic fibroblast identity genes MGP (−5.70), NDUFA4L2 (−7.85), IGFBP3 (−6.65), PPP1R14A (−6.46), COX4I2 (−6.17), GPX3 (−5.07), MCAM (−5.05), SPARCL1 (−4.14) and ADAMTS1 (−4.06). Pathway activities (PROGENy, IAF minus other fibroblasts) were dominated by TNFα (+4.82), NFκB (+4.04) and JAK-STAT (+3.10), with TGFβ activity reduced (−0.75) — an inflammatory rather than classic fibrotic state. Transcription-factor activities (DoRothEA) agreed: RELA (+3.03), NFKB1 (+2.81), FOSL1 (+1.94), IRF1 (+1.78), STAT2 (+1.65) and JUN (+1.37) up; the homeostatic regulators SRF (−1.44), ELK1 (−1.36) and MEF2C (−1.21) down. TWIST2 was up-regulated both as a gene and as an inferred activity.

### The stromal state is independent of the GIMATS module

Because the Martin cohort defined the GIMATS resistance module, we tested whether IAF burden simply tracks that module within patients. Across the 11 involved samples, IAF burden (share of all cells) correlated weakly or not at all with plasma-cell share (Pearson r = 0.19), myeloid share (r = 0.03) and T-cell share (r = −0.18). The fibroblast state varies independently of the module's immune components and constitutes a separate, stromal axis of variation in the same lesions.

### Cross-species validation in mouse DSS colitis

In the Kinchen mouse arm, the same rule recovered the state: the four fibroblast clusters with the highest program scores (clusters 1, 2, 6, 7; scores 0.51–1.05) were entirely DSS-specific (0.000 share in healthy colon) and together accounted for 73.82% of DSS fibroblasts versus 0.00% in healthy. Il11 expression peaked in clusters 1 and 2; Osm was confined to the small immune cluster, with Osmr distributed across fibroblasts — the same source–target anatomy as in human. Two caveats apply: Chi3l1 is absent from the deposited mouse matrices (the score used eight of nine genes), and in this deposit condition is confounded with file (healthy and DSS are separate matrices from the same study and pipeline).

### The IAF program predicts anti-TNF non-response in pretreatment bulk biopsies

Finally, we asked whether the state has a measurable consequence in bulk tissue: whether its program, scored in pretreatment biopsies, separates patients who will fail infliximab from those who will respond. Both signatures scored higher in non-responders in every cohort stratum tested (Fig. 4). The 25-gene DE signature gave AUCs of 0.93 in GSE12251 (UC, n = 23; one-sided p = 2.5 × 10⁻⁴), 0.74 in GSE16879 UC (n = 24; p = 0.030), 0.83 in GSE16879 CD (n = 37; p = 3.4 × 10⁻⁴), 0.71 in GSE23597 (UC, n = 32; p = 0.051) and 0.73 in GSE73661 (UC, n = 23; p = 0.038). The compact 9-gene program performed at least as well in four of six strata (AUC 0.86, 0.95, 0.82, 0.77 and 0.87 respectively; GSE16879 pooled AUC 0.88, p = 2.1 × 10⁻⁷). The single-cell burden result therefore has a bulk-level counterpart: pretreatment mucosal abundance of this one fibroblast state's program tracks anti-TNF non-response across independent UC and CD cohorts. Two strata are small (GSE23597 includes only 7 non-responders) and GSE16879 is analysed at biopsy level; the claim rests on the consistency of direction and magnitude across cohorts, not on any single one.

**Table 1.** Inflammatory fibroblast share of fibroblasts by cohort and condition (uniform rule; see Methods). Source: results/cross_cohort_summary.csv.

**Figure 1.** IAF share of fibroblasts across cohorts (figures/cross_cohort_iaf_share.png).
**Figure 2.** Per-patient IAF burden distributions (figures/per_patient_iaf_burden.png).
**Figure 3.** Pseudobulk DE, IAF versus other fibroblasts, Martin involved (figures/de_iaf_volcano.png).
**Figure 4.** Bulk response validation: signature AUCs for anti-TNF non-response by cohort (figures/bulk_response_auc.png; values in results/bulk_response_validation.csv).
