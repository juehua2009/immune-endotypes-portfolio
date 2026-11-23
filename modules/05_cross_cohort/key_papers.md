# Five similar-idea papers (verified Nov 22, 2026) — positioning set for the manuscript

1. **West NR et al. Oncostatin M drives intestinal inflammation and predicts
   response to tumor necrosis factor-neutralizing therapy in patients with
   inflammatory bowel disease. Nature Medicine 2017;23(5):579-589.**
   PMID 28368383; DOI 10.1038/nm.4307.
   Idea: OSM/OSMR high in inflamed IBD mucosa; OSMR sits on stromal cells;
   high pretreatment OSM predicts anti-TNF failure in >200 patients incl. two
   phase 3 trial cohorts. Relation: our circuit's founding paper and the
   archetype of bulk-cohort response validation. We add the single-cell,
   cross-cohort, cross-species circuit replication (OSM-OSMR switch at
   communication level) that 2017 bulk data could not resolve.

2. **Martin JC et al. Single-cell analysis of Crohn's disease lesions
   identifies a pathogenic cellular module associated with resistance to
   anti-TNF therapy. Cell 2019;178(6):1493-1508.e20.**
   Idea: GIMATS cellular module (IgG plasma + inflammatory mononuclear
   phagocytes + activated T + stromal) in involved lesions, validated against
   anti-TNF resistance in bulk cohorts. Relation: our Martin data source and
   the module-to-bulk template. Our addition: the stromal arm resolved to one
   fibroblast state, shown to be an INDEPENDENT axis (burden does not
   correlate with GIMATS components within patients, r <= 0.19).

3. **Friedrich M et al. IL-1-driven stromal-neutrophil interactions define a
   subset of patients with inflammatory bowel disease that does not respond
   to therapies. Nature Medicine 2021;27:1970-1981.**
   DOI 10.1038/s41591-021-01520-5.
   Idea: bulk + single-cell + histopathology across 3 cohorts (n=376) define
   pathotypes; the ulcer-bed pathotype = neutrophil infiltration + activated
   fibroblasts (neutrophil-chemoattractant, IL-1R-dependent, not TNF);
   signatures increased in nonresponders across 4 independent cohorts (n=343).
   Relation: the closest conceptual sibling — our IAF signature (CXCL5, CSF3,
   IL1R1) is the same biology at state resolution; we replicate per cohort
   with exact burden numbers and add the mouse arm.

4. **Multi-omics reveals fibroblast subtypes with distinct roles in
   inflammatory bowel disease. Biology Direct 2026.**
   DOI 10.1186/s13062-026-00853-w.
   Idea: integration of 4 GEO scRNA datasets (incl. GSE134809; 165,618 cells)
   + spatial; ICAM1+ inflammatory and FAP+ matrix fibroblasts; TF activity;
   CellChat; fibroblast signatures in bulk predict anti-TNF poor response
   (ML), with mIHC validation. Relation: CLOSEST COMPETITOR — overlapping
   data and claims. Our differentiation: per-cohort replication with one
   uniform rule (not one integrated atlas), exact per-patient burden
   distributions (endotype), the OSM/IL11/SIGIRR circuit readouts at
   communication level, cross-species DSS validation, GIMATS-independence.

5. **A reproducible pretreatment mucosal-inflammatory-remodeling state is
   associated with induction-phase anti-TNF non-response in ulcerative
   colitis. Frontiers in Genetics 2026.**
   Idea: 68-gene consensus signature from three pretreatment infliximab UC
   cohorts (GSE12251, GSE16879, GSE23597; total n=79), leave-one-dataset-out
   validation, pooled Hedges g = 1.36. Relation: occupies the generic bulk
   signature space on exactly the cohorts planned for our notebook 09. Our
   differentiation: the signature tested is one defined cell state's program
   (25-gene IAF DE signature), derived at single-cell resolution, not a bulk
   DE consensus.

Related (6th): Multi-omics analysis identifies fibroblast IGFBP5 as a key
target of anti-TNF therapy in ulcerative colitis. Frontiers in Immunology
2026 (fimmu.2026.1841149) — fibroblast gene + spatial enrichment + bulk
validation in GSE12251/GSE16879/GSE73661; cite in discussion.
