# TAURUS resolution-control and residual-state analysis plan
Locked before any computation under this plan. Companion to taurus_analysis_plan_2026-11-24.md (primary plan, commit b00f686, tag taurus-plan-v1). This plan covers NEW questions raised after the primary result; nothing here alters the primary plan or its reported result.

## Question
Does single-cell resolution add information over a transcript score computed from the SAME cells of the SAME biopsies? All comparisons below are within-patient, within-biopsy: single-cell state fractions versus pseudobulk scores aggregated from the identical cell sets. The only thing that varies is resolution.

## Data and frozen definitions
- Cohort: TAURUS (Thomas et al., Nat Immunol 2024), adalimumab, biologic-naive; primary analysis set as in the primary plan (36 patients with Pre samples and remission labels; 34 paired).
- Compartments: authors' per-compartment h5ad files (Zenodo 10.5281/zenodo.14007626): fibroblast/pericyte, vascular, ileal epithelial, colonic epithelial, myeloid, ILC, CD4 T, CD8 T, B cell, plasma cell. Cell-state labels = authors' annotations (major/minor/sub_bucket), used as deposited.
- Frozen ruler: SIG8 AUCell (tmin = 1), universe 16,632, T = 0.107229, no recalibration — IAF fraction exactly as in the primary plan.
- Pseudobulk: per patient per timepoint, counts summed over ALL cells of the biopsy (all compartments pooled), CPM -> log1p; SIG8 score = mean over the 8 genes; INFLAM6 score = mean over the 6 genes. Compartment-restricted pseudobulks are sensitivity analyses only.
- Post-treatment timing is heterogeneous (8 weeks to 1.5 years after treatment start, per the published record). All change analyses under this plan are therefore ASSOCIATIONS WITH REMISSION STATUS, never predictions of later outcome. This limitation is stated wherever a change result is reported.

## Analysis A — baseline prediction at two resolutions (one primary comparison)
- A1 (primary): AUROC of the full-biopsy pseudobulk SIG8 score for Non_Remission vs Remission at baseline, with the single-cell IAF-fraction AUROC (0.4875, from the primary plan) as the fixed comparator. Report both with bootstrap 95% CIs (10,000 resamples, seed fixed).
- Interpretation fixed in advance: if pseudobulk AUROC is materially higher (CI excludes 0.4875), the baseline predictive signal in bulk-type data is a property of aggregate transcript abundance, not of the fibroblast-state fraction. If comparable, the bulk-cohort association (pooled AUC 0.826) does not reproduce within TAURUS at either resolution and is attributed to cohort differences (drug mix, population, sampling site).
- Disclosure: a fibroblast/pericyte-restricted pseudobulk version was computed on 2026-11-25 (step7d) BEFORE this plan was locked, as a post-hoc bridge requested in review. Its output has NOT been examined at lock time. It is reported, labelled post hoc and compartment-restricted, as a sensitivity analysis only; A1 above is the primary test.

## Analysis B — which states change with remission (exploratory, one primary comparison)
- For every author-annotated state with >= 20 cells in >= 20 paired patients: per-patient delta = Post fraction - Pre fraction (fraction within the state's parent compartment); two-sided Mann-Whitney on delta, Remission vs Non_Remission; AUROC of delta for remission. Benjamini-Hochberg across states within this analysis.
- B1 (primary comparison): the best state-level delta AUROC versus the pseudobulk INFLAM6 delta AUROC (same patients). Single-cell resolution is claimed to add change information ONLY if a state's delta AUROC exceeds the pseudobulk inflammation delta AUROC with the difference's bootstrap 95% CI excluding 0. Otherwise the conclusion is that, in paired data, resolution adds no discrimination beyond an aggregate inflammation change score in this cohort.
- The IAF delta result from the primary plan (AUROC 0.767) and the fibroblast INFLAM6 delta (0.785) are already known and are carried in as fixed reference values, not re-derived as findings.

## Analysis C — residual states in non-remitters (descriptive, pre-specified)
- In Non_Remission patients only: per-patient Post-treatment state fractions (states as in B). Report (i) the dominant residual state per patient (largest Post fraction among states whose Pre->Post delta >= 0), (ii) the distribution of dominant states across patients, (iii) pairwise patient similarity (cosine on Post state-fraction vectors) summarised as median within-group similarity.
- Interpretation fixed in advance: if a majority of non-remitters share one dominant residual state, the residual biology reads as common; if dominant states are dispersed across patients, non-remission is compositionally heterogeneous in a way a pseudobulk score cannot represent. No clustering, no new state discovery, no inferential test in C.

## Conduct rules
- No recalibration, no threshold changes, no gene substitutions. Missing universe genes handled as in the primary plan (available intersection).
- Every number computed under this plan is reported whichever way it lands, including a null B1.
- Lock = public GitHub push of this file + tag taurus-resolution-plan-v1 before any A1/B/C computation begins. Local commit timestamps are not evidence.
