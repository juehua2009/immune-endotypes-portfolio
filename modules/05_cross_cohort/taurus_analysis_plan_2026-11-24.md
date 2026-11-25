# TAURUS validation — locked analysis plan

**Locked: 2026-11-24 (EST), before any TAURUS scoring output was examined.**
This plan is committed to the repository at lock time; the TAURUS scoring run
(step7) was launched only after this commit.

## Cohort and treatment
Thomas et al. 2024 (TAURUS), fibroblast/pericyte object `fibperi_final.h5ad`
(Zenodo 10.5281/zenodo.14007626, v3 release of 2024-10-30). Biologic-naïve IBD
patients sampled before and after **adalimumab (anti-TNF)** — the same drug
class as the bulk anti-TNF cohorts, so the directional hypothesis below is a
direct single-cell counterpart of the bulk result, not a new one.

## Frozen ruler (no parameter is re-estimated in TAURUS)
- Gene set: SIG8 = IL11, IL13RA2, CHI3L1, MMP3, MMP10, IL6, CXCL1, CXCL5.
- Scoring: AUCell (decoupler, tmin = 1) on the frozen gene universe — the
  16,632-gene intersection of the Smillie / Martin / Kinchen matrices
  (calibration universe; TAURUS matrices are restricted to it before scoring).
- Threshold: T = 0.107229, calibrated once in Smillie 2019 against its
  published inflammatory-fibroblast annotations. No recalibration.
- Denominator: fibroblasts = the TAURUS authors' `major` annotation ≠ Pericyte.
- Statistical unit: the **patient**. Pre-treatment burden = mean of per-sample
  IAF shares (fibroblasts ≥ T) over the patient's Pre samples.

## Primary test (one test, direction fixed in advance)
Pre-treatment IAF burden vs remission status — one-sided Mann–Whitney U,
direction: **burden higher in Non_Remission** (fixed by the bulk result: a high
signature predicts anti-TNF non-response).

Analysis set: patients with ≥1 Pre sample and Remission_status ∈ {Remission,
Non_Remission} — **n = 36** (Remission 16, Non_Remission 20; CD 10 R / 6 NR;
UC 6 R / 14 NR). Power note: at one-sided α = 0.05, n = 16 vs 20 gives 80%
power only for AUC ≳ 0.74; a smaller true effect will read as inconclusive,
not negative.

## Pre-specified sensitivity (confounding by baseline inflammation)
Logistic regression Remission ~ burden (z) + baseline inflammation, where
baseline inflammation is the patient's mean sample Inflammation_score over
Pre samples. A second sensitivity replaces it with the per-sample single-cell
INFLAM6 AUCell mean over fibroblasts. The adjusted odds ratio for burden is
reported alongside the primary test.

## Secondary analyses (exploratory, interpretation fixed in advance)
Paired Pre → Post change in burden (Wilcoxon), overall and by remission.
A post-treatment decline in remitters is expected as inflammation resolves
and is **not** evidence of prediction; only the primary test speaks to
prediction. Shares by inflammation stratum and the Healthy arm (n = 3) are
descriptive only.

## Exclusions
Remission_status Not_avail / None: excluded from the primary test.
Pericytes: excluded from every denominator.

## Amendment 1 — 2026-11-24 (same lock session, still before any scoring)

Added after review of the cohort's design table (patient counts by disease and
remission only; no burden value had been computed or seen):

1. **Disease-stratified sensitivity.** The outcome groups are imbalanced by
   disease (Non_Remission is UC-heavy, 14/20; Remission is CD-heavy, 10/16),
   and sampling site (ileum vs colon) differs by disease. The primary test is
   therefore accompanied by (a) a van Elteren stratified Wilcoxon test of
   pre-treatment burden vs remission, stratified by Disease, and (b) the
   adjusted logistic model of the sensitivity section with Disease added as
   a covariate. If the unstratified and stratified results disagree, the
   stratified result governs the interpretation.
2. **Minimum fibroblast count.** A patient's share is only as stable as its
   denominator. Patients whose Pre samples contain fewer than 20 fibroblasts
   in total are flagged in every table, and the primary test is repeated
   excluding them as a named sensitivity.
3. **Gene coverage note.** The TAURUS annotation carries 14,537 of the 16,632
   frozen-universe genes (2,095 absent, predominantly lowly expressed
   non-coding loci); all 8 SIG8 genes are present. Scoring restricts TAURUS
   to the available intersection — the frozen assay applied to the genes the
   deposited matrix measures. No gene is substituted.

## Reporting commitment
The primary result is reported whichever way it lands. The manuscript will
cite this plan's commit hash and state the public (GitHub push) timestamp
honestly, including the order in which the scoring run and the push occurred.
