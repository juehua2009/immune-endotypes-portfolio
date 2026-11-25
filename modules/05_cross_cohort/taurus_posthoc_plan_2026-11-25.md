# TAURUS / Martin post hoc analysis plan — 2026-11-25

**Status.** Post hoc. Written after the locked TAURUS primary analysis
(tag `taurus-plan-v1`), the resolution companion plan (tag
`taurus-resolution-plan-v1`) and the frozen-implementation re-scoring
(Stage 2, `frozen_aucell`, T_frozen = 0.108157) were all complete.
Locked by git commit and tag `taurus-posthoc-plan-v1` before any of the
analyses below are run. Nothing in this plan alters the locked primary
analysis or its interpretation; all three analyses are reported as
post hoc.

**Scoring.** All scores are the Stage 2 frozen scores (`frozen_aucell`:
fixed 16,632-gene human universe, n_up = 832 absolute, absent genes
zero-ranked last, universe-order tie-breaking; threshold
T_frozen = 0.108157). Patient unit and burden estimand as in the locked
plan (mean over a patient's samples of the per-sample share of
fibroblasts at or above T).

## Analysis 1 — TAURUS inflamed-only baseline sensitivity

Restrict to pre-treatment samples whose authors' annotation is
`Inflammation == "Inflamed"`. A patient's burden is the mean over their
inflamed pre-treatment samples only. Patients with no inflamed
pre-treatment sample are excluded; the reduced n (and the excluded
patients) is reported. On the remaining patients, report the AUROC of
pre-treatment burden against remission (positive class: non-remission)
with a bootstrap 95% CI (10,000 stratified resamples) and the one-sided
Mann–Whitney p in the locked direction (burden higher in
non-remitters). Labelled a post hoc sensitivity in all reporting.

## Analysis 2 — TAURUS biopsy-site description

Collapse the deposit's `Site` values to ileum vs colon using the
authors' `Ileum_vs_Colon` field. Report, per site, the number of
patients and samples contributing pre-treatment burdens and the median
pre-treatment burden by remission status, with bootstrap 95% CIs for
the medians. **Descriptive only: no hypothesis test and no conclusion
is drawn from this analysis.**

## Analysis 3 — Martin minimum-fibroblast sensitivity

Repeat the Martin paired involved-vs-uninvolved burden contrast
(per-patient shares at T_frozen, paired Wilcoxon) after excluding
samples with fewer than 20 fibroblasts — the same minimum-count rule
the locked TAURUS plan applies as a flag. Report how many samples and
donors the exclusion removes (donor_122's involved sample, 9
fibroblasts, is expected to be among them) and the resulting medians
and paired p alongside the full-set values.

## Reporting

Results are appended to the module's results memo under a heading
identifying them as post hoc under this plan, with this plan's tag
cited. Any manuscript use labels them post hoc.
