# PCD cohort analysis code for Figures 5A and 5B

Version 1.0.0, 2026-11-03. Companion to the existing R scripts. This release contains the Python estimator and analysis code used in revision, with portable input/output paths. Python 3.12.14 and the pinned packages in requirements.txt were used for verification. The PDF output is a statistical plot. Final manuscript artwork was arranged in Adobe Illustrator.

## Run

From this directory, in a Python 3.12 virtual environment:

    python -m pip install -r requirements.txt
    python prepare_inputs.py /path/to/confirmed-workbook.xlsx private_inputs
    python figure5a.py private_inputs/a.json results
    python figure5b.py private_inputs/b.json results

Access to the controlled clinical data is required for participant-level reanalysis. Keep the generated private_inputs directory outside public repositories. No participant-level records are included in this release. Both analysis scripts also accept equivalently formatted authorized JSON inputs.

## Preparation and units

prepare_inputs.py maps the frozen Confirmed worksheet headers, including the source spelling EFV1%. Enrollment age is used as the confirmed age at spirometry for these measurements. Yes and No become 1 and 0; undocumented phenotypes remain missing. Rows 52 and 53 refer to one individual with DNAH5 and DNAH11 findings. Their covariates and phenotypes are checked before merging. Run this row-specific importer only on the original layout.

Input A is a JSON array with one object per confirmed participant: patient_id (unique coded string), genes (list of gene symbols, empty for genetically unresolved cases), enrollment_age (years), sex (M/F), and phenotypes (mapping from phenotype names in model_pairs.csv to 1, 0 or null). Input B contains unique patients with paired spirometry and age: patient_id, gene, sex, age_at_spirometry_years, FEV1_percent_predicted (percentage points), multigene_patient (JSON boolean). For the dual-gene record use a combined gene label distinct from either single gene. Four columns from the primary sample alone cannot reproduce the sex and dual-gene sensitivity analyses.

## Figure 5A and S3 File

Model pairs and plotted rows are specified in separate CSV files. Each model compares carriers of the gene to other genetically characterized patients using available cases. The dual-gene individual belongs to each relevant carrier group in its respective model and appears only once per model. Nonvarying sex is omitted. The primary models adjust for continuous enrollment age and sex. Unadjusted fits are also reported. firth.py estimates the gene odds ratio by penalized likelihood, its 95% profile penalized-likelihood interval, and a two-sided penalized likelihood-ratio p-value. Age is centered and scaled for numerical conditioning. The all-tests and within-phenotype Benjamini–Hochberg corrections are separate outputs.

figure5a_results.json contains four-cell counts (carrier yes/no, comparator yes/no), sample sizes, estimates and multiplicity corrections for all 431 gene–phenotype pairs, including the nonestimable pair. figure5a_display.csv contains the 12 displayed comparisons, in figure order. Selection reproduces the manuscript panel. These 12 rows are a subset of the estimable associations.

## Figure 5B and S4 File

Primary model: FEV1 percent predicted as a function of gene group and continuous age at spirometry, with a common age slope. The partial F-test compares this model to an age-only model. The primary analysis excludes the dual-gene participant. Group means and 95% confidence intervals are estimated at the analysis sample's mean age. The 21 contrasts in S4 are differences between age-adjusted means; their confidence intervals are pointwise. p_unadjusted is the age-adjusted contrast p-value before correction for multiplicity; p_Holm is corrected over all 21 contrasts.

Outputs include group means (raw and age-adjusted), the full pairwise CSV, a vector PDF, and figure5b_results.json. The latter includes the primary model, separate reassignment of the dual-gene participant to DNAH5 or DNAH11, further adjustment for sex, and the HC3 robust omnibus Wald F-test. Residual normality and maximum leverage are printed to the console. Contrast direction follows the S4 File order; the display order of HYDIN and RSPH4A follows the final figure.
