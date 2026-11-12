# Clinical and genetic characterization of primary ciliary dyskinesia in a Chinese cohort
Analysis release 1.1.0. The reference directory contains aggregate results. The input template contains headers only. The clinical registry and intermediate participant-level files are held under controlled access. To request de-identified study data, contact the Institutional Ethics Committee of the Children’s Hospital, Zhejiang University School of Medicine (chyb@zju.edu.cn; reference 2022-IRB-019-A2). Access is subject to review and the applicable data-use agreement.

## Public example using only published counts
After installing the Python dependencies, run `python comparisons.py --registry-only results` to reproduce the six S5 Fig results directly from the published aggregate counts. Compare the generated s5_figure.csv with reference/s5_figure.csv.

## Run locally with authorized data
Use Python 3.12 and install the pinned requirements. The study used Python 3.12.14. The subgroup R analyses used R 4.6.0 with logistf 1.26.1; R 4.4.1 was used for other study analyses. Run the commands in order from this directory. Keep the entire private and results directories on authorized storage; the R input CSV in results contains participant-level information.

```sh
python -m pip install -r requirements.txt
python prepare_inputs.py /authorized/path/confirmed.xlsx private
python figure5a.py private/a.json results
python figure5b.py private/b.json results
python subgroup_models.py private/a.json results
Rscript subgroup_firth.R results/private_subgroups.csv subgroup_pairs.csv results/R_subgroups.csv
python comparisons.py private/a.json results
```

| Paper item | Analysis and output |
|---|---|
| Fig 5A / S2 File | figure5a.py. All 431 listed gene–phenotype comparisons in model_pairs.csv; display_pairs.csv defines the selected rows. The output contains counts, unadjusted and adjusted profile penalized-likelihood intervals and likelihood-ratio p-values. BH adjustments across genes within a phenotype and across all estimable models are distinct. |
| Fig 5B / S3 File | figure5b.py. Unadjusted group-only and common-slope age-adjusted models; estimates at the cohort mean measurement age, 21 contrasts per model with Holm correction, partial F-tests, robust HC3 inference, age-by-group, sex-adjusted and dual-gene sensitivity analyses. S3_all_contrasts.csv includes both models. In the age-adjusted output p_unadjusted denotes a nominal p before multiplicity correction. |
| Fig 5C / S4 and S5 Files | subgroup_firth.R for inheritance-pattern analyses except maternal miscarriage; subgroup_models.py for all variant-type analyses and maternal miscarriage. S5 includes only genetically characterized individuals. Complete model rows can include groups too small for display. |
| Fig 6C and 6D | comparisons.py. Mutually exclusive gene or beating-pattern groups with at least three evaluable patients. Two-sided asymptotic Mann–Whitney tests with tie and continuity corrections; Holm correction over 15 or six contrasts, respectively. |
| Fig 6E / S6 File | subgroup_models.py for the 16 updated phenotype analyses and maternal miscarriage; subgroup_firth.R for infertility. Compare each beating pattern with the remaining eligible patterns for each clinical phenotype. BH correction within each phenotype. |
| Fig 6F | subgroup_models.py. For each beating pattern, test each eligible gene group against the other genetically characterized participants with pattern data. The dual-gene group contributes once to each comparator. |
| S5 Fig | comparisons.py. Two-sided Fisher exact tests for six comparisons between the Chinese cohort and the ERN-LUNG registry; conditional odds ratios and exact conditional confidence intervals. Registry counts are supplied in the script. |

## Input schema
Create a worksheet named Confirmed with the columns in input_template.csv. Its first column contains a stable coded participant identifier. Keep one row per patient–gene record; repeated participant identifiers represent the same individual and their clinical information must agree. Use the adjudicated patient-level variant category. Use blank or NA for unavailable values. The phenotype mapping in prepare_inputs.py identifies the source column for every binary endpoint; positive/negative values are yes/no. The historical field names are preserved for compatibility. Data checks distinguish 130 unique confirmed patients from the 125 with an identified PCD-associated gene. Age is continuous in years; enrollment age was the measurement age for spirometry in this cohort. EFV1% is the original field for percent-predicted FEV1. CBF is in hertz; CBP is categorical. Sex is coded M/F. Use numeric values at original precision.

The ascertainment of fetal lateral ventricular enlargement uses the original prenatal ultrasound report. Maternal miscarriage denotes reported spontaneous miscarriage and excludes ectopic pregnancy and induced abortion. Covariates and outcomes must be available for the relevant adjusted analysis. Unrepresented outcomes/groups and nonestimable results are not evidence of an odds ratio of zero. The adjusted heterozygosity–bronchiectasis result is blank in S4 because the original analysis identified competing solutions. For Fig 5A the dual-gene patient enters each relevant carrier group in separate models. For Fig 5B and Fig 6 the gene groups are mutually exclusive.

## Aggregate reference files
S2/S4/S5/S6_models.csv are tabular transcriptions of the corresponding final supplementary workbooks, preserving published numeric precision and unavailable values. S3_contrasts.csv and S3_tests.csv contain the two worksheets of S3 File. Additional figure-level reference outputs are listed in reference/. The release excludes original patient rows, personal identifiers, clinical reports, and generated private inputs. Underlying individual-level information is obtained through the access process above. The blank template documents the layout; it is not a substitute for the authorized cohort in analyses with study-specific sample size checks.

## Reproduction notes
All Fig 6E and Fig 6F models use enrollment age at its original recorded precision.

Independent reruns matched all 431 S2 models, 69 S4 models, 158 S5 models and 87 S6 models against the corresponding aggregate tables (relative numerical tolerance 0.0001, absolute tolerance 0.0000001). All 21 raw and age-adjusted FEV1 contrasts were checked. The R branch was additionally tested with R 4.5.3 and logistf 1.26.1. To install the study R dependency, use logistf version 1.26.1 in R.
