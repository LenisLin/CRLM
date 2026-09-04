# Figure 6

`02_cohort_and_rfs.R` generates the Figure 6B clinical-distribution components
and Figure 6C-D recurrence-free survival panels from prepared patient tables.
`04_continuous_isr_multivariable_cox.R` is the single Figure 6E-F
continuous-ISR multivariable Cox entry point. It applies the same model
specification separately to the two cohorts using the Supplementary Tables
workbook and the FDZS-1 patient-level H&E-derived ISR table.
`03_prepare_public_isr_tables.py` generates the public FDZS-1 and FDZS-2 ISR/RFS
tables, discovery/test assignment, and score-definition/cutoff table.
The R analysis uses host Conda environment `Spatial` with R 4.2.2, and the
Python table export uses `Spatial_py`; package versions are listed in
`envs/package_versions.yml`.

Edit the top-level input and output paths in the relevant script, enable its run
switch, and run from the repository root without arguments:

```text
Rscript figures/Figure6/02_cohort_and_rfs.R
python figures/Figure6/03_prepare_public_isr_tables.py
Rscript figures/Figure6/04_continuous_isr_multivariable_cox.R
```

## Required prepared tables

Provide discovery and test TSVs with these columns:

| Table | Required columns | Rows and analysis unit |
|---|---|---|
| Discovery RFS | `patient_id`, `wsi_id`, `isr`, `rfs_time_months`, `rfs_event`, `treatment` | 34 patient rows |
| Test RFS | Discovery columns plus `TBS`, `CRLM_number`, `CRLM_size`, `fong_score`, `KRAS_mutation` | 95 patient rows |

`patient_id` occurs once per patient row. The analysis uses the patient-to-WSI
and ISR mapping recorded in each table; `wsi_id` can be shared across patient
rows. `rfs_event` uses conventional coding: `1=event` and
`0=censored`. Treatment values are `Chemo` and `Combo`.

## Analysis and outputs

The discovery-table median ISR is the cutoff: `isr > discovery median` defines
High ISR and values equal to the median define Low ISR. The same numeric cutoff
is applied to the test table. Figure 6C-D display the four treatment-by-ISR
strata with Kaplan-Meier curves, confidence bands, risk tables, and nominal
two-sided log-rank comparisons. Panel annotations present, in order: Combo
High versus Chemo High, Combo High versus Combo Low, Chemo High versus Chemo
Low, and Combo Low versus Chemo Low.

The entry point writes:

- `Figure6B_test_clinical_distributions.pdf`
- `Figure6B_distribution_summary.tsv`
- `Figure6C_discovery_rfs.pdf`
- `Figure6D_test_rfs.pdf`
- `Figure6C_group_counts.tsv`
- `Figure6D_group_counts.tsv`
- `Figure6C_risk_counts.tsv`
- `Figure6D_risk_counts.tsv`
- `Figure6C_pairwise_logrank.tsv`
- `Figure6D_pairwise_logrank.tsv`
- `Figure6_isr_cutoff.tsv`

The entry point reads the editable in-script configuration and consumes
patient-level ISR tables supplied by the SpMap workflow.

## Continuous-ISR multivariable Cox model

The single Figure 6E-F entry point reads the 35 FDZS-1 clinical rows in Supplementary
Table 1, merges the 34 rows marked `figure6_included == 1` in
`FDZS1_patient_level_ISR.tsv`, and reads all 95 FDZS-2 rows directly from
Supplementary Table 5. It fits the same model separately in each cohort:

```text
Surv(rfs_time_months, rfs_event) ~ isr * treatment + age + gender + kras +
  fong + cea + ca199 + tbs + crlm_number + crlm_size +
  differentiation_grade + t_stage + lymph_positive
```

ISR is retained as a continuous score, and `Chemo` is the treatment reference.
Age, TBS, CRLM number, largest CRLM size, CEA, and CA19-9 are standardized
within cohort before fitting. Fong score and pathological T stage remain
unstandardized numeric covariates; gender is a two-level cohort-specific factor,
and the remaining binary fields retain their displayed coding. Each model has
15 coefficients: the ISR and treatment main effects, their interaction, and 12
clinicopathological covariates.

The entry point writes:

- `Figure6EF_continuous_ISR_multivariable_Cox.tsv`
- `Figure6E_continuous_ISR_multivariable_Cox.pdf`
- `Figure6F_continuous_ISR_multivariable_Cox.pdf`

## Public patient-level outputs

Edit the workbook, saved IMC score, discovery H&E score, and output paths in
`03_prepare_public_isr_tables.py`. The script writes:

- `FDZS1_patient_level_ISR.tsv`: 35 FDZS-1 patients with IMC-derived ISR and
  available H&E-derived ISR; 34 rows are marked for the Figure 6 discovery set.
- `FDZS2_patient_level_ISR_RFS.tsv`: 95 pseudonymized independent-test patient
  rows matching the Figure 6 input schema.
- `SpMap_patient_cohort_assignment.tsv`: 34 discovery and 95 independent-test
  patient assignments.
- `ISR_definition_and_cutoff.tsv`: operational formulas and the discovery
  median H&E-derived ISR cutoff.

For the discovery input to `02_cohort_and_rfs.R`, select the 34 rows with
`figure6_included == 1` from `FDZS1_patient_level_ISR.tsv` and retain the
required analysis columns listed above.
