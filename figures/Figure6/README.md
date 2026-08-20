# Figure 6

`02_cohort_and_rfs.R` generates the Figure 6B clinical-distribution components
and Figure 6C-D recurrence-free survival panels from prepared patient tables.
It uses host Conda environment `Spatial` with R 4.2.2; package versions are
listed under `imc` in `envs/package_versions.yml`.

Edit `FIGURE6_CONFIG$discovery_rfs`, `test_rfs`, and `output_dir`, set
`RUN_FIGURE6` to `TRUE`, and run from the repository root without arguments:

```text
Rscript figures/Figure6/02_cohort_and_rfs.R
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
