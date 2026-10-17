# Figure reports

These reports organize the R and Python analysis entry points by figure and display the
current publication figures alongside the related analysis sections. Complete
figure layouts retain the panel labels, legends, and imaging examples. Captions
identify the displayed observations and statistical units.

| Figure | Analyses | Publication figures displayed | Read on GitHub | R Markdown source |
|---|---|---|---|---|
| 1 | Study framework and treatment-stratified RFS | Figure 1A–D; Figures S1–S2 | [Figure 1](../figures/Figure1/Figure1.md) | [Figure1.Rmd](../figures/Figure1/Figure1.Rmd) |
| 2 | Cell phenotypes, composition, and clinical associations | Figure 2A–D; Figures S3–S5 | [Figure 2](../figures/Figure2/Figure2.md) | [Figure2.Rmd](../figures/Figure2/Figure2.Rmd) |
| 3 | Cellular neighborhoods and spatial interactions | Figure 3A–F; Figures S6–S8 | [Figure 3](../figures/Figure3/Figure3.md) | [Figure3.Rmd](../figures/Figure3/Figure3.Rmd) |
| 4 | Peritumoral niches, epithelial features, and ISR | Figure 4A–J; Figures S9–S10 | [Figure 4](../figures/Figure4/Figure4.md) | [Figure4.Rmd](../figures/Figure4/Figure4.Rmd) |
| 5 | SpMap workflow, evaluation, and WSI predictions | Figure 5A–C; Figures S11–S12 | [Figure 5](../figures/Figure5/Figure5.md) | [Figure5.Rmd](../figures/Figure5/Figure5.Rmd) |
| 6 | Clinical distributions, survival, and continuous/binary ISR Cox models | Figure 6A–F; Figures S13–S14 | [Figure 6](../figures/Figure6/Figure6.md) | [Figure6.Rmd](../figures/Figure6/Figure6.Rmd) |
| 7 | Spatial composition, differential expression, pathways, and communication | Figure 7A–H; Figures S15–S18 | [Figure 7](../figures/Figure7/Figure7.md) | [Figure7.Rmd](../figures/Figure7/Figure7.Rmd) |

## Python notebooks

The notebooks present all nine Figure-owned Python entry points in workflow
sections, with current publication figures embedded directly in the files.

| Notebook | Python workflow | Embedded figures |
|---|---|---|
| [Figure 5](../figures/Figure5/Figure5_python.ipynb) | SpMap workflow, reference counts, public artifact preparation, OOF metrics, and ISR correlation | Figure 5; Figures S11–S12 |
| [Figure 6](../figures/Figure6/Figure6_python.ipynb) | Patient-level ISR and survival-table preparation | Figure 6 |
| [Figure 7](../figures/Figure7/Figure7_python.ipynb) | scRNA-seq integration, annotation, and processed-object plotting | Figure 7; Figure S15 |

Open the `.ipynb` files in GitHub or Jupyter to read the code and saved results.
Code cells are stored unexecuted and tagged `skip-execution`; saved images are
Markdown attachments. Script execution instructions are in each Figure README.

## Reading the reports

Each report introduces its inputs and analysis units, presents the figure-owned
R code in logical sections, and displays the publication results. The owning
Figure README provides execution instructions and input schemas. Shared
functions and companion Python notebooks are linked from this index.

Figure 6 provides [continuous and binary ISR models for both cohorts](../figures/Figure6/05_continuous_and_binary_isr_cox.R), alongside the existing
[continuous ISR entry](../figures/Figure6/04_continuous_isr_multivariable_cox.R).

Figures use relative paths into each Figure directory's `report_assets/` folder.
PDF panels have PNG previews and links to the original vector files. Processed
data are described in the [public data contract](data_contract.md) and available
through [Zenodo](https://doi.org/10.5281/zenodo.22300873).

## Rendering

The R Markdown documents use `github_document` with Markdown image links and
a single caption per image. All analysis chunks have
`eval=FALSE`; rendering assembles the displayed code, prose, tables, and stored
images. The rendering environment uses R, knitr, rmarkdown, and Pandoc.

From the repository root, render one report with:

```sh
Rscript -e 'rmarkdown::render("figures/Figure1/Figure1.Rmd")'
```

Change the figure number to render another report. The resulting `.md` file is
saved beside its `.Rmd` source.

## Maintaining the reports

The `.R` entry points are the source for the displayed code. Source comments
above the R Markdown chunks identify the script and line range. Update those
chunks and their ranges when the owning script changes, then render the reading
version again. Keep each panel paired with its generating analysis and retain
its figure number, title, and analysis unit.

The `.py` entry points supply the notebook code. Each code cell records
`source_file`, `source_start_line`, and `source_end_line` in its metadata.
Update the displayed code and ranges together when a source script changes.
Replace embedded image attachments when the publication figure is updated.
