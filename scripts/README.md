# Optional reruns

The main `analysis.R` uses the final published database-dependent results by default. The scripts in this folder are for users who want to repeat those stages with current online resources.

Run them from the project root.

## `01_refresh_GO_KEGG.R`

Rebuilds the protein-coding gene lists for the significant chromosomal bands from current Ensembl, then repeats GO, KEGG and reciprocal chromosomal enrichment.

Each GO term is processed independently in the reciprocal GO loop. The output workbook follows the final Supplementary numbering:

```text
Table S1  GO
Table S2  KEGG
Table S3  reciprocal GO
Table S4  reciprocal KEGG
```

The new gene lists and workbook are written under `data/refresh/`. The next run of `analysis.R` uses them automatically.

## `02_refresh_TCGA_mutations.R`

Downloads current open-access TCGA-COAD and TCGA-READ masked somatic mutation data from GDC, retrieves current transcript CDS lengths from Ensembl and repeats the mutation-burden tests used in the original script.

If refreshed band gene lists from script 01 are available, they are used; otherwise the archived band gene lists are used. The mutation tables are written under `data/refresh/mutation/`.

## `03_prepare_cSurvival_candidates.R`

cSurvival is a web service, so the survival model itself is not rerun in R. This script collects the current reciprocal GO/KEGG candidates into a workbook for a new cSurvival run.

After the web analysis, save the final result as:

```text
data/refresh/csurvival/Table_2_survival.xlsx
```

The next run of `analysis.R` will use it instead of the published survival table.

## Returning to the published run

Delete the generated files under `data/refresh/`. The next run of `analysis.R` will fall back to the final published files under `data/archived/`.
