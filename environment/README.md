# Tested R environment

`sessionInfo_tested.txt` records the R session from a successful publication-matched validation run.

Key versions used in that session included:

| Component | Version |
| --- | --- |
| R | 4.4.1 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| ggplot2 | 3.5.1 |
| openxlsx | 4.2.7.1 |
| patchwork | 1.3.0 |
| scales | 1.3.0 |
| chromPlot | 1.33.0 |
| biomaRt | 2.61.3 |
| clusterProfiler | 4.14.4 |
| org.Hs.eg.db | 3.20.0 |
| AnnotationDbi | 1.67.0 |
| KEGGREST | 1.45.1 |
| msigdbr | 7.5.1 |
| TCGAbiolinks | 2.33.0 |
| maftools | 2.20.0 |
| BiocManager | 1.30.25 |

## Why there is no hand-written `renv.lock`

Using `renv` is a good way to make an R environment easier to recreate, but a lockfile should be generated from the actual tested project library with `renv::snapshot()`. The validation session contained several Bioconductor development-snapshot versions (for example `biomaRt 2.61.3` and `AnnotationDbi 1.67.0`). Reconstructing a lockfile manually from `sessionInfo()` would risk producing a file that looks authoritative but cannot be restored exactly from today’s repositories.

For that reason this repository retains the complete tested `sessionInfo()` rather than shipping an unverified synthetic lockfile.

To create a true lockfile on the machine/environment used for a future validated run:

```r
install.packages("renv")
renv::init(bare = TRUE)
renv::snapshot()
```

Commit the resulting `renv.lock` only after `renv::restore()` has been tested in a clean R installation.

The publication-matched default analysis is additionally protected from database drift by using the archived publication-associated GO/KEGG, cSurvival, band-gene and mutation outputs. The optional refresh scripts intentionally query current online resources and can therefore return updated results even when package versions are pinned.
