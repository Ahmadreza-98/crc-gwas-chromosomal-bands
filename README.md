# CRC GWAS chromosomal-band analysis

Repository accompanying the published article:

**The possibility of prognostic and functional values of the 8q24 and 20q13 chromosomal bands in colorectal cancer**

- **Author:** Seyed Ahmadreza Siadat
- **ORCID:** https://orcid.org/0009-0004-1439-6306
- **Journal:** *Molecular Biology Research Communications* (2026), 15(1):3–10
- **DOI:** https://doi.org/10.22099/mbrc.2025.54114.2202
- **PubMed:** https://pubmed.ncbi.nlm.nih.gov/41346766/
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12673626/
- **Language:** R
- **Status:** publication-matched default workflow with optional current-database refresh scripts

## Quick start

Open the project folder in R or RStudio and run:

```r
source("analysis.R")
```

No local path needs to be edited. The script locates the project directory, checks the required files, and writes generated figures and tables to `results/`.

The default run uses the publication-associated files retained in this repository: the final Supplementary workbook for GO/KEGG and cSurvival, the archived publication-era band-gene lists, and the archived mutation tables. This is intentional. Ensembl, GO, KEGG, MSigDB, GDC and R/Bioconductor packages change over time, so a current online rerun is not expected to reproduce every historical annotation or p-value exactly.

## Optional database refresh

The scripts in `scripts/` are for users who want to repeat database-dependent stages with current online resources:

```r
source("scripts/01_refresh_GO_KEGG.R")
source("scripts/02_refresh_TCGA_mutations.R")
source("scripts/03_prepare_cSurvival_candidates.R")
```

They write new files under `data/refresh/`. `analysis.R` checks that folder first. If a complete refreshed result exists for a stage, the refreshed file is used; otherwise the publication-associated archived file is used.

After running `01_refresh_GO_KEGG.R`, the next run of `analysis.R` uses refreshed Ensembl band-gene lists and the refreshed GO/KEGG workbook. After running `02_refresh_TCGA_mutations.R`, it uses the refreshed mutation tables. cSurvival itself is a web service, so the third script prepares the current reciprocal GO/KEGG candidates; a completed current cSurvival result can be saved as:

```text
data/refresh/csurvival/Table_2_survival.xlsx
```

and will then be used automatically.

To return to the publication-matched run, remove generated files under `data/refresh/` and run `analysis.R` again.

## Repository layout

```text
crc-gwas-chromosomal-bands/
├── analysis.R
├── README.md
├── DATA_PROVENANCE.md
├── CITATION.cff
├── LICENSE
├── environment/
│   ├── README.md
│   └── sessionInfo_tested.txt
├── scripts/
│   ├── 01_refresh_GO_KEGG.R
│   ├── 02_refresh_TCGA_mutations.R
│   ├── 03_prepare_cSurvival_candidates.R
│   └── README.md
├── data/
│   ├── raw/
│   │   ├── colorectal_cancer_GWAS_MONDO_0005575.xlsx
│   │   ├── cytoBand_hg38.txt
│   │   ├── COAD_UALCAN_top250_upregulated.txt
│   │   ├── COAD_UALCAN_top250_downregulated.txt
│   │   └── TCGA_COAD_READ_GDC_manifest.txt
│   ├── archived/
│   │   ├── band_genes/
│   │   ├── enrichment/
│   │   └── published_tables/
│   └── refresh/
└── archive/
    ├── README.md
    └── original_analysis_redacted.R
```

`results/` is created automatically and is ignored by Git.

## Main analysis

The main script:

1. reads the 1,796 CRC GWAS associations and removes duplicate SNP identifiers;
2. calculates major- and minor-band polymorphism distributions;
3. identifies the 13 significant major chromosomal bands;
4. creates the Manhattan and chromosomal distribution figures;
5. summarizes the coding/non-coding variant classes;
6. reads the archived or refreshed GO/KEGG enrichment results and creates enrichment figures;
7. repeats the DEG enrichment test using the archived or refreshed band-gene lists;
8. reads the archived or refreshed cSurvival and mutation results; and
9. writes `sessionInfo()` and `data_sources_used.tsv` for the run.

The 13 significant major bands are:

```text
10p14, 10q25, 11q12, 12p13, 15q13, 18q21, 19q13,
1q41, 20p12, 20q13, 6p21, 8q24, 9q34
```

The DEG enrichment analysis identifies **20q13** and **8q24** as significant.

## Notes on the published calculations

The chromosomal-band analysis preserves the fixed 3.2 Gb genome-length denominator and the multiplicative p-value correction used in the original working script and published band table. In the cleaned code this denominator is named `PUBLISHED_GENOME_LENGTH_BP` so that the methodological choice is explicit rather than hidden as a numeric literal.

Among the 1,346 unique SNPs, 637 map to at least one protein-coding gene. The variant-class figure contains 632 coding-consequence variants because five variants map to protein-coding genes but have downstream-only consequence annotations.

### Reproducibility note on Table 3

For the DEG enrichment test, the original working R code uses Benjamini-Hochberg adjustment. The adjusted values printed in the published Table 3 correspond to the raw p-values multiplied by 13. Both approaches identify the same two significant bands, **8q24** and **20q13**. `analysis.R` therefore reports both values rather than silently replacing either one.

## cSurvival

Survival analysis was performed with the cSurvival web tool, not calculated locally in R. For the publication-matched run, `analysis.R` reads Table S5 directly from the final Supplementary workbook. The optional cSurvival script only prepares current reciprocal GO/KEGG candidates; the web analysis still has to be run separately.

## R environment

A successful validation run used **R 4.4.1**. The exact session information from that run is retained in `environment/sessionInfo_tested.txt`; a concise dependency summary and reproducibility note are in `environment/README.md`.

A hand-written `renv.lock` is deliberately not included. The tested machine contained several Bioconductor development-snapshot package versions, so reconstructing a lockfile by hand would imply a level of restore fidelity that has not been verified. A true lockfile should be generated with `renv::snapshot()` in the exact tested R library. This is documented in `environment/README.md`.

## Historical working script

The historical working script is retained at:

```text
archive/original_analysis_redacted.R
```

Its analysis logic is preserved, but the original machine-specific `setwd()` path has been replaced by a placeholder before public release. The file is provided for provenance and is **not** the recommended entry point. See `archive/README.md`.

## Citation

If you use this repository, please cite the associated article:

> Siadat SA. The possibility of prognostic and functional values of the 8q24 and 20q13 chromosomal bands in colorectal cancer. *Molecular Biology Research Communications*. 2026;15(1):3–10. doi:10.22099/mbrc.2025.54114.2202. PMID: 41346766; PMCID: PMC12673626.

A machine-readable citation is also provided in `CITATION.cff` so GitHub can display **Cite this repository**.

## License

The repository code and documentation are released under the MIT License; see `LICENSE`. Third-party datasets and database-derived files remain subject to the terms of their original sources. The published article and its Supplementary material are available under the license stated by the journal/PMC.

For file-level provenance and source limitations, see `DATA_PROVENANCE.md`.
