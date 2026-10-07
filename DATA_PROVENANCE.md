# Data provenance

This file distinguishes starting inputs, publication-associated archived outputs, and files that can be regenerated with current online resources. When an exact historical retrieval date or database release was not preserved in the retained project files, that fact is stated explicitly rather than inferred.

## Published article

**Siadat SA.** The possibility of prognostic and functional values of the 8q24 and 20q13 chromosomal bands in colorectal cancer. *Molecular Biology Research Communications*. 2026;15(1):3–10. DOI: `10.22099/mbrc.2025.54114.2202`. PMID: `41346766`. PMCID: `PMC12673626`.

## Starting inputs

| File | Source | Retrieval/build information retained in the project | Use |
| --- | --- | --- | --- |
| `data/raw/colorectal_cancer_GWAS_MONDO_0005575.xlsx` | NHGRI-EBI GWAS Catalog, colorectal cancer trait `MONDO_0005575` | The article records retrieval on **2024-07-12** from a dataset representing 75 CRC GWAS studies; the retained workbook contains 1,796 association rows | Starting GWAS associations and annotations |
| `data/raw/cytoBand_hg38.txt` | hg38/GRCh38 cytoband coordinates, retained from the original analysis | Genome build is **GRCh38/hg38**; the exact historical download date and source release were not recorded in the retained file | Major/minor band coordinates and lengths |
| `data/raw/COAD_UALCAN_top250_upregulated.txt` | UALCAN, TCGA-COAD | Top 250 upregulated genes retained from the publication-era analysis; exact retrieval date was not recorded | DEG enrichment |
| `data/raw/COAD_UALCAN_top250_downregulated.txt` | UALCAN, TCGA-COAD | Top 250 downregulated genes retained from the publication-era analysis; exact retrieval date was not recorded | DEG enrichment |
| `data/raw/TCGA_COAD_READ_GDC_manifest.txt` | NCI Genomic Data Commons | Manifest retained from the TCGA-COAD/TCGA-READ masked somatic mutation download; exact query date was not recorded in the retained project files | Provenance for the mutation-data download |

The article reports 1,796 CRC-associated polymorphisms from 75 GWAS studies. After duplicate SNP identifiers were removed, 1,346 unique variants remained.

## Archived analysis files

### Band gene lists

`data/archived/band_genes/` contains publication-era protein-coding gene symbols for the 13 significant major bands. These lists were obtained through Ensembl/`biomaRt` during the original analysis. The exact Ensembl release was not recorded in the retained files, so the repository does not assign a release retrospectively.

These archived lists are used by default for the DEG calculation. `scripts/01_refresh_GO_KEGG.R` can rebuild them from the Ensembl service available at the time of rerun. Complete refreshed lists are written to `data/refresh/band_genes/` and take precedence on the next run of `analysis.R`.

### GO/KEGG enrichment

`data/archived/enrichment/Supplementary_Tables_published.xlsx` is the final Supplementary workbook distributed with the published article. `analysis.R` reads Tables S1-S4 from this workbook for the publication-matched GO, KEGG and reciprocal-enrichment results, and Table S5 for the publication-matched cSurvival results.

The original workflow used `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `KEGGREST`, Ensembl/`biomaRt`, and the MSigDB C1 positional collection through `msigdbr`. Because those resources are versioned or updated over time, a current rerun can legitimately differ from the archived publication-era results.

`scripts/01_refresh_GO_KEGG.R` repeats the database-dependent enrichment stages with current resources and writes:

```text
data/refresh/enrichment/Supplementary_Tables_current.xlsx
```

If that workbook exists, `analysis.R` uses it instead of the archived enrichment workbook.

### cSurvival

Survival analysis was performed outside R with the cSurvival web tool using TCGA-COAD/TCGA-READ data. The final publication-associated survival results are Table S5 of:

```text
data/archived/enrichment/Supplementary_Tables_published.xlsx
```

No separate duplicate survival workbook is kept in the archived tree.

`scripts/03_prepare_cSurvival_candidates.R` prepares current reciprocal GO/KEGG candidates for a new web-based cSurvival run. A completed refreshed survival table can be saved as:

```text
data/refresh/csurvival/Table_2_survival.xlsx
```

and will then be used by `analysis.R`.

### Mutation analysis

The original mutation analysis queried open-access masked somatic mutation data for TCGA-COAD and TCGA-READ through GDC/TCGAbiolinks, summarized mutations with `maftools`, and used Ensembl transcript CDS lengths to calculate adjusted mutation burden.

The publication-associated result tables retained in the repository are:

```text
data/archived/published_tables/Table_4_mutation_enrichment.xlsx
data/archived/published_tables/Table_5_mutation_ttest.xlsx
```

`scripts/02_refresh_TCGA_mutations.R` repeats the mutation workflow with the current GDC and Ensembl services. It writes replacement tables under `data/refresh/mutation/`. The main analysis switches to the refreshed mutation results only when both required tables are present.

### Other archived tables

`data/archived/published_tables/Table_1_band_statistics.xlsx` is retained as the publication-associated band table, although the main script recalculates band statistics directly from the GWAS and cytoband inputs.

`data/archived/published_tables/Table_3_DEG_enrichment.xlsx` is retained as the published DEG table. The main script also recalculates DEG enrichment locally and exposes both the BH adjustment used in the working R code and the `p × 13` values shown in the published Table 3.

## Fixed publication parameter

The chromosomal-band test in the original working script used a genome-length denominator of **3.2 × 10^9 bp**. The cleaned workflow preserves that value as `PUBLISHED_GENOME_LENGTH_BP <- 3.2e9` for publication reproduction rather than silently substituting a newly calculated genome size.

## Online resources

- GWAS Catalog trait page: <https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005575>
- UCSC Genome Browser: <https://genome.ucsc.edu/>
- Ensembl: <https://www.ensembl.org/>
- Gene Ontology: <https://geneontology.org/>
- KEGG: <https://www.genome.jp/kegg/>
- MSigDB: <https://www.gsea-msigdb.org/gsea/msigdb/>
- UALCAN: <https://ualcan.path.uab.edu/>
- NCI Genomic Data Commons: <https://portal.gdc.cancer.gov/>
- cSurvival publication: <https://doi.org/10.1093/bib/bbac090>

## Reproducibility note

Ensembl annotations, GO/KEGG membership, MSigDB collections, GDC files and package versions change over time. The repository therefore uses publication-associated archived files by default and keeps current-database reruns separate under `data/refresh/`. A refreshed run is intended to answer “what does the same workflow return using current resources?”, not to guarantee byte-for-byte recreation of historical database outputs.

Exact package versions from a successful validation run are retained in `environment/sessionInfo_tested.txt`. See `environment/README.md` for the reason a synthetic hand-written `renv.lock` is not included.
