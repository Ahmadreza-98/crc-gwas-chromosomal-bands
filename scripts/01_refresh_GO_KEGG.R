# Refresh the Ensembl band gene lists and repeat the GO/KEGG enrichment.
#
# This is optional. The main analysis uses the archived results from the paper
# unless this script has been run successfully. When refreshed files are present
# under data/refresh/, analysis.R will use them instead.

cran_packages <- c("dplyr", "tidyr", "openxlsx", "msigdbr")
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "biomaRt", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "KEGGREST"
)
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

find_project_dir <- function() {
  candidates <- unique(c(getwd(), dirname(getwd())))
  marker <- file.path("data", "raw", "colorectal_cancer_GWAS_MONDO_0005575.xlsx")
  hits <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hits) == 0) {
    stop("Run this script from the project root.", call. = FALSE)
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

PROJECT_DIR <- find_project_dir()
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
REFRESH_DIR <- file.path(PROJECT_DIR, "data", "refresh")
BAND_GENE_DIR <- file.path(REFRESH_DIR, "band_genes")
ENRICHMENT_DIR <- file.path(REFRESH_DIR, "enrichment")

dir.create(BAND_GENE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICHMENT_DIR, recursive = TRUE, showWarnings = FALSE)

# Recreate the significant major bands from the same local GWAS/cytoband inputs
# and the same band test used by the main analysis. The fixed 3.2 Gb denominator
# is preserved because it is part of the published calculation.
PUBLISHED_GENOME_LENGTH_BP <- 3.2e9
gwas <- openxlsx::read.xlsx(
  file.path(RAW_DIR, "colorectal_cancer_GWAS_MONDO_0005575.xlsx")
)
names(gwas) <- make.names(names(gwas), unique = TRUE)
gwas <- gwas[!duplicated(gwas$SNPS), , drop = FALSE]
gwas$major.band <- sub("\\..*$", "", gwas$REGION)

cytoband <- read.delim(
  file.path(RAW_DIR, "cytoBand_hg38.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cytoband <- cytoband[cytoband$Band != "", , drop = FALSE]
cytoband$major.band <- paste0(
  sub("^chr", "", cytoband$Chr),
  sub("\\..*$", "", cytoband$Band)
)

major_counts <- as.data.frame(table(gwas$major.band), stringsAsFactors = FALSE)
colnames(major_counts) <- c("major.band", "ni")
major_counts$ni <- as.numeric(major_counts$ni)

major_lengths <- cytoband |>
  group_by(major.band) |>
  summarise(Start = min(Start), End = max(End), .groups = "drop") |>
  mutate(Length = End - Start)

major_stats <- merge(major_counts, major_lengths, by = "major.band")
major_stats$n <- nrow(gwas)
major_stats$theta <- major_stats$Length / PUBLISHED_GENOME_LENGTH_BP
major_stats$df1 <- 2 * (major_stats$n - major_stats$ni + 1)
major_stats$df2 <- 2 * major_stats$ni
major_stats$F.test <- ((1 - major_stats$theta) / major_stats$theta) *
  (major_stats$ni / (major_stats$n - major_stats$ni + 1))
major_stats$pvalue <- pf(
  major_stats$F.test,
  major_stats$df1,
  major_stats$df2,
  lower.tail = FALSE
)
major_stats$adj.pvalue <- major_stats$pvalue * nrow(major_stats)

significant_bands <- major_stats[major_stats$adj.pvalue < 0.05, , drop = FALSE]
significant_bands$Chr <- sub("[a-z].*$", "", significant_bands$major.band)

message("Retrieving current Ensembl protein-coding genes for ", nrow(significant_bands), " bands...")
mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

band_genes <- list()
for (i in seq_len(nrow(significant_bands))) {
  band <- significant_bands$major.band[i]
  genes <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = c("chromosome_name", "start", "end", "biotype"),
    values = list(
      significant_bands$Chr[i],
      significant_bands$Start[i],
      significant_bands$End[i],
      "protein_coding"
    ),
    mart = mart
  )

  symbols <- unique(sub("-.*", "", genes$hgnc_symbol[genes$hgnc_symbol != ""]))
  band_genes[[band]] <- symbols
  write.table(
    symbols,
    file.path(BAND_GENE_DIR, paste0(band, "_genes.txt")),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

safe_go <- function(entrez, ontology) {
  if (length(entrez) == 0) return(data.frame())
  as.data.frame(
    clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  )
}

safe_kegg <- function(entrez) {
  if (length(entrez) == 0) return(data.frame())
  as.data.frame(
    clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )
  )
}

go_all <- list()
kegg_all <- list()

for (band in names(band_genes)) {
  mapped <- suppressMessages(
    clusterProfiler::bitr(
      band_genes[[band]],
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  entrez <- unique(mapped$ENTREZID)

  go_parts <- lapply(c("BP", "MF", "CC"), function(ontology) {
    x <- safe_go(entrez, ontology)
    if (nrow(x) > 0) {
      x$GO <- ontology
      x$major.band <- band
    }
    x
  })
  go_all[[band]] <- bind_rows(go_parts)

  x <- safe_kegg(entrez)
  if (nrow(x) > 0) x$major.band <- band
  kegg_all[[band]] <- x
}

go_results <- bind_rows(go_all)
kegg_results <- bind_rows(kegg_all)

# Reciprocal enrichment against the MSigDB C1 positional collection.
c1 <- tryCatch(
  msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "C1"
  ),
  error = function(e) {
    msigdbr::msigdbr(species = "Homo sapiens", category = "C1")
  }
)

gene_col <- if ("ncbi_gene" %in% names(c1)) "ncbi_gene" else "entrez_gene"
term2gene <- c1 |>
  transmute(gs_name = gs_name, entrez_gene = as.character(.data[[gene_col]])) |>
  filter(!is.na(entrez_gene), nzchar(entrez_gene)) |>
  distinct()

clean_band_name <- function(x) {
  x <- sub("^CHR", "", toupper(x))
  x <- sub("_.*$", "", x)
  tolower(x)
}

reciprocal_kegg <- list()
if (nrow(kegg_results) > 0) {
  for (i in seq_len(nrow(kegg_results))) {
    pathway_id <- as.character(kegg_results$ID[i])
    band <- as.character(kegg_results$major.band[i])

    linked <- tryCatch(
      KEGGREST::keggLink("hsa", pathway_id),
      error = function(e) character(0)
    )
    entrez <- sub("^hsa:", "", as.character(linked))
    if (length(entrez) == 0) next

    enriched <- as.data.frame(
      clusterProfiler::enricher(
        gene = entrez,
        TERM2GENE = term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
      )
    )
    if (nrow(enriched) == 0) next

    enriched$band_from_C1 <- clean_band_name(enriched$Description)
    hit <- enriched[
      enriched$p.adjust < 0.05 &
        enriched$band_from_C1 == tolower(band),
      ,
      drop = FALSE
    ]
    if (nrow(hit) > 0) {
      hit$KEGG.ID <- pathway_id
      reciprocal_kegg[[length(reciprocal_kegg) + 1]] <- hit
    }
  }
}

reciprocal_go <- list()
if (nrow(go_results) > 0) {
  for (i in seq_len(nrow(go_results))) {
    go_id <- as.character(go_results$ID[i])
    band <- as.character(go_results$major.band[i])

    genes <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go_id,
        columns = "ENTREZID",
        keytype = "GOALL"
      )
    )
    entrez <- unique(na.omit(as.character(genes$ENTREZID)))
    if (length(entrez) == 0) next

    enriched <- as.data.frame(
      clusterProfiler::enricher(
        gene = entrez,
        TERM2GENE = term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
      )
    )
    if (nrow(enriched) == 0) next

    enriched$band_from_C1 <- clean_band_name(enriched$Description)
    hit <- enriched[
      enriched$p.adjust < 0.05 &
        enriched$band_from_C1 == tolower(band),
      ,
      drop = FALSE
    ]
    if (nrow(hit) > 0) {
      hit$GO.ID <- go_id
      hit$GO.Type <- as.character(go_results$GO[i])
      reciprocal_go[[length(reciprocal_go) + 1]] <- hit
    }
  }
}

reciprocal_kegg <- bind_rows(reciprocal_kegg)
reciprocal_go <- bind_rows(reciprocal_go)

# Keep the same four-sheet layout expected by analysis.R and the final paper:
# Table S1 = GO, Table S2 = KEGG, Table S3 = reciprocal GO, Table S4 = reciprocal KEGG.
# The first four rows are intentionally blank because analysis.R reads headers from row 5.
outfile <- file.path(ENRICHMENT_DIR, "Supplementary_Tables_current.xlsx")
out_wb <- openxlsx::createWorkbook()
for (sheet in c("Table S1", "Table S2", "Table S3", "Table S4")) {
  openxlsx::addWorksheet(out_wb, sheet)
}
openxlsx::writeData(out_wb, "Table S1", go_results, startRow = 5, rowNames = FALSE)
openxlsx::writeData(out_wb, "Table S2", kegg_results, startRow = 5, rowNames = FALSE)
openxlsx::writeData(out_wb, "Table S3", reciprocal_go, startRow = 5, rowNames = FALSE)
openxlsx::writeData(out_wb, "Table S4", reciprocal_kegg, startRow = 5, rowNames = FALSE)
openxlsx::saveWorkbook(out_wb, outfile, overwrite = TRUE)

capture.output(
  sessionInfo(),
  file = file.path(ENRICHMENT_DIR, "sessionInfo.txt")
)

message("Refresh complete.")
message("Band gene lists: ", BAND_GENE_DIR)
message("GO/KEGG workbook: ", outfile)
message("Run source(\"analysis.R\") again to use these refreshed files.")
