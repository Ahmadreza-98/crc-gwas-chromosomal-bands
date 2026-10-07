# Repeat the TCGA-COAD / TCGA-READ mutation analysis with the current GDC and
# Ensembl data. This is optional and can take a while because the mutation files
# are downloaded again.
#
# If refreshed Ensembl band gene lists exist, they are used. Otherwise the gene
# lists archived with the paper are used. The main analysis will automatically
# pick up the refreshed mutation tables after this script finishes.

cran_packages <- c("dplyr", "openxlsx")
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c("TCGAbiolinks", "maftools", "biomaRt")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages(library(dplyr))

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
ARCHIVED_GENE_DIR <- file.path(PROJECT_DIR, "data", "archived", "band_genes")
REFRESH_GENE_DIR <- file.path(PROJECT_DIR, "data", "refresh", "band_genes")
OUT_DIR <- file.path(PROJECT_DIR, "data", "refresh", "mutation")
GDC_DIR <- file.path(OUT_DIR, "gdc_download")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GDC_DIR, recursive = TRUE, showWarnings = FALSE)

bands <- c(
  "10p14", "10q25", "11q12", "12p13", "15q13", "18q21", "19q13",
  "1q41", "20p12", "20q13", "6p21", "8q24", "9q34"
)

refreshed_gene_files <- file.path(REFRESH_GENE_DIR, paste0(bands, "_genes.txt"))
GENE_DIR <- if (all(file.exists(refreshed_gene_files))) {
  REFRESH_GENE_DIR
} else {
  ARCHIVED_GENE_DIR
}

read_band_genes <- function(band) {
  path <- file.path(GENE_DIR, paste0(band, "_genes.txt"))
  x <- trimws(readLines(path, warn = FALSE))
  unique(x[nzchar(x)])
}

query_mutation <- function(project) {
  TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
}

message("Querying GDC...")
q_coad <- query_mutation("TCGA-COAD")
q_read <- query_mutation("TCGA-READ")

TCGAbiolinks::GDCdownload(q_coad, method = "api", directory = GDC_DIR)
TCGAbiolinks::GDCdownload(q_read, method = "api", directory = GDC_DIR)

coad <- TCGAbiolinks::GDCprepare(
  q_coad,
  directory = GDC_DIR,
  summarizedExperiment = FALSE
)
read <- TCGAbiolinks::GDCprepare(
  q_read,
  directory = GDC_DIR,
  summarizedExperiment = FALSE
)

mutation_data <- bind_rows(coad, read)
if (!"Hugo_Symbol" %in% names(mutation_data)) {
  stop("Prepared GDC mutation data do not contain Hugo_Symbol.", call. = FALSE)
}

maf <- maftools::read.maf(mutation_data)
gene_summary <- maftools::getGeneSummary(maf)

message("Retrieving current Ensembl transcript CDS lengths...")
mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

gene_map <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id"),
  filters = "hgnc_symbol",
  values = unique(mutation_data$Hugo_Symbol),
  mart = mart
)

cds_lengths <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters = "ensembl_transcript_id",
  values = unique(gene_map$ensembl_transcript_id),
  mart = mart
)

longest_transcripts <- merge(
  gene_map,
  cds_lengths,
  by = "ensembl_transcript_id"
) |>
  filter(!is.na(cds_length), nzchar(hgnc_symbol)) |>
  group_by(hgnc_symbol) |>
  slice_max(order_by = cds_length, n = 1) |>
  ungroup() |>
  transmute(Hugo_Symbol = hgnc_symbol, cds_length = cds_length)

mutation_summary <- merge(
  gene_summary,
  longest_transcripts,
  by = "Hugo_Symbol"
)

q1_length <- unname(quantile(mutation_summary$cds_length, probs = 0.25, na.rm = TRUE))
mutation_summary$total.adj <- (
  mutation_summary$total / (mutation_summary$cds_length + q1_length)
) * log2(mutation_summary$cds_length)
mutation_summary <- mutation_summary[!duplicated(mutation_summary), ]

# The original script used quantile(total.adj)[[3]]. With R's default quantile
# probabilities, that is the 50th percentile (the median).
high_cutoff <- unname(quantile(mutation_summary$total.adj, na.rm = TRUE)[[3]])
high_mutation <- mutation_summary |>
  filter(total.adj >= high_cutoff)

hyper_results <- list()
ttest_results <- list()

for (band in bands) {
  genes <- read_band_genes(band)
  in_universe <- intersect(genes, mutation_summary$Hugo_Symbol)
  high_overlap <- intersect(genes, high_mutation$Hugo_Symbol)

  hyper_p <- phyper(
    q = length(high_overlap) - 1,
    m = length(in_universe),
    n = nrow(mutation_summary) - length(in_universe),
    k = nrow(high_mutation),
    lower.tail = FALSE
  )

  hyper_results[[band]] <- data.frame(
    band = band,
    bands_gene = length(genes),
    bands_gene_in_mutataion_data = length(in_universe),
    high_mutation = length(high_overlap),
    hyper_p_value = hyper_p,
    stringsAsFactors = FALSE
  )

  band_values <- mutation_summary$total.adj[
    mutation_summary$Hugo_Symbol %in% genes
  ]
  background_values <- mutation_summary$total.adj[
    !mutation_summary$Hugo_Symbol %in% genes
  ]

  band_values <- band_values[is.finite(band_values)]
  background_values <- background_values[is.finite(background_values)]

  if (length(band_values) >= 2 && length(background_values) >= 2) {
    variance_test <- var.test(band_values, background_values)
    equal_variance <- is.finite(variance_test$p.value) && variance_test$p.value > 0.05

    tt <- t.test(
      band_values,
      background_values,
      alternative = "greater",
      var.equal = equal_variance
    )

    ttest_results[[band]] <- data.frame(
      band = band,
      f = unname(variance_test$statistic),
      f_test_p_value = variance_test$p.value,
      t = unname(tt$statistic),
      df = unname(tt$parameter),
      band_mean = mean(band_values),
      background_mean = mean(background_values),
      P_value = tt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

hyper_results <- bind_rows(hyper_results)
hyper_results$adj_p_value <- p.adjust(hyper_results$hyper_p_value, method = "BH")

ttest_results <- bind_rows(ttest_results)
if (nrow(ttest_results) > 0) {
  ttest_results$adj_p_value <- p.adjust(ttest_results$P_value, method = "BH")
}

openxlsx::write.xlsx(
  mutation_summary,
  file.path(OUT_DIR, "TCGA_COAD_READ_adjusted_mutation_burden_current.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  hyper_results,
  file.path(OUT_DIR, "Table_4_mutation_enrichment.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  ttest_results,
  file.path(OUT_DIR, "Table_5_mutation_ttest.xlsx"),
  rowNames = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "sessionInfo.txt")
)

message("Mutation refresh complete: ", OUT_DIR)
message("Band gene source: ", GENE_DIR)
message("Run source(\"analysis.R\") again to use the refreshed mutation tables.")
