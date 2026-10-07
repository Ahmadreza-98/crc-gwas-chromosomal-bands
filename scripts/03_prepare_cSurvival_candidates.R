# Prepare the GO/KEGG candidate list for a new cSurvival run.
#
# cSurvival is a web service, so this script does not perform survival analysis
# itself. It collects the reciprocally enriched GO/KEGG terms from whichever
# enrichment workbook is currently active and writes a simple workbook that can
# be used as a checklist when repeating the cSurvival step.
#
# After a new cSurvival analysis, save the final survival table as:
#   data/refresh/csurvival/Table_2_survival.xlsx
# analysis.R will then use it instead of the archived survival table.

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}

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
ARCHIVED_ENRICHMENT <- file.path(
  PROJECT_DIR,
  "data", "archived", "enrichment", "Supplementary_Tables_published.xlsx"
)
REFRESH_ENRICHMENT <- file.path(
  PROJECT_DIR,
  "data", "refresh", "enrichment", "Supplementary_Tables_current.xlsx"
)
OUT_DIR <- file.path(PROJECT_DIR, "data", "refresh", "csurvival")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ENRICHMENT_FILE <- if (file.exists(REFRESH_ENRICHMENT)) {
  REFRESH_ENRICHMENT
} else {
  ARCHIVED_ENRICHMENT
}

# The archived file is the original published Supplementary workbook. To keep
# this script independent of any header sanitization performed by openxlsx, the
# published Tables S3 and S4 are read by column position. A refreshed workbook
# already uses the internal R column names written by script 01.
is_refreshed <- grepl(
  "/data/refresh/",
  normalizePath(ENRICHMENT_FILE, winslash = "/", mustWork = TRUE)
)

ratio_count <- function(x) {
  suppressWarnings(as.numeric(sub("/.*$", "", as.character(x))))
}

if (!is_refreshed) {
  read_published_sheet <- function(sheet_name, n_columns) {
    x <- openxlsx::read.xlsx(
      ENRICHMENT_FILE,
      sheet = sheet_name,
      startRow = 4,
      colNames = FALSE,
      skipEmptyRows = TRUE,
      skipEmptyCols = TRUE
    )
    if (nrow(x) < 2 || ncol(x) < n_columns) {
      stop("Could not read ", sheet_name, " from the published Supplementary workbook.", call. = FALSE)
    }
    x <- x[-1, seq_len(n_columns), drop = FALSE]
    rownames(x) <- NULL
    x
  }

  s3 <- read_published_sheet("Table S3", 7)
  reciprocal_go <- data.frame(
    GO.ID = as.character(s3[[1]]),
    GO.Type = as.character(s3[[2]]),
    Description = as.character(s3[[3]]),
    p.adjust = suppressWarnings(as.numeric(s3[[6]])),
    geneID = as.character(s3[[7]]),
    Count = ratio_count(s3[[4]]),
    stringsAsFactors = FALSE
  )

  s4 <- read_published_sheet("Table S4", 6)
  reciprocal_kegg <- data.frame(
    KEGG.ID = as.character(s4[[1]]),
    Description = as.character(s4[[2]]),
    p.adjust = suppressWarnings(as.numeric(s4[[5]])),
    geneID = as.character(s4[[6]]),
    Count = ratio_count(s4[[3]]),
    stringsAsFactors = FALSE
  )
} else {
  reciprocal_go <- openxlsx::read.xlsx(
    ENRICHMENT_FILE,
    sheet = "Table S3",
    startRow = 5,
    check.names = FALSE
  )
  reciprocal_kegg <- openxlsx::read.xlsx(
    ENRICHMENT_FILE,
    sheet = "Table S4",
    startRow = 5,
    check.names = FALSE
  )

  if (!"Count" %in% names(reciprocal_go) && "GeneRatio" %in% names(reciprocal_go)) {
    reciprocal_go$Count <- ratio_count(reciprocal_go$GeneRatio)
  }
  if (!"Count" %in% names(reciprocal_kegg) && "GeneRatio" %in% names(reciprocal_kegg)) {
    reciprocal_kegg$Count <- ratio_count(reciprocal_kegg$GeneRatio)
  }
}

kegg_candidates <- unique(reciprocal_kegg[, intersect(
  c("KEGG.ID", "Description", "p.adjust", "geneID", "Count"),
  names(reciprocal_kegg)
), drop = FALSE])

go_candidates <- unique(reciprocal_go[, intersect(
  c("GO.ID", "GO.Type", "Description", "p.adjust", "geneID", "Count"),
  names(reciprocal_go)
), drop = FALSE])

out_file <- file.path(OUT_DIR, "cSurvival_candidates.xlsx")
openxlsx::write.xlsx(
  list(
    GO = go_candidates,
    KEGG = kegg_candidates
  ),
  out_file,
  rowNames = FALSE
)

writeLines(
  c(
    "cSurvival is run outside R at the cSurvival web server.",
    "Use cSurvival_candidates.xlsx as a checklist of the current reciprocal GO/KEGG terms.",
    "After the web analysis, save the final survival table as Table_2_survival.xlsx in this folder.",
    "The main analysis will automatically prefer that refreshed table over the archived one."
  ),
  file.path(OUT_DIR, "README.txt")
)

message("cSurvival candidate workbook written to: ", out_file)
message("Enrichment source: ", ENRICHMENT_FILE)
