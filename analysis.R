# ==============================================================================
# CRC GWAS chromosomal-band analysis
# ==============================================================================
# Reproducible analysis accompanying:
# "The possibility of prognostic and functional values of the 8q24 and 20q13
#  chromosomal bands in colorectal cancer"
#
# The default run uses the archived files retained with the study. Optional
# scripts under scripts/ can refresh database-dependent steps; when refreshed
# files are present under data/refresh/, this script uses them automatically.
#
# Run from R or RStudio with:
#   source("analysis.R")
# ==============================================================================

# ---- Package setup ------------------------------------------------------------

cran_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "patchwork",
  "scales"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("chromPlot", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  message("Installing Bioconductor package: chromPlot")
  BiocManager::install("chromPlot", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(scales)
})


# ---- Project paths ------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_file <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) {
      value <- x$ofile
      if (is.null(value)) NA_character_ else as.character(value)[1]
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]

  if (length(frame_files) > 0) {
    return(dirname(normalizePath(tail(frame_files, 1), winslash = "/", mustWork = TRUE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_file <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(active_file)) {
      return(dirname(normalizePath(active_file, winslash = "/", mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

find_project_dir <- function() {
  script_dir <- get_script_dir()
  candidates <- unique(c(
    script_dir,
    getwd(),
    dirname(script_dir),
    dirname(getwd())
  ))

  marker <- file.path(
    "data", "raw", "colorectal_cancer_GWAS_MONDO_0005575.xlsx"
  )

  hits <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hits) == 0) {
    stop(
      "Could not find the project data folder.\n",
      "Keep analysis.R in the project root, next to the data/ folder.",
      call. = FALSE
    )
  }

  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

PROJECT_DIR <- find_project_dir()
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
ARCHIVED_DIR <- file.path(PROJECT_DIR, "data", "archived")
REFRESH_DIR <- file.path(PROJECT_DIR, "data", "refresh")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
SUPP_FIGURES_DIR <- file.path(FIGURES_DIR, "supplementary")

invisible(lapply(
  c(RESULTS_DIR, FIGURES_DIR, TABLES_DIR, SUPP_FIGURES_DIR),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

prefer_refresh <- function(refreshed, archived) {
  if (file.exists(refreshed)) refreshed else archived
}

refresh_mutation_dir <- file.path(REFRESH_DIR, "mutation")
refresh_mutation_files <- c(
  file.path(refresh_mutation_dir, "Table_4_mutation_enrichment.xlsx"),
  file.path(refresh_mutation_dir, "Table_5_mutation_ttest.xlsx")
)
refresh_mutation_ready <- all(file.exists(refresh_mutation_files))

PUBLISHED_SUPPLEMENT <- file.path(
  ARCHIVED_DIR, "enrichment", "Supplementary_Tables_published.xlsx"
)

files <- list(
  gwas = file.path(RAW_DIR, "colorectal_cancer_GWAS_MONDO_0005575.xlsx"),
  cytoband = file.path(RAW_DIR, "cytoBand_hg38.txt"),
  deg_up = file.path(RAW_DIR, "COAD_UALCAN_top250_upregulated.txt"),
  deg_down = file.path(RAW_DIR, "COAD_UALCAN_top250_downregulated.txt"),
  supplement = prefer_refresh(
    file.path(REFRESH_DIR, "enrichment", "Supplementary_Tables_current.xlsx"),
    PUBLISHED_SUPPLEMENT
  ),
  table1 = file.path(
    ARCHIVED_DIR, "published_tables", "Table_1_band_statistics.xlsx"
  ),
  table2 = prefer_refresh(
    file.path(REFRESH_DIR, "csurvival", "Table_2_survival.xlsx"),
    PUBLISHED_SUPPLEMENT
  ),
  table3 = file.path(
    ARCHIVED_DIR, "published_tables", "Table_3_DEG_enrichment.xlsx"
  ),
  table4 = if (refresh_mutation_ready) {
    refresh_mutation_files[1]
  } else {
    file.path(ARCHIVED_DIR, "published_tables", "Table_4_mutation_enrichment.xlsx")
  },
  table5 = if (refresh_mutation_ready) {
    refresh_mutation_files[2]
  } else {
    file.path(ARCHIVED_DIR, "published_tables", "Table_5_mutation_ttest.xlsx")
  }
)

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop(
    "Required project file(s) are missing:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    call. = FALSE
  )
}


# ---- Small helpers ------------------------------------------------------------

save_ggplot <- function(plot, name, width, height, directory = FIGURES_DIR, dpi = 600) {
  ggsave(
    file.path(directory, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    file.path(directory, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

save_base_plot <- function(draw_fun, name, width, height, directory = FIGURES_DIR, dpi = 600) {
  pdf(
    file.path(directory, paste0(name, ".pdf")),
    width = width,
    height = height
  )
  draw_fun()
  dev.off()

  png(
    file.path(directory, paste0(name, ".png")),
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  draw_fun()
  dev.off()
}

write_tsv <- function(x, filename) {
  write.table(
    x,
    file.path(TABLES_DIR, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

clean_names <- function(x) {
  names(x) <- make.names(names(x), unique = TRUE)
  x
}

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE)
  genes <- trimws(genes)
  unique(genes[nzchar(genes)])
}


# ---- 1. Read the archived input data -----------------------------------------

gwas_all <- openxlsx::read.xlsx(files$gwas) |>
  clean_names()

required_gwas_columns <- c(
  "REGION", "CHR_ID", "CHR_POS", "SNPS", "Gene.Type", "Gene.Name",
  "Maped.downe.gene", "Distance", "Type", "P.VALUE"
)

missing_columns <- setdiff(required_gwas_columns, names(gwas_all))
if (length(missing_columns) > 0) {
  stop(
    "The GWAS file is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# The manuscript started from 1,796 associations and retained one row per SNP.
gwas <- gwas_all[!duplicated(gwas_all$SNPS), , drop = FALSE]
gwas$major.band <- sub("\\..*$", "", gwas$REGION)
gwas$minor.band <- gwas$REGION

if (nrow(gwas_all) != 1796) {
  warning("Expected 1,796 archived GWAS rows; found ", nrow(gwas_all), ".")
}
if (nrow(gwas) != 1346) {
  warning("Expected 1,346 unique SNPs; found ", nrow(gwas), ".")
}

cytoband <- read.delim(
  files$cytoband,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cytoband_columns <- c("Chr", "Start", "End", "Band")
if (!all(required_cytoband_columns %in% names(cytoband))) {
  stop("cytoBand_hg38.txt does not have the expected columns.", call. = FALSE)
}

cytoband <- cytoband[cytoband$Band != "", , drop = FALSE]


# ---- 2. Chromosomal-band distribution test -----------------------------------

# Fixed denominator used in the original working script and the published
# chromosomal-band test. It is kept as a publication-reproduction constant
# rather than being recomputed from the cytoband file.
PUBLISHED_GENOME_LENGTH_BP <- 3.2e9

major_counts <- as.data.frame(table(gwas$major.band), stringsAsFactors = FALSE)
colnames(major_counts) <- c("major.band", "ni")
major_counts$ni <- as.numeric(major_counts$ni)

minor_counts <- as.data.frame(table(gwas$minor.band), stringsAsFactors = FALSE)
colnames(minor_counts) <- c("minor.band", "ni")
minor_counts$ni <- as.numeric(minor_counts$ni)

cytoband$major.band <- paste0(
  sub("^chr", "", cytoband$Chr),
  sub("\\..*$", "", cytoband$Band)
)
cytoband$minor.band <- paste0(sub("^chr", "", cytoband$Chr), cytoband$Band)

major_lengths <- cytoband |>
  group_by(major.band) |>
  summarise(
    Start = min(Start),
    End = max(End),
    .groups = "drop"
  ) |>
  mutate(Length = End - Start)

minor_lengths <- cytoband |>
  transmute(
    minor.band = minor.band,
    Start = Start,
    End = End,
    Length = End - Start
  )

calculate_band_statistics <- function(counts, lengths, band_column) {
  result <- merge(counts, lengths, by = band_column)
  result$n <- nrow(gwas)
  result$theta <- result$Length / PUBLISHED_GENOME_LENGTH_BP
  result$df1 <- 2 * (result$n - result$ni + 1)
  result$df2 <- 2 * result$ni
  result$F.test <- ((1 - result$theta) / result$theta) *
    (result$ni / (result$n - result$ni + 1))
  result$pvalue <- pf(result$F.test, result$df1, result$df2, lower.tail = FALSE)

  # This is the same multiplicative correction used in the published script.
  result$adj.pvalue <- result$pvalue * nrow(result)
  result
}

major_stats <- calculate_band_statistics(
  major_counts,
  major_lengths,
  "major.band"
)
minor_stats <- calculate_band_statistics(
  minor_counts,
  minor_lengths,
  "minor.band"
)

significant_major <- major_stats |>
  filter(adj.pvalue < 0.05) |>
  arrange(major.band)

significant_minor <- minor_stats |>
  filter(adj.pvalue < 0.05) |>
  arrange(minor.band)

expected_major_bands <- c(
  "10p14", "10q25", "11q12", "12p13", "15q13", "18q21", "19q13",
  "1q41", "20p12", "20q13", "6p21", "8q24", "9q34"
)

if (!setequal(significant_major$major.band, expected_major_bands)) {
  warning(
    "The significant major-band list differs from the published 13-band result."
  )
}

openxlsx::write.xlsx(
  major_stats,
  file.path(TABLES_DIR, "major_band_statistics.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  minor_stats,
  file.path(TABLES_DIR, "minor_band_statistics.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  significant_major,
  file.path(TABLES_DIR, "Table_1_significant_major_bands.xlsx"),
  rowNames = FALSE
)

write_tsv(significant_major, "Table_1_significant_major_bands.tsv")


# ---- 3. Figure 1A: Manhattan plot of the 1,796 archived associations ---------

manhattan_data <- gwas_all |>
  transmute(
    SNP = SNPS,
    CHR = as.character(CHR_ID),
    BP = as.character(CHR_POS),
    P = suppressWarnings(as.numeric(P.VALUE))
  )

# These two coordinates were corrected manually in the original analysis.
manhattan_data$BP[manhattan_data$SNP == "rs71167281"] <- "222034660"
manhattan_data$BP[manhattan_data$SNP == "rs67052019"] <- "109822839"
manhattan_data$BP <- suppressWarnings(as.numeric(manhattan_data$BP))
manhattan_data$logP <- -log10(manhattan_data$P)

chromosome_lengths <- cytoband |>
  group_by(Chr) |>
  summarise(chr_length = max(End), .groups = "drop") |>
  mutate(CHR = sub("^chr", "", Chr)) |>
  filter(CHR %in% c(as.character(1:22), "X")) |>
  mutate(order = match(CHR, c(as.character(1:22), "X"))) |>
  arrange(order) |>
  mutate(
    cumulative_start = c(0, head(cumsum(as.numeric(chr_length)), -1)),
    chr_center = cumulative_start + chr_length / 2
  )

manhattan_data <- manhattan_data |>
  left_join(
    chromosome_lengths |>
      dplyr::select(CHR, cumulative_start),
    by = "CHR"
  ) |>
  filter(!is.na(BP), !is.na(P), P > 0, !is.na(cumulative_start)) |>
  mutate(
    CHR_POS = BP + cumulative_start,
    CHR = factor(CHR, levels = c(as.character(1:22), "X"))
  )

manhattan_colors <- rep(c("blue", "red"), length.out = 23)
names(manhattan_colors) <- c(as.character(1:22), "X")

p_manhattan <- ggplot(
  manhattan_data,
  aes(x = CHR_POS, y = logP, color = CHR)
) +
  geom_point(alpha = 0.65, size = 0.8) +
  scale_color_manual(values = manhattan_colors, guide = "none") +
  scale_x_continuous(
    breaks = chromosome_lengths$chr_center,
    labels = chromosome_lengths$CHR,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Chromosomes",
    y = expression(-log[10](italic(P) * "-value"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  )

save_ggplot(p_manhattan, "Figure_1A_Manhattan", width = 14, height = 6)


# ---- 4. Figure 1B and Figure 3: chromosomal views -----------------------------

utils::data("hg_cytoBandIdeo", package = "chromPlot", envir = environment())
utils::data("hg_gap", package = "chromPlot", envir = environment())

all_band_plot_data <- major_stats |>
  dplyr::select(major.band, ni, Start, End) |>
  mutate(Chrom = paste0("chr", sub("[a-z].*$", "", major.band)))

all_band_annotations <- all_band_plot_data[
  rep(seq_len(nrow(all_band_plot_data)), all_band_plot_data$ni),
  c("major.band", "Chrom", "Start", "End")
]
rownames(all_band_annotations) <- NULL
colnames(all_band_annotations) <- c("Name", "Chrom", "Start", "End")
all_band_annotations <- all_band_annotations[, c("Chrom", "Start", "End", "Name")]

save_base_plot(
  draw_fun = function() {
    chromPlot::chromPlot(
      gaps = hg_gap,
      bands = hg_cytoBandIdeo,
      annot1 = all_band_annotations,
      figCols = 8
    )
  },
  name = "Figure_1B_GWAS_band_histogram",
  width = 10,
  height = 15
)

significant_band_plot_data <- significant_major |>
  dplyr::select(major.band, Start, End) |>
  mutate(Chrom = sub("[a-z].*$", "", major.band))

colnames(significant_band_plot_data) <- c("ID", "Start", "End", "Chrom")
significant_band_stat <- significant_band_plot_data[, c("Chrom", "Start", "End", "ID")]
significant_band_stat$Value <- 1

save_base_plot(
  draw_fun = function() {
    chromPlot::chromPlot(
      gaps = hg_gap,
      bands = hg_cytoBandIdeo,
      stat = significant_band_stat,
      statCol = "Value",
      statName = "Value",
      noHist = TRUE,
      figCols = 8,
      cex = 0.7,
      statTyp = "n",
      chrSide = c(1, 1, 1, 1, 1, 1, -1, 1)
    )
  },
  name = "Figure_3_significant_chromosomal_bands",
  width = 10,
  height = 15
)


# ---- 5. Figure 2 and additional band-level variant distributions -------------

summarise_variant_classes <- function(df) {
  gene_type <- as.character(df$Gene.Type)
  consequence <- as.character(df$Type)

  has_protein_coding <- grepl("protein coding", gene_type, fixed = TRUE)

  # Some variants carry more than one consequence label. For Figure 2 the
  # coding subclasses are mutually exclusive, matching the published counts.
  other_coding <- has_protein_coding &
    grepl(
      "Synonymous Variant|Splice Acceptor Variant|Stop Gained",
      consequence
    )
  missense <- has_protein_coding &
    grepl("Missense Variant", consequence, fixed = TRUE)
  utr <- has_protein_coding &
    grepl("5 Prime UTR Variant|3 Prime UTR Variant", consequence)
  intronic <- has_protein_coding &
    grepl("Intron Variant", consequence, fixed = TRUE) &
    !(other_coding | missense | utr)

  functional_coding <- intronic | missense | utr | other_coding
  uncharacterized <- trimws(gene_type) == "0"
  annotated_noncoding <- !functional_coding & !uncharacterized

  data.frame(
    total = nrow(df),
    mapped_to_protein_coding = sum(has_protein_coding),
    not_mapped_to_protein_coding = sum(!has_protein_coding),
    coding = sum(functional_coding),
    noncoding = sum(!functional_coding),
    intronic = sum(intronic),
    missense = sum(missense),
    utr = sum(utr),
    other_coding = sum(other_coding),
    annotated_noncoding = sum(annotated_noncoding),
    uncharacterized_noncoding = sum(uncharacterized)
  )
}

variant_summary <- summarise_variant_classes(gwas)

# Two related definitions occur in the published material:
# 637 variants map to at least one protein-coding gene, while 632 have one of the
# coding consequence classes shown in Figure 2. The remaining five are annotated
# as downstream variants and are shown with the non-coding consequence group.
write_tsv(variant_summary, "variant_classification_summary.tsv")

variant_colors <- c(
  "Polymorphisms" = "#FFFFFF",
  "Coding" = "#5E0075",
  "Intronic" = "#8B00A3",
  "Missense" = "#B300E0",
  "Others (Synonymous, Splice Acceptor, and Stop Gained)" = "#D633FF",
  "UTRs (5 prime UTR and 3 prime UTR)" = "#F07CFF",
  "Non-Coding" = "#005F73",
  "Non-coding (annotated)" = "#008C99",
  "Non-coding (uncharacterized)" = "#00C4CC"
)

make_nested_variant_plot <- function(
    summary_row,
    title = NULL,
    show_labels = TRUE,
    legend_breaks = names(variant_colors)
) {
  total <- summary_row$total

  plot_data <- data.frame(
    Category = c(
      "Polymorphisms",
      "Coding",
      "Non-Coding",
      "Intronic",
      "Missense",
      "Others (Synonymous, Splice Acceptor, and Stop Gained)",
      "UTRs (5 prime UTR and 3 prime UTR)",
      "Non-coding (annotated)",
      "Non-coding (uncharacterized)"
    ),
    Count = c(
      total,
      summary_row$coding,
      summary_row$noncoding,
      summary_row$intronic,
      summary_row$missense,
      summary_row$other_coding,
      summary_row$utr,
      summary_row$annotated_noncoding,
      summary_row$uncharacterized_noncoding
    ),
    Ring = c(1, 2, 2, 3, 3, 3, 3, 3, 3),
    stringsAsFactors = FALSE
  )

  coding_children <- c(
    "Intronic",
    "Missense",
    "Others (Synonymous, Splice Acceptor, and Stop Gained)",
    "UTRs (5 prime UTR and 3 prime UTR)"
  )

  # Keep the same factor levels in every panel so patchwork can collect one
  # consistent legend even when a category is absent from a particular band.
  plot_data$Category <- factor(
    plot_data$Category,
    levels = names(variant_colors)
  )

  plot_data <- plot_data[plot_data$Count > 0, , drop = FALSE]
  plot_data$Percentage <- ifelse(
    plot_data$Ring == 1,
    100,
    ifelse(
      plot_data$Ring == 2,
      100 * plot_data$Count / total,
      ifelse(
        plot_data$Category %in% coding_children,
        100 * plot_data$Count / summary_row$coding,
        100 * plot_data$Count / summary_row$noncoding
      )
    )
  )
  plot_data$Label <- ifelse(
    show_labels & plot_data$Ring > 1 & plot_data$Percentage >= 2.5,
    paste0(round(plot_data$Percentage, 1), "%"),
    ""
  )

  p <- ggplot(plot_data, aes(x = Ring, y = Count, fill = Category)) +
    geom_col(width = 0.98, color = NA) +
    geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = if (show_labels) 4 else 3
    ) +
    coord_polar(theta = "y") +
    xlim(0.25, 3.55) +
    scale_fill_manual(
      values = variant_colors,
      breaks = legend_breaks,
      drop = FALSE
    ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  p
}

p_variant_distribution <- make_nested_variant_plot(variant_summary, show_labels = TRUE)
save_ggplot(
  p_variant_distribution,
  "Figure_2_variant_distribution",
  width = 10,
  height = 7
)

band_variant_plots <- list()
for (band in expected_major_bands) {
  band_summary <- summarise_variant_classes(
    gwas[gwas$major.band == band, , drop = FALSE]
  )

  # Legends are suppressed in every individual panel. A single, independent
  # legend is drawn below the complete multi-band figure instead.
  band_variant_plots[[band]] <- make_nested_variant_plot(
    band_summary,
    title = band,
    show_labels = FALSE
  ) +
    theme(legend.position = "none")
}

# Thirteen panels are arranged as 5 + 5 + 3. Two spacers center the last row.
supplementary_panels <- c(
  band_variant_plots[1:10],
  list(patchwork::plot_spacer()),
  band_variant_plots[11:13],
  list(patchwork::plot_spacer())
)

p_supplementary_panels <- patchwork::wrap_plots(
  supplementary_panels,
  ncol = 5
)

# One manual legend for the whole figure. Using a separate legend panel avoids
# repeated or clipped guides when some chromosomal bands do not contain every
# variant class. This is an additional repository figure, not the published Figure S1.
supplementary_legend <- data.frame(
  key = c(
    "Coding",
    "Non-Coding",
    "Intronic",
    "Missense",
    "UTRs (5 prime UTR and 3 prime UTR)",
    "Others (Synonymous, Splice Acceptor, and Stop Gained)",
    "Non-coding (annotated)",
    "Non-coding (uncharacterized)"
  ),
  label = c(
    "Coding",
    "Non-coding",
    "Intronic",
    "Missense",
    "UTR variants (5' and 3' UTR)",
    "Other coding variants\n(synonymous, splice acceptor, stop gained)",
    "Non-coding (annotated)",
    "Non-coding (uncharacterized)"
  ),
  x = rep(c(1.0, 5.2, 9.4, 13.6), 2),
  y = c(rep(2, 4), rep(1, 4)),
  stringsAsFactors = FALSE
)
supplementary_legend$x_text <- supplementary_legend$x + 0.30

p_supplementary_legend <- ggplot(
  supplementary_legend,
  aes(x = x, y = y)
) +
  geom_point(
    aes(fill = key),
    shape = 22,
    size = 5.2,
    stroke = 0
  ) +
  geom_text(
    aes(x = x_text, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = variant_colors,
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0.65, 17.7),
    ylim = c(0.45, 2.55),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 4, 2, 4)
  )

p_supplementary <- (
  p_supplementary_panels /
    p_supplementary_legend
) +
  patchwork::plot_layout(
    heights = c(10, 1.25)
  )

save_ggplot(
  p_supplementary,
  "Additional_Figure_band_variant_distribution",
  width = 14,
  height = 10.8,
  directory = SUPP_FIGURES_DIR,
  dpi = 400
)


# ---- 6. GO/KEGG enrichment results and Figures 4-5 ---------------------------

# The default source is the original Supplementary workbook published with the
# article. Its Tables S1-S4 use a presentation-oriented layout. To avoid relying
# on how a particular openxlsx version sanitizes those header labels, the
# published sheets are read by column position and then given explicit internal
# names. A refreshed workbook produced by scripts/01_refresh_GO_KEGG.R already
# uses the internal R column names and is read directly.

using_refreshed_enrichment <- grepl(
  "/data/refresh/",
  gsub("\\\\", "/", files$supplement)
)

ratio_count <- function(x) {
  suppressWarnings(as.numeric(sub("/.*$", "", as.character(x))))
}

read_published_sheet <- function(sheet_name, n_columns) {
  x <- openxlsx::read.xlsx(
    files$supplement,
    sheet = sheet_name,
    startRow = 4,
    colNames = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = TRUE
  )

  if (nrow(x) < 2 || ncol(x) < n_columns) {
    stop(
      "Could not read the expected layout from ", sheet_name,
      " in the published Supplementary workbook.",
      call. = FALSE
    )
  }

  # Row 4 contains the printed column labels. With colNames = FALSE it is read
  # as the first data row, so remove it and keep only the known table columns.
  x <- x[-1, seq_len(n_columns), drop = FALSE]
  rownames(x) <- NULL
  x
}

if (!using_refreshed_enrichment) {
  s1 <- read_published_sheet("Table S1", 8)
  go_results <- data.frame(
    ID = as.character(s1[[2]]),
    Description = as.character(s1[[4]]),
    GeneRatio = as.character(s1[[5]]),
    BgRatio = as.character(s1[[6]]),
    p.adjust = suppressWarnings(as.numeric(s1[[7]])),
    geneID = as.character(s1[[8]]),
    Count = ratio_count(s1[[5]]),
    GO = as.character(s1[[3]]),
    major.band = as.character(s1[[1]]),
    stringsAsFactors = FALSE
  )

  s2 <- read_published_sheet("Table S2", 9)
  kegg_results <- data.frame(
    category = as.character(s2[[3]]),
    subcategory = as.character(s2[[4]]),
    ID = as.character(s2[[2]]),
    Description = as.character(s2[[5]]),
    GeneRatio = as.character(s2[[6]]),
    BgRatio = as.character(s2[[7]]),
    p.adjust = suppressWarnings(as.numeric(s2[[8]])),
    geneID = as.character(s2[[9]]),
    Count = ratio_count(s2[[6]]),
    major.band = as.character(s2[[1]]),
    stringsAsFactors = FALSE
  )

  s3 <- read_published_sheet("Table S3", 7)
  mutual_go <- data.frame(
    GO.ID = as.character(s3[[1]]),
    GO.Type = as.character(s3[[2]]),
    Description = as.character(s3[[3]]),
    GeneRatio = as.character(s3[[4]]),
    BgRatio = as.character(s3[[5]]),
    p.adjust = suppressWarnings(as.numeric(s3[[6]])),
    geneID = as.character(s3[[7]]),
    Count = ratio_count(s3[[4]]),
    major.band = as.character(s3[[3]]),
    stringsAsFactors = FALSE
  )

  s4 <- read_published_sheet("Table S4", 6)
  mutual_kegg <- data.frame(
    KEGG.ID = as.character(s4[[1]]),
    Description = as.character(s4[[2]]),
    GeneRatio = as.character(s4[[3]]),
    BgRatio = as.character(s4[[4]]),
    p.adjust = suppressWarnings(as.numeric(s4[[5]])),
    geneID = as.character(s4[[6]]),
    Count = ratio_count(s4[[3]]),
    major.band = as.character(s4[[2]]),
    stringsAsFactors = FALSE
  )
} else {
  read_refreshed_sheet <- function(sheet_name) {
    openxlsx::read.xlsx(
      files$supplement,
      sheet = sheet_name,
      startRow = 5,
      check.names = FALSE
    )
  }

  go_results <- read_refreshed_sheet("Table S1")
  kegg_results <- read_refreshed_sheet("Table S2")
  mutual_go <- read_refreshed_sheet("Table S3")
  mutual_kegg <- read_refreshed_sheet("Table S4")

  if (!"Count" %in% names(go_results) && "GeneRatio" %in% names(go_results)) {
    go_results$Count <- ratio_count(go_results$GeneRatio)
  }
  if (!"Count" %in% names(kegg_results) && "GeneRatio" %in% names(kegg_results)) {
    kegg_results$Count <- ratio_count(kegg_results$GeneRatio)
  }
  if (!"Count" %in% names(mutual_go) && "GeneRatio" %in% names(mutual_go)) {
    mutual_go$Count <- ratio_count(mutual_go$GeneRatio)
  }
  if (!"Count" %in% names(mutual_kegg) && "GeneRatio" %in% names(mutual_kegg)) {
    mutual_kegg$Count <- ratio_count(mutual_kegg$GeneRatio)
  }

  required_go <- c("ID", "Description", "p.adjust", "GO", "major.band", "Count")
  required_kegg <- c("ID", "Description", "p.adjust", "major.band", "Count")
  required_mutual_go <- c("GO.ID", "GO.Type", "Description", "p.adjust", "Count")
  required_mutual_kegg <- c("KEGG.ID", "Description", "p.adjust", "Count")

  missing_refresh_columns <- c(
    setdiff(required_go, names(go_results)),
    setdiff(required_kegg, names(kegg_results)),
    setdiff(required_mutual_go, names(mutual_go)),
    setdiff(required_mutual_kegg, names(mutual_kegg))
  )
  if (length(missing_refresh_columns) > 0) {
    stop(
      "The refreshed enrichment workbook is missing expected columns: ",
      paste(unique(missing_refresh_columns), collapse = ", "),
      call. = FALSE
    )
  }
}

# The final publication contains 207 GO, 10 KEGG, 179 reciprocal GO and 7
# reciprocal KEGG rows. A current database refresh is allowed to differ.
if (!using_refreshed_enrichment) {
  expected_enrichment_rows <- c(
    GO = 207L,
    KEGG = 10L,
    reciprocal_GO = 179L,
    reciprocal_KEGG = 7L
  )
  observed_enrichment_rows <- c(
    GO = nrow(go_results),
    KEGG = nrow(kegg_results),
    reciprocal_GO = nrow(mutual_go),
    reciprocal_KEGG = nrow(mutual_kegg)
  )
  if (!all(observed_enrichment_rows == expected_enrichment_rows)) {
    warning(
      "The published enrichment workbook does not have the expected row counts: ",
      paste(names(observed_enrichment_rows), observed_enrichment_rows, sep = "=", collapse = ", ")
    )
  }
}

make_enrichment_plot <- function(df, title, y_title, n = 10) {
  top <- df |>
    filter(!is.na(p.adjust), !is.na(Description)) |>
    arrange(p.adjust) |>
    slice_head(n = n) |>
    mutate(Description = factor(Description, levels = rev(Description)))

  ggplot(
    top,
    aes(x = p.adjust, y = Description, size = Count, color = p.adjust)
  ) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = title,
      x = "Adjusted p-value",
      y = y_title,
      size = "Count",
      color = "p.adjust"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p_bp <- make_enrichment_plot(
  go_results |> filter(GO == "BP"),
  "Top 10 Biological Process (BP) Enrichment Dot Plot",
  "Biological Process"
)
p_mf <- make_enrichment_plot(
  go_results |> filter(GO == "MF"),
  "Top 10 Molecular Function (MF) Enrichment Dot Plot",
  "Molecular Function"
)
p_cc <- make_enrichment_plot(
  go_results |> filter(GO == "CC"),
  "Top 10 Cellular Component (CC) Enrichment Dot Plot",
  "Cellular Component"
)

p_go <- p_bp / p_mf / p_cc
save_ggplot(p_go, "Figure_4_GO_enrichment", width = 11, height = 16)

p_kegg <- make_enrichment_plot(
  kegg_results,
  "Top 10 KEGG Pathway Enrichment Dot Plot",
  "Pathway"
)
save_ggplot(p_kegg, "Figure_5_KEGG_enrichment", width = 10, height = 7)

openxlsx::write.xlsx(
  go_results,
  file.path(TABLES_DIR, "Table_S1_GO_enrichment.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  kegg_results,
  file.path(TABLES_DIR, "Table_S2_KEGG_enrichment.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  mutual_go,
  file.path(TABLES_DIR, "Table_S3_mutual_GO_band_enrichment.xlsx"),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  mutual_kegg,
  file.path(TABLES_DIR, "Table_S4_mutual_KEGG_band_enrichment.xlsx"),
  rowNames = FALSE
)


# ---- 7. DEG enrichment in significant chromosomal bands ----------------------

deg_up <- read.delim(
  files$deg_up,
  stringsAsFactors = FALSE,
  check.names = TRUE
)
deg_down <- read.delim(
  files$deg_down,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

if (!"Gene" %in% names(deg_up) || !"Gene" %in% names(deg_down)) {
  stop("The UALCAN DEG files do not contain a Gene column.", call. = FALSE)
}

coad_deg <- rbind(deg_up, deg_down)
coad_genes <- as.character(coad_deg$Gene)

refresh_band_gene_dir <- file.path(REFRESH_DIR, "band_genes")
refresh_band_gene_files <- file.path(
  refresh_band_gene_dir,
  paste0(expected_major_bands, "_genes.txt")
)

band_gene_dir <- if (all(file.exists(refresh_band_gene_files))) {
  refresh_band_gene_dir
} else {
  file.path(ARCHIVED_DIR, "band_genes")
}

band_genes <- setNames(
  lapply(
    expected_major_bands,
    function(band) {
      path <- file.path(band_gene_dir, paste0(band, "_genes.txt"))
      if (!file.exists(path)) {
        stop("Missing band gene list: ", path, call. = FALSE)
      }
      read_gene_list(path)
    }
  ),
  expected_major_bands
)

deg_enrichment <- do.call(
  rbind,
  lapply(expected_major_bands, function(band) {
    genes <- band_genes[[band]]
    overlap <- length(intersect(genes, coad_genes))

    p_value <- phyper(
      q = overlap - 1,
      m = length(genes),
      n = 25000 - length(genes),
      k = length(coad_genes),
      lower.tail = FALSE
    )

    data.frame(
      band = band,
      bands_gene = length(genes),
      DEGs = overlap,
      hyper_p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
)

# The Methods and the working R script specify BH correction. The archived
# published Table 3, however, contains the simple p * 13 correction. We keep
# both values so the difference is visible rather than silently replacing one.
deg_enrichment$adj_p_value_BH <- p.adjust(
  deg_enrichment$hyper_p_value,
  method = "BH"
)
deg_enrichment$adj_p_value_published <- pmin(
  deg_enrichment$hyper_p_value * length(expected_major_bands),
  1
)

deg_enrichment <- deg_enrichment |>
  arrange(match(band, expected_major_bands))

openxlsx::write.xlsx(
  deg_enrichment,
  file.path(TABLES_DIR, "Table_3_DEG_enrichment_recalculated.xlsx"),
  rowNames = FALSE
)
write_tsv(deg_enrichment, "Table_3_DEG_enrichment_recalculated.tsv")

published_deg_table <- openxlsx::read.xlsx(files$table3)
openxlsx::write.xlsx(
  published_deg_table,
  file.path(TABLES_DIR, "Table_3_DEG_enrichment_published.xlsx"),
  rowNames = FALSE
)

published_deg_sig <- c("20q13", "8q24")
bh_sig <- deg_enrichment$band[deg_enrichment$adj_p_value_BH < 0.05]
published_adjustment_sig <- deg_enrichment$band[
  deg_enrichment$adj_p_value_published < 0.05
]

if (!setequal(bh_sig, published_deg_sig) ||
    !setequal(published_adjustment_sig, published_deg_sig)) {
  warning("The DEG-enriched band list differs from the published result.")
}


# ---- 8. Survival and mutation results -----------------------------------------

# Survival was performed with the cSurvival web tool. Mutation analysis depended
# on TCGA/GDC plus Ensembl transcript lengths. Published/archived results are
# used unless refreshed files are present under data/refresh/.

survival_is_refreshed <- grepl(
  "/data/refresh/",
  gsub("\\\\", "/", files$table2)
)

# The publication result is read directly from Table S5 of the original
# Supplementary workbook. A refreshed cSurvival file is read from its first
# worksheet.
survival_table <- if (survival_is_refreshed) {
  openxlsx::read.xlsx(files$table2, check.names = FALSE)
} else {
  s5 <- openxlsx::read.xlsx(
    files$table2,
    sheet = "Table S5",
    startRow = 4,
    colNames = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = TRUE
  )
  if (nrow(s5) < 2 || ncol(s5) < 9) {
    stop("Could not read Table S5 from the published Supplementary workbook.", call. = FALSE)
  }
  s5 <- s5[-1, seq_len(9), drop = FALSE]
  names(s5) <- c(
    "Category",
    "Term ID",
    "Description",
    "Chromosomal Band",
    "Gene Ratio (Term Genes in Band / Total Term Genes)",
    "Adjusted p-value (Term Genes in Band)",
    "Gene Ratio (Band Genes in Term / Total Band Genes)",
    "Adjusted p-value (Band Genes in Term)",
    "Adjusted Survival p-value"
  )
  rownames(s5) <- NULL
  s5
}

mutation_hyper_table <- openxlsx::read.xlsx(files$table4)
mutation_ttest_table <- openxlsx::read.xlsx(files$table5)

survival_suffix <- if (survival_is_refreshed) {
  "refreshed"
} else {
  "published"
}
mutation_suffix <- if (grepl("/data/refresh/", gsub("\\\\", "/", files$table4))) {
  "refreshed"
} else {
  "archived"
}

openxlsx::write.xlsx(
  survival_table,
  file.path(TABLES_DIR, paste0("Table_2_survival_", survival_suffix, ".xlsx")),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  mutation_hyper_table,
  file.path(TABLES_DIR, paste0("Table_4_mutation_enrichment_", mutation_suffix, ".xlsx")),
  rowNames = FALSE
)
openxlsx::write.xlsx(
  mutation_ttest_table,
  file.path(TABLES_DIR, paste0("Table_5_mutation_ttest_", mutation_suffix, ".xlsx")),
  rowNames = FALSE
)


# ---- 9. Analysis summary and reproducibility ---------------------------------

# Use the BH-adjusted result for the recalculated summary. The published
# p-value x 13 correction identifies the same two bands (20q13 and 8q24).
current_deg_sig <- bh_sig

summary_table <- data.frame(
  metric = c(
    "Archived GWAS associations",
    "Unique SNPs after duplicate removal",
    "Variants mapped to at least one protein-coding gene",
    "Variants not mapped to a protein-coding gene",
    "Coding-consequence variants shown in Figure 2",
    "Non-coding-consequence variants shown in Figure 2",
    "Significant major chromosomal bands",
    "Significant DEG-enriched bands"
  ),
  value = c(
    nrow(gwas_all),
    nrow(gwas),
    variant_summary$mapped_to_protein_coding,
    variant_summary$not_mapped_to_protein_coding,
    variant_summary$coding,
    variant_summary$noncoding,
    nrow(significant_major),
    length(current_deg_sig)
  )
)

write_tsv(summary_table, "analysis_summary.tsv")

data_sources_used <- data.frame(
  item = c(
    "GO/KEGG enrichment workbook",
    "Band gene lists",
    "cSurvival table",
    "Mutation enrichment table",
    "Mutation t-test table"
  ),
  source = c(
    files$supplement,
    band_gene_dir,
    files$table2,
    files$table4,
    files$table5
  ),
  stringsAsFactors = FALSE
)
write_tsv(data_sources_used, "data_sources_used.tsv")

capture.output(
  sessionInfo(),
  file = file.path(RESULTS_DIR, "sessionInfo.txt")
)

message("\n============================================================")
message("Analysis complete")
message("Project: ", PROJECT_DIR)
message("Figures: ", FIGURES_DIR)
message("Tables:  ", TABLES_DIR)
message("============================================================\n")
