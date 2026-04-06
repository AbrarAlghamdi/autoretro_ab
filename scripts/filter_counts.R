args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 5) {

  stop("Usage: Rscript filter_counts.R <counts.tsv> <samples.tsv> <comparison.tsv> <output.tsv> <min_count>")

}


counts_file <- args[1]

samples_file <- args[2]

comparison_file <- args[3]

output_file <- args[4]

min_count <- as.numeric(args[5])


counts <- read.delim(counts_file, check.names = FALSE)

samples <- read.delim(samples_file, check.names = FALSE)

comparison <- read.delim(comparison_file, check.names = FALSE)


if (!"locus_id" %in% colnames(counts)) {

  stop("counts.tsv must contain column: locus_id")

}

if (!all(c("sample_id", "condition") %in% colnames(samples))) {

  stop("samples.tsv must contain columns: sample_id, condition")

}

if (!all(c("group1", "group2") %in% colnames(comparison))) {

  stop("comparison.tsv must contain columns: group1, group2")

}


group1 <- comparison$group1[1]

group2 <- comparison$group2[1]


group1_samples <- samples$sample_id[samples$condition == group1]

group2_samples <- samples$sample_id[samples$condition == group2]


missing_cols <- setdiff(c(group1_samples, group2_samples), colnames(counts))

if (length(missing_cols) > 0) {

  stop(paste("Missing sample columns in counts:", paste(missing_cols, collapse = ", ")))

}


keep1 <- apply(counts[, group1_samples, drop = FALSE], 1, function(x) all(x >= min_count))

keep2 <- apply(counts[, group2_samples, drop = FALSE], 1, function(x) all(x >= min_count))

keep <- keep1 & keep2


filtered <- counts[keep, , drop = FALSE]


write.table(filtered, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved filtered counts:", output_file, "\n")

cat("Input rows:", nrow(counts), "\n")

cat("Kept rows:", nrow(filtered), "\n")
