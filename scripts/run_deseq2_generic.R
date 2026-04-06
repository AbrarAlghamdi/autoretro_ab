args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 7) {

  stop("Usage: Rscript run_deseq2_generic.R <counts.tsv> <samples.tsv> <comparison.tsv> <output.csv> <volcano.png> <alpha> <log2fc_threshold>")

}


counts_file <- args[1]

samples_file <- args[2]

comparison_file <- args[3]

output_csv <- args[4]

volcano_png <- args[5]

alpha <- as.numeric(args[6])

log2fc_threshold <- as.numeric(args[7])


suppressPackageStartupMessages({

  library(DESeq2)

  library(ggplot2)

})


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


# keep only samples present in counts

sample_ids <- intersect(samples$sample_id, colnames(counts))

samples <- samples[samples$sample_id %in% sample_ids, , drop = FALSE]

samples <- samples[match(sample_ids, samples$sample_id), , drop = FALSE]


# subset counts

count_mat <- counts[, c("locus_id", sample_ids), drop = FALSE]

rownames(count_mat) <- count_mat$locus_id

count_mat$locus_id <- NULL


count_mat <- as.matrix(count_mat)

mode(count_mat) <- "numeric"

count_mat <- round(count_mat)


coldata <- data.frame(

  row.names = samples$sample_id,

  condition = factor(samples$condition)

)


# set reference level to group2 so log2FC = group1 vs group2

coldata$condition <- relevel(coldata$condition, ref = group2)


dds <- DESeqDataSetFromMatrix(

  countData = count_mat,

  colData = coldata,

  design = ~ condition

)


dds <- DESeq(dds)


res <- results(dds, contrast = c("condition", group1, group2))

res_df <- as.data.frame(res)

res_df$locus_id <- rownames(res_df)


res_df$significant <- ifelse(

  !is.na(res_df$padj) &

    res_df$padj < alpha &

    abs(res_df$log2FoldChange) >= log2fc_threshold,

  "yes",

  "no"

)


# reorder columns

res_df <- res_df[, c(

  "baseMean", "log2FoldChange", "lfcSE", "stat",

  "pvalue", "padj", "locus_id", "significant"

)]


write.csv(res_df, output_csv, row.names = FALSE)


# Volcano plot

plot_df <- res_df

plot_df$neglog10padj <- -log10(plot_df$padj)

plot_df$neglog10padj[is.infinite(plot_df$neglog10padj)] <- NA


p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10padj)) +

  geom_point(aes(color = significant), alpha = 0.7) +

  theme_minimal() +

  xlab("log2 Fold Change") +

  ylab("-log10 adjusted p-value") +

  ggtitle(paste0("Volcano plot: ", group1, " vs ", group2)) +

  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed") +

  geom_hline(yintercept = -log10(alpha), linetype = "dashed")


ggsave(volcano_png, plot = p, width = 7, height = 5, dpi = 150)


cat("Saved DESeq2 results:", output_csv, "\n")

cat("Saved volcano plot:", volcano_png, "\n")


