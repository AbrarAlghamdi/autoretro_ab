args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 5) {

  stop("Usage: Rscript chromosomal_distribution.R <de_results.csv> <metadata.tsv> <output.csv> <alpha> <log2fc_threshold>")

}


de_file <- args[1]

meta_file <- args[2]

out_file <- args[3]

alpha <- as.numeric(args[4])

log2fc_threshold <- as.numeric(args[5])


de <- read.csv(de_file, check.names = FALSE)

meta <- read.delim(meta_file, check.names = FALSE)


if (!all(c("locus_id", "padj", "log2FoldChange") %in% colnames(de))) {

  stop("DE results must contain locus_id, padj, log2FoldChange")

}

if (!all(c("locus_id", "chromosome") %in% colnames(meta))) {

  stop("Metadata must contain locus_id, chromosome")

}


sig <- de[!is.na(de$padj) & de$padj < alpha & abs(de$log2FoldChange) >= log2fc_threshold, , drop = FALSE]

if (nrow(sig) == 0) {

  out <- data.frame(

    chromosome = character(),

    de_total = integer(),

    de_up = integer(),

    de_down = integer()

  )

  write.csv(out, out_file, row.names = FALSE)

  cat("Saved empty result:", out_file, "\n")

  quit(save = "no")

}


meta2 <- unique(meta[, c("locus_id", "chromosome")])

sig2 <- merge(sig[, c("locus_id", "log2FoldChange")], meta2, by = "locus_id", all.x = TRUE)


chroms <- sort(unique(sig2$chromosome))


out <- lapply(chroms, function(chr) {

  sub <- sig2[sig2$chromosome == chr, , drop = FALSE]

  data.frame(

    chromosome = chr,

    de_total = nrow(sub),

    de_up = sum(sub$log2FoldChange > 0, na.rm = TRUE),

    de_down = sum(sub$log2FoldChange < 0, na.rm = TRUE)

  )

})


out <- do.call(rbind, out)

out <- out[order(out$chromosome), , drop = FALSE]


write.csv(out, out_file, row.names = FALSE)

cat("Saved:", out_file, "\n")
