args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 5) {

  stop("Usage: Rscript family_enrichment.R <de_results.csv> <metadata.tsv> <output.csv> <alpha> <log2fc_threshold>")

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

if (!all(c("locus_id", "family") %in% colnames(meta))) {

  stop("Metadata must contain locus_id, family")

}


tested <- unique(de$locus_id)

de_loci <- unique(de$locus_id[

  !is.na(de$padj) & de$padj < alpha & abs(de$log2FoldChange) >= log2fc_threshold

])


meta2 <- meta[meta$locus_id %in% tested, , drop = FALSE]

meta2 <- unique(meta2[, c("locus_id", "family")])


families <- sort(unique(meta2$family))


results <- lapply(families, function(fam) {

  fam_tested_loci <- unique(meta2$locus_id[meta2$family == fam])


  a <- sum(de_loci %in% fam_tested_loci)                  # DE in family

  b <- sum(!(de_loci %in% fam_tested_loci))              # DE not in family

  c <- sum((tested %in% fam_tested_loci) & !(tested %in% de_loci))  # tested non-DE in family

  d <- sum(!(tested %in% fam_tested_loci) & !(tested %in% de_loci)) # tested non-DE not in family


  mat <- matrix(c(a, b, c, d), nrow = 2)

  ft <- fisher.test(mat)


  data.frame(

    family = fam,

    de_count = a,

    tested_count = length(fam_tested_loci),

    odds_ratio = unname(ft$estimate),

    pvalue = ft$p.value

  )

})


res <- do.call(rbind, results)

res$padj <- p.adjust(res$pvalue, method = "BH")

res <- res[order(res$padj, res$pvalue), , drop = FALSE]


write.csv(res, out_file, row.names = FALSE)

cat("Saved:", out_file, "\n")

