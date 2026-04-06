(بدون موضوع)
Abrar alghamdi
​
Abrar alghamdi​

Yes — here are the full files.





scripts/merge_telescope_counts.py


#!/usr/bin/env python3


import os

import sys

import glob

import pandas as pd



def find_report_file(sample_dir):

    candidates = [

        os.path.join(sample_dir, "telescope-telescope_report.tsv"),

        os.path.join(sample_dir, "telescope_report.tsv"),

    ]

    for c in candidates:

        if os.path.exists(c):

            return c

    return None



def load_counts(report_path):

    df = pd.read_csv(report_path, sep="\t", comment="#")

    required = {"transcript", "final_count"}

    missing = required - set(df.columns)

    if missing:

        raise ValueError(f"{report_path} missing columns: {missing}")


    df = df[["transcript", "final_count"]].copy()

    df = df[df["transcript"] != "__no_feature"]

    df.rename(columns={"transcript": "locus_id"}, inplace=True)

    return df



def main(input_root, output_tsv):

    sample_dirs = sorted(

        d for d in glob.glob(os.path.join(input_root, "*")) if os.path.isdir(d)

    )


    if not sample_dirs:

        raise FileNotFoundError(f"No sample directories found in: {input_root}")


    merged = None

    sample_names = []


    for sample_dir in sample_dirs:

        sample_name = os.path.basename(sample_dir)

        report = find_report_file(sample_dir)

        if report is None:

            print(f"Skipping {sample_name}: no telescope report found", file=sys.stderr)

            continue


        counts = load_counts(report)

        counts.rename(columns={"final_count": sample_name}, inplace=True)


        sample_names.append(sample_name)

        if merged is None:

            merged = counts

        else:

            merged = merged.merge(counts, on="locus_id", how="outer")


    if merged is None:

        raise RuntimeError("No valid Telescope reports found")


    merged.fillna(0, inplace=True)


    for c in sample_names:

        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)


    merged.sort_values("locus_id", inplace=True)

    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

    merged.to_csv(output_tsv, sep="\t", index=False)


    print(f"Saved merged counts to: {output_tsv}")

    print(f"Rows: {merged.shape[0]}, Samples: {len(sample_names)}")



if __name__ == "__main__":

    if len(sys.argv) != 3:

        print("Usage: python merge_telescope_counts.py <input_root> <output_tsv>")

        sys.exit(1)


    main(sys.argv[1], sys.argv[2])





scripts/build_locus_metadata.py


#!/usr/bin/env python3


import sys

import os

import gzip

import re

import pandas as pd



ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')



def open_maybe_gzip(path):

    if path.endswith(".gz"):

        return gzip.open(path, "rt")

    return open(path, "r")



def parse_attributes(attr_string):

    attrs = {}

    for key, value in ATTR_RE.findall(attr_string):

        attrs[key] = value

    return attrs



def family_from_attrs(attrs):

    # HERV annotation often has repFamily or category

    if "repFamily" in attrs and attrs["repFamily"]:

        return attrs["repFamily"]

    if "category" in attrs and attrs["category"]:

        return attrs["category"]

    if "family" in attrs and attrs["family"]:

        return attrs["family"]

    return "unknown"



def main(gtf_path, output_tsv):

    loci = {}


    with open_maybe_gzip(gtf_path) as fh:

        for line in fh:

            if not line.strip():

                continue

            if line.startswith("#"):

                continue


            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9:

                continue


            chrom, source, feature, start, end, score, strand, frame, attrs_raw = parts

            attrs = parse_attributes(attrs_raw)


            locus_id = attrs.get("locus")

            if not locus_id:

                locus_id = attrs.get("gene_id")

            if not locus_id:

                locus_id = attrs.get("transcript_id")

            if not locus_id:

                continue


            start = int(start)

            end = int(end)

            fam = family_from_attrs(attrs)


            if locus_id not in loci:

                loci[locus_id] = {

                    "locus_id": locus_id,

                    "family": fam,

                    "chromosome": chrom,

                    "start": start,

                    "end": end,

                    "strand": strand,

                }

            else:

                loci[locus_id]["start"] = min(loci[locus_id]["start"], start)

                loci[locus_id]["end"] = max(loci[locus_id]["end"], end)


                # prefer known family over unknown

                if loci[locus_id]["family"] == "unknown" and fam != "unknown":

                    loci[locus_id]["family"] = fam


    df = pd.DataFrame(list(loci.values()))

    df.sort_values(["chromosome", "start", "end", "locus_id"], inplace=True)


    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

    df.to_csv(output_tsv, sep="\t", index=False)


    print(f"Saved: {output_tsv}")

    print(f"Total loci: {df.shape[0]}")



if __name__ == "__main__":

    if len(sys.argv) != 3:

        print("Usage: python build_locus_metadata.py <input.gtf> <output.tsv>")

        sys.exit(1)


    main(sys.argv[1], sys.argv[2])





scripts/filter_counts.R



This is the strict version you wanted:


    keep locus only if all samples in both groups are >= min_count


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





scripts/run_deseq2_generic.R


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





scripts/extract_top_hits.R


args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 5) {

  stop("Usage: Rscript extract_top_hits.R <de_results.csv> <alpha> <log2fc_threshold> <top_up.csv> <top_down.csv>")

}


de_file <- args[1]

alpha <- as.numeric(args[2])

log2fc_threshold <- as.numeric(args[3])

up_file <- args[4]

down_file <- args[5]


df <- read.csv(de_file, check.names = FALSE)


required <- c("locus_id", "log2FoldChange", "padj")

missing <- setdiff(required, colnames(df))

if (length(missing) > 0) {

  stop(paste("Missing columns in DE results:", paste(missing, collapse = ", ")))

}


sig <- df[!is.na(df$padj) & df$padj < alpha & abs(df$log2FoldChange) >= log2fc_threshold, , drop = FALSE]


up <- sig[sig$log2FoldChange > 0, , drop = FALSE]

down <- sig[sig$log2FoldChange < 0, , drop = FALSE]


up <- up[order(up$padj, -up$log2FoldChange), , drop = FALSE]

down <- down[order(down$padj, down$log2FoldChange), , drop = FALSE]


top_up <- head(up, 10)

top_down <- head(down, 10)


write.csv(top_up, up_file, row.names = FALSE)

write.csv(top_down, down_file, row.names = FALSE)


cat("Saved:", up_file, "\n")

cat("Saved:", down_file, "\n")
