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


