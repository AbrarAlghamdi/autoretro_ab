
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


