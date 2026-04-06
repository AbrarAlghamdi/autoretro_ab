#!/usr/bin/env python3

import os
import sys
import subprocess


def run_cmd(cmd, log_file=None):
    print("RUN:", " ".join(cmd), flush=True)

    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    ) as proc:
        output_lines = []
        for line in proc.stdout:
            print(line, end="")
            output_lines.append(line)
        ret = proc.wait()

    if log_file:
        with open(log_file, "a") as f:
            f.write("\n$ " + " ".join(cmd) + "\n")
            f.writelines(output_lines)
            f.write(f"\n[exit_code={ret}]\n")

    if ret != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def parse_key_value_file(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def parse_yaml_like_config(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                continue
            k, v = line.split(":", 1)
            d[k.strip()] = v.strip()
    return d


def find_fastqs(fastq_dir, sample_id):
    candidates = {
        "single_gz": os.path.join(fastq_dir, f"{sample_id}.fastq.gz"),
        "single": os.path.join(fastq_dir, f"{sample_id}.fastq"),
        "pair1_gz": os.path.join(fastq_dir, f"{sample_id}_1.fastq.gz"),
        "pair2_gz": os.path.join(fastq_dir, f"{sample_id}_2.fastq.gz"),
        "pair1": os.path.join(fastq_dir, f"{sample_id}_1.fastq"),
        "pair2": os.path.join(fastq_dir, f"{sample_id}_2.fastq"),
    }

    if os.path.exists(candidates["pair1_gz"]) and os.path.exists(candidates["pair2_gz"]):
        return "paired", [candidates["pair1_gz"], candidates["pair2_gz"]]

    if os.path.exists(candidates["pair1"]) and os.path.exists(candidates["pair2"]):
        return "paired", [candidates["pair1"], candidates["pair2"]]

    if os.path.exists(candidates["single_gz"]):
        return "single", [candidates["single_gz"]]

    if os.path.exists(candidates["single"]):
        return "single", [candidates["single"]]

    raise FileNotFoundError(f"No FASTQ found for {sample_id}")


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def read_samples(samples_file):
    samples = []
    with open(samples_file) as f:
        header = f.readline().strip().split("\t")
        if "sample_id" not in header:
            raise ValueError("samples.tsv must contain a 'sample_id' column")
        i_sample = header.index("sample_id")

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            samples.append(parts[i_sample])

    return samples


def align_sample(sample_id, fastq_dir, bam_dir, bowtie2_index, bowtie2_bin, threads, log_file):
    mode, fq_files = find_fastqs(fastq_dir, sample_id)

    sample_bam = os.path.join(bam_dir, f"{sample_id}.sorted.bam")

    if os.path.exists(sample_bam):
        print(f"Skipping alignment for {sample_id}: sorted BAM already exists", flush=True)
        return sample_bam

    sam_path = os.path.join(bam_dir, f"{sample_id}.sam")
    bam_path = os.path.join(bam_dir, f"{sample_id}.bam")

    if mode == "paired":
        cmd = [
            bowtie2_bin,
            "-x", bowtie2_index,
            "-1", fq_files[0],
            "-2", fq_files[1],
            "-k", "100",
            "-p", str(threads),
            "-S", sam_path,
        ]
    else:
        cmd = [
            bowtie2_bin,
            "-x", bowtie2_index,
            "-U", fq_files[0],
            "-k", "100",
            "-p", str(threads),
            "-S", sam_path,
        ]

    run_cmd(cmd, log_file=log_file)
    run_cmd(["samtools", "view", "-@", str(threads), "-bS", sam_path, "-o", bam_path], log_file=log_file)
    run_cmd(["samtools", "sort", "-@", str(threads), "-o", sample_bam, bam_path], log_file=log_file)
    run_cmd(["samtools", "index", sample_bam], log_file=log_file)

    if os.path.exists(sam_path):
        os.remove(sam_path)
    if os.path.exists(bam_path):
        os.remove(bam_path)

    return sample_bam


def run_telescope(sample_id, bam_file, gtf_file, out_root, telescope_bin, threads, log_file):
    sample_out = os.path.join(out_root, sample_id)
    ensure_dir(sample_out)

    report1 = os.path.join(sample_out, "telescope_report.tsv")
    report2 = os.path.join(sample_out, "telescope-telescope_report.tsv")

    if os.path.exists(report1) or os.path.exists(report2):
        print(f"Skipping Telescope for {sample_id}: report already exists", flush=True)
        return

    cmd = [
        telescope_bin,
        "assign",
        bam_file,
        gtf_file,
        "--outdir", sample_out,
        "--exp_tag", "telescope",
        "--theta_prior", "200000",
        "--ncpu", str(threads),
    ]

    run_cmd(cmd, log_file=log_file)


def main():
    if len(sys.argv) != 2:
        print("Usage: python run_alignment_and_telescope.py <project_root>")
        sys.exit(1)

    project_root = os.path.abspath(sys.argv[1])

    config_yaml = os.path.join(project_root, "config.yaml")
    config_dir = os.path.join(project_root, "config")
    results_dir = os.path.join(project_root, "results")
    logs_dir = os.path.join(project_root, "logs")

    cfg = parse_yaml_like_config(config_yaml)
    samples = read_samples(os.path.join(config_dir, "samples.tsv"))

    bowtie2_index = cfg["bowtie2_index"]
    herv_gtf = cfg["herv_gtf"]
    l1_gtf = cfg["l1_gtf"]
    bowtie2_bin = cfg.get("bowtie2_bin", "bowtie2")
    telescope_bin = cfg.get("telescope_bin", "telescope")
    threads = int(cfg.get("threads", 8))

    fastq_dir = os.path.join(results_dir, "fastq")
    bam_dir = os.path.join(results_dir, "bam")
    herv_dir = os.path.join(results_dir, "telescope_herv")
    l1_dir = os.path.join(results_dir, "telescope_l1")

    ensure_dir(bam_dir)
    ensure_dir(herv_dir)
    ensure_dir(l1_dir)
    ensure_dir(logs_dir)

    log_file = os.path.join(logs_dir, "pipeline.log")

    if not os.path.exists(bowtie2_bin) and bowtie2_bin != "bowtie2":
        raise RuntimeError(f"Bowtie2 binary not found: {bowtie2_bin}")
    if not os.path.exists(telescope_bin) and telescope_bin != "telescope":
        raise RuntimeError(f"Telescope binary not found: {telescope_bin}")

    for sample_id in samples:
        print(f"\n=== Processing {sample_id} ===", flush=True)

        bam = align_sample(
            sample_id=sample_id,
            fastq_dir=fastq_dir,
            bam_dir=bam_dir,
            bowtie2_index=bowtie2_index,
            bowtie2_bin=bowtie2_bin,
            threads=threads,
            log_file=log_file
        )

        run_telescope(
            sample_id=sample_id,
            bam_file=bam,
            gtf_file=herv_gtf,
            out_root=herv_dir,
            telescope_bin=telescope_bin,
            threads=threads,
            log_file=log_file
        )

        run_telescope(
            sample_id=sample_id,
            bam_file=bam,
            gtf_file=l1_gtf,
            out_root=l1_dir,
            telescope_bin=telescope_bin,
            threads=threads,
            log_file=log_file
        )

    print("\nDONE: alignment + Telescope complete", flush=True)


if __name__ == "__main__":
    main()
