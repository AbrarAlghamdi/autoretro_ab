#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_DIR = PACKAGE_ROOT / "workflow"
ANNOTATIONS_DIR = PACKAGE_ROOT / "annotations"


def run_cmd(cmd, cwd=None):
    print("RUN:", " ".join(str(x) for x in cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def parse_sample_args(sample_args):
    samples = []
    for item in sample_args:
        if ":" not in item:
            raise ValueError(f"Invalid sample format: {item}. Use SRR_ID:condition")
        srr_id, condition = item.split(":", 1)
        srr_id = srr_id.strip()
        condition = condition.strip()
        if not srr_id or not condition:
            raise ValueError(f"Invalid sample format: {item}. Use SRR_ID:condition")
        samples.append((srr_id, condition))
    return samples


def write_samples_tsv(path: Path, samples):
    lines = ["sample_id\tcondition"]
    for srr_id, condition in samples:
        lines.append(f"{srr_id}\t{condition}")
    path.write_text("\n".join(lines) + "\n")


def write_download_json(path: Path, samples):
    srr_ids = [srr_id for srr_id, _ in samples]
    payload = {
        "ids": srr_ids,
        "output_dir": "results/fastq",
        "max_concurrent": 2
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def write_comparison_tsv(path: Path, samples):
    conditions = []
    for _, condition in samples:
        if condition not in conditions:
            conditions.append(condition)

    if len(conditions) >= 2:
        text = f"group1\tgroup2\n{conditions[1]}\t{conditions[0]}\n"
    else:
        text = "group1\tgroup2\ngroup1\tgroup2\n"

    path.write_text(text)


def init_project(project_dir: str, samples=None):
    project = Path(project_dir).resolve()
    config_dir = project / "config"

    (project / "config").mkdir(parents=True, exist_ok=True)
    (project / "results" / "fastq").mkdir(parents=True, exist_ok=True)
    (project / "results" / "bam").mkdir(parents=True, exist_ok=True)
    (project / "results" / "telescope_herv").mkdir(parents=True, exist_ok=True)
    (project / "results" / "telescope_l1").mkdir(parents=True, exist_ok=True)
    (project / "results" / "merged").mkdir(parents=True, exist_ok=True)
    (project / "results" / "metadata").mkdir(parents=True, exist_ok=True)
    (project / "results" / "deseq_herv").mkdir(parents=True, exist_ok=True)
    (project / "results" / "deseq_l1").mkdir(parents=True, exist_ok=True)
    (project / "logs").mkdir(parents=True, exist_ok=True)
    (project / "reference").mkdir(parents=True, exist_ok=True)

    # write SRR IDs file
    srr_file = config_dir / "srr_ids.txt"

    config_text = f"""project_root: {project}

samples_tsv: {project}/config/samples.tsv
comparison_tsv: {project}/config/comparison.tsv
srr_ids: {project}/config/srr_ids.txt

reference_fasta: resources/genome/hg38.fa
bowtie2_index: resources/index/hg38
herv_gtf: {ANNOTATIONS_DIR}/HERV_rmsk.hg38.v2.gtf
l1_gtf: {ANNOTATIONS_DIR}/transcripts_with_coords.gtf

threads: 8
min_count: 10
alpha: 0.05
log2fc_threshold: 1

bowtie2_bin: bowtie2
telescope_bin: telescope
"""
    (project / "config.yaml").write_text(config_text)

    if samples:
        write_samples_tsv(project / "config" / "samples.tsv", samples)
        write_comparison_tsv(project / "config" / "comparison.tsv", samples)

        with open(srr_file, "w") as f:
            for srr_id, condition in samples:
                f.write(f"{srr_id}\n")
    else:
        default_samples = [
            ("SRR1", "control"),
            ("SRR2", "treated"),
        ]
        write_samples_tsv(project / "config" / "samples.tsv", default_samples)
        write_comparison_tsv(project / "config" / "comparison.tsv", default_samples)

        with open(srr_file, "w") as f:
            for srr_id, condition in default_samples:
                f.write(f"{srr_id}\n")

    print(f"Project initialized at: {project}")
    if samples:
        print("Samples and conditions added automatically.")
        print("Now edit only:")
        print(f"  {project}/config.yaml")
    else:
        print("Edit config.yaml and files in config/ before running.")

threads: 8
min_count: 10
alpha: 0.05
log2fc_threshold: 1

bowtie2_bin: bowtie2
telescope_bin: telescope
"""
    (project / "config.yaml").write_text(config_text)

    if samples:
        write_samples_tsv(project / "config" / "samples.tsv", samples)
        write_download_json(project / "config" / "download.json", samples)
        write_comparison_tsv(project / "config" / "comparison.tsv", samples)
    else:
        default_samples = [
            ("SRR1", "control"),
            ("SRR2", "treated"),
        ]
        write_samples_tsv(project / "config" / "samples.tsv", default_samples)
        write_download_json(project / "config" / "download.json", default_samples)
        write_comparison_tsv(project / "config" / "comparison.tsv", default_samples)

    print(f"Project initialized at: {project}")
    if samples:
        print("Samples and conditions added automatically.")
        print("Now edit only:")
        print(f"  {project}/config.yaml")
    else:
        print("Edit config.yaml and files in config/ before running.")


def run_workflow(project_dir: str, cores: int):
    import sys

    project = Path(project_dir).resolve()

    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--directory", str(project),
        "--configfile", str(project / "config.yaml"),
        "--cores", str(cores),
    ]
    run_cmd(cmd)

def main():
    parser = argparse.ArgumentParser(prog="autoretro_ab")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_init = subparsers.add_parser("init", help="Initialize a new autoretro_ab project")
    p_init.add_argument("project_dir", help="Path to the new project directory")
    p_init.add_argument(
        "--sample",
        nargs="+",
        help="Samples as SRR_ID:condition, e.g. SRR1553606:control SRR1553607:treated"
    )

    p_run = subparsers.add_parser("run", help="Run the pipeline with Snakemake")
    p_run.add_argument("project_dir", help="Path to the project directory")
    p_run.add_argument("--cores", type=int, default=4, help="Number of CPU cores")

    args = parser.parse_args()

    if args.command == "init":
        samples = parse_sample_args(args.sample) if args.sample else None
        init_project(args.project_dir, samples=samples)
    elif args.command == "run":
        run_workflow(args.project_dir, args.cores)


if __name__ == "__main__":
    main()
