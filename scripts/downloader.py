
#!/usr/bin/env python3

import asyncio
import json
import os
import shutil
import sys
import time
from pathlib import Path


def detect_sra_download_tool():
    """
    Prefer fasterq-dump if available, otherwise fall back to prefetch + fasterq-dump if needed.
    """
    fasterq = shutil.which("fasterq-dump")
    prefetch = shutil.which("prefetch")
    pigz = shutil.which("pigz")
    gzip_bin = shutil.which("gzip")

    return {
        "fasterq_dump": fasterq,
        "prefetch": prefetch,
        "pigz": pigz,
        "gzip": gzip_bin,
    }


def compress_fastq_files(output_dir, srr_id, pigz_bin=None, gzip_bin=None):
    """
    Compress any FASTQ files produced for this SRR.
    """
    patterns = [
        f"{srr_id}.fastq",
        f"{srr_id}_1.fastq",
        f"{srr_id}_2.fastq",
    ]

    for name in patterns:
        fq = os.path.join(output_dir, name)
        if os.path.exists(fq):
            if pigz_bin:
                cmd = [pigz_bin, "-f", fq]
            elif gzip_bin:
                cmd = [gzip_bin, "-f", fq]
            else:
                raise RuntimeError("Neither pigz nor gzip found for compression")

            proc = asyncio.create_subprocess_exec(*cmd)
            # This function is sync, so compression is handled outside async
            raise RuntimeError(
                "compress_fastq_files should not be called directly in sync mode"
            )


async def compress_one_file(fq, pigz_bin=None, gzip_bin=None):
    if pigz_bin:
        cmd = [pigz_bin, "-f", fq]
    elif gzip_bin:
        cmd = [gzip_bin, "-f", fq]
    else:
        raise RuntimeError("Neither pigz nor gzip found for compression")

    proc = await asyncio.create_subprocess_exec(*cmd)
    await proc.communicate()

    if proc.returncode != 0:
        raise RuntimeError(f"Compression failed for {fq}")


async def compress_fastqs_for_srr(output_dir, srr_id, pigz_bin=None, gzip_bin=None):
    patterns = [
        f"{srr_id}.fastq",
        f"{srr_id}_1.fastq",
        f"{srr_id}_2.fastq",
    ]

    tasks = []
    for name in patterns:
        fq = os.path.join(output_dir, name)
        if os.path.exists(fq):
            tasks.append(compress_one_file(fq, pigz_bin, gzip_bin))

    if tasks:
        await asyncio.gather(*tasks)


async def run_cmd(cmd, cwd=None):
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=cwd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()
    return proc.returncode, stdout.decode(errors="ignore"), stderr.decode(errors="ignore")


async def download_one_srr(srr_id, output_dir, semaphore, tools):
    async with semaphore:
        print(f"Starting download for {srr_id}", flush=True)

        fasterq = tools["fasterq_dump"]
        prefetch = tools["prefetch"]
        pigz_bin = tools["pigz"]
        gzip_bin = tools["gzip"]

        if fasterq is None:
            raise RuntimeError("fasterq-dump not found in PATH")

        # skip if already present
        already = [
            os.path.join(output_dir, f"{srr_id}.fastq.gz"),
            os.path.join(output_dir, f"{srr_id}_1.fastq.gz"),
            os.path.join(output_dir, f"{srr_id}_2.fastq.gz"),
        ]
        if any(os.path.exists(x) for x in already):
            print(f"Skipping {srr_id}: already downloaded", flush=True)
            return

        # direct fasterq-dump
        cmd = [
            fasterq,
            srr_id,
            "-O", output_dir,
            "--split-files",
            "-e", "4",
        ]

        code, out, err = await run_cmd(cmd)

        if code != 0:
            # fallback: prefetch + fasterq-dump
            if prefetch is None:
                print(f"Failed: {srr_id}", flush=True)
                print(err, file=sys.stderr)
                raise RuntimeError(f"Download failed for {srr_id}")

            tmp_sra_dir = os.path.join(output_dir, "_prefetch_cache")
            os.makedirs(tmp_sra_dir, exist_ok=True)

            prefetch_cmd = [
                prefetch,
                srr_id,
                "-O", tmp_sra_dir,
            ]
            code2, out2, err2 = await run_cmd(prefetch_cmd)

            if code2 != 0:
                print(f"Failed: {srr_id}", flush=True)
                print(err2, file=sys.stderr)
                raise RuntimeError(f"prefetch failed for {srr_id}")

            sra_path = os.path.join(tmp_sra_dir, srr_id, f"{srr_id}.sra")
            if not os.path.exists(sra_path):
                # some versions place it differently
                sra_path = os.path.join(tmp_sra_dir, f"{srr_id}.sra")

            fasterq_cmd = [
                fasterq,
                sra_path,
                "-O", output_dir,
                "--split-files",
                "-e", "4",
            ]
            code3, out3, err3 = await run_cmd(fasterq_cmd)

            if code3 != 0:
                print(f"Failed: {srr_id}", flush=True)
                print(err3, file=sys.stderr)
                raise RuntimeError(f"fasterq-dump fallback failed for {srr_id}")

        # compress FASTQ files
        await compress_fastqs_for_srr(output_dir, srr_id, pigz_bin, gzip_bin)

        print(f"Success: {srr_id}", flush=True)


async def main_async(config):
    ids = config.get("ids", [])
    output_dir = config.get("output_dir")
    max_concurrent = int(config.get("max_concurrent", 4))

    if not ids:
        raise ValueError("No SRR IDs found in input JSON under 'ids'")
    if not output_dir:
        raise ValueError("No output_dir found in input JSON")

    os.makedirs(output_dir, exist_ok=True)

    tools = detect_sra_download_tool()

    if tools["fasterq_dump"] is None:
        raise RuntimeError(
            "fasterq-dump not found. Please install sra-tools in your environment."
        )
    if tools["pigz"] is None and tools["gzip"] is None:
        raise RuntimeError(
            "Neither pigz nor gzip found. Please install one for FASTQ compression."
        )

    semaphore = asyncio.Semaphore(max_concurrent)

    start = time.time()
    tasks = [download_one_srr(srr, output_dir, semaphore, tools) for srr in ids]
    await asyncio.gather(*tasks)
    end = time.time()

    print(f"\nFinished in {end - start:.2f} seconds", flush=True)
    print("All downloads successful", flush=True)


def main(json_path):
    with open(json_path) as f:
        config = json.load(f)

    asyncio.run(main_async(config))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python downloader.py <input.json>")
        sys.exit(1)

    main(sys.argv[1])
