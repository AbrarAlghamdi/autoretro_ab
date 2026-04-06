#!/usr/bin/env python3

import aiohttp
import aiofiles
import asyncio
import os
import time
import argparse
import shutil
from aiohttp import TCPConnector
from asyncio import Semaphore

# -----------------------------
# GEO download settings
# -----------------------------
BASE_HTTP_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{0}nnn/{1}/{2}/"

FILE_TYPES = {
    "matrix": {"suffix": "series_matrix.txt.gz", "folder": "matrix"},
    "miniml": {"suffix": "family.xml.tgz", "folder": "miniml"},
    "soft": {"suffix": "family.soft.gz", "folder": "soft"}
}


# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def detect_tool(preferred=None, fallback_name=None):
    if preferred:
        return preferred
    if fallback_name:
        return shutil.which(fallback_name)
    return None


async def run_cmd(cmd):
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()
    return proc.returncode, stdout.decode(errors="ignore"), stderr.decode(errors="ignore")


async def compress_fastqs(output_dir, srr_id, pigz_bin=None, gzip_bin=None):
    candidates = [
        os.path.join(output_dir, f"{srr_id}.fastq"),
        os.path.join(output_dir, f"{srr_id}_1.fastq"),
        os.path.join(output_dir, f"{srr_id}_2.fastq"),
    ]

    tasks = []
    for fq in candidates:
        if os.path.exists(fq):
            if pigz_bin:
                cmd = [pigz_bin, "-f", fq]
            elif gzip_bin:
                cmd = [gzip_bin, "-f", fq]
            else:
                raise RuntimeError("Neither pigz nor gzip found for FASTQ compression")
            tasks.append(run_cmd(cmd))

    if tasks:
        results = await asyncio.gather(*tasks)
        for code, out, err in results:
            if code != 0:
                raise RuntimeError(f"Compression failed:\n{err}")


def fastq_outputs_exist(output_dir, srr_id):
    candidates = [
        os.path.join(output_dir, f"{srr_id}.fastq.gz"),
        os.path.join(output_dir, f"{srr_id}_1.fastq.gz"),
        os.path.join(output_dir, f"{srr_id}_2.fastq.gz"),
        os.path.join(output_dir, f"{srr_id}.fastq"),
        os.path.join(output_dir, f"{srr_id}_1.fastq"),
        os.path.join(output_dir, f"{srr_id}_2.fastq"),
    ]
    return any(os.path.exists(x) for x in candidates)


def find_local_sra(cache_dir, srr_id):
    candidates = [
        os.path.join(cache_dir, srr_id, f"{srr_id}.sra"),
        os.path.join(cache_dir, f"{srr_id}.sra"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


# -----------------------------
# GEO download
# -----------------------------
async def download_gse_file(session, semaphore, gse_id, output_dir, file_type):
    gse_prefix = gse_id[:-3]
    file_config = FILE_TYPES[file_type]
    url = BASE_HTTP_URL.format(gse_prefix, gse_id, file_config["folder"]) + f"{gse_id}_{file_config['suffix']}"
    local_filename = f"{gse_id}_{file_config['suffix']}"
    local_file = os.path.join(output_dir, local_filename)

    print(f"Starting GEO download: {local_filename}", flush=True)

    async with semaphore:
        try:
            async with session.get(url) as response:
                if response.status == 200:
                    async with aiofiles.open(local_file, "wb") as f:
                        await f.write(await response.read())
                    print(f"Downloaded {local_filename}", flush=True)
                    return local_file
                else:
                    print(f"Failed {local_filename} (HTTP {response.status})", flush=True)
                    return None
        except Exception as e:
            print(f"Error downloading {local_filename}: {e}", flush=True)
            return None


# -----------------------------
# SRA download
# -----------------------------
async def download_fastq_file(
    semaphore,
    srr_id,
    output_dir,
    sra_cache_dir,
    fasterq_dump_bin=None,
    prefetch_bin=None,
    pigz_bin=None,
    gzip_bin=None,
    threads=4,
    max_spots=None,
):
    print(f"Starting SRR download: {srr_id}", flush=True)

    async with semaphore:
        try:
            if fastq_outputs_exist(output_dir, srr_id):
                print(f"Skipping {srr_id}: FASTQ already exists", flush=True)
                return True

            ensure_dir(output_dir)
            ensure_dir(sra_cache_dir)

            fasterq = detect_tool(fasterq_dump_bin, "fasterq-dump")
            prefetch = detect_tool(prefetch_bin, "prefetch")
            pigz = detect_tool(pigz_bin, "pigz")
            gzip_tool = detect_tool(gzip_bin, "gzip")

            if not fasterq:
                raise RuntimeError("fasterq-dump not found")
            if not prefetch:
                raise RuntimeError("prefetch not found")
            if not pigz and not gzip_tool:
                raise RuntimeError("Neither pigz nor gzip found")

            # -----------------------------
            # Try direct fasterq-dump first
            # -----------------------------
            direct_cmd = [
                fasterq,
                srr_id,
                "-O", output_dir,
                "--split-files",
                "-e", str(threads),
            ]
            if max_spots is not None:
                direct_cmd.extend(["-X", str(max_spots)])

            code, out, err = await run_cmd(direct_cmd)
            if code == 0:
                await compress_fastqs(output_dir, srr_id, pigz, gzip_tool)
                print(f"Success: {srr_id} (direct fasterq-dump)", flush=True)
                return True

            print(f"Direct fasterq-dump failed for {srr_id}, trying prefetch fallback...", flush=True)
            if err:
                print(err, flush=True)

            # -----------------------------
            # Fallback: prefetch + local fasterq-dump
            # -----------------------------
            prefetch_cmd = [prefetch, srr_id, "-O", sra_cache_dir]
            code2, out2, err2 = await run_cmd(prefetch_cmd)
            if code2 != 0:
                print(f"prefetch failed for {srr_id}", flush=True)
                if err2:
                    print(err2, flush=True)
                return False

            sra_path = find_local_sra(sra_cache_dir, srr_id)
            if not sra_path:
                print(f"Could not find local .sra file for {srr_id}", flush=True)
                return False

            fallback_cmd = [
                fasterq,
                sra_path,
                "-O", output_dir,
                "--split-files",
                "-e", str(threads),
            ]
            if max_spots is not None:
                fallback_cmd.extend(["-X", str(max_spots)])

            code3, out3, err3 = await run_cmd(fallback_cmd)
            if code3 != 0:
                print(f"Fallback fasterq-dump failed for {srr_id}", flush=True)
                if err3:
                    print(err3, flush=True)
                return False

            await compress_fastqs(output_dir, srr_id, pigz, gzip_tool)
            print(f"Success: {srr_id} (prefetch fallback)", flush=True)
            return True

        except Exception as e:
            print(f"Error downloading {srr_id}: {e}", flush=True)
            return False


# -----------------------------
# Main workflow
# -----------------------------
async def download_all_files(
    file_ids,
    output_dir,
    max_concurrent=4,
    db_type="sra",
    sra_cache_dir="sra_cache",
    fasterq_dump_bin=None,
    prefetch_bin=None,
    pigz_bin=None,
    gzip_bin=None,
    threads=4,
    max_spots=None,
):
    start_time = time.perf_counter()
    ensure_dir(output_dir)

    semaphore = Semaphore(max_concurrent)
    failed_ids = []

    print(f"Starting downloads for {len(file_ids)} IDs ({db_type.upper()})", flush=True)

    if db_type == "geo":
        connector = TCPConnector(limit=max_concurrent)
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = [
                download_gse_file(session, semaphore, gse_id, output_dir, file_type)
                for gse_id in file_ids
                for file_type in FILE_TYPES
            ]
            results = await asyncio.gather(*tasks)
            for i, result in enumerate(results):
                if result is None:
                    failed_ids.append(file_ids[i // len(FILE_TYPES)])

    else:
        tasks = [
            download_fastq_file(
                semaphore=semaphore,
                srr_id=srr_id,
                output_dir=output_dir,
                sra_cache_dir=sra_cache_dir,
                fasterq_dump_bin=fasterq_dump_bin,
                prefetch_bin=prefetch_bin,
                pigz_bin=pigz_bin,
                gzip_bin=gzip_bin,
                threads=threads,
                max_spots=max_spots,
            )
            for srr_id in file_ids
        ]
        results = await asyncio.gather(*tasks)
        failed_ids = [file_ids[i] for i, ok in enumerate(results) if not ok]

    end_time = time.perf_counter()
    total_time = end_time - start_time
    minutes = int(total_time // 60)
    seconds = total_time % 60

    print(f"All downloads finished in {minutes} min {seconds:.2f} sec", flush=True)

    return failed_ids


def main(
    file_ids,
    output_dir="Downloads",
    max_concurrent=4,
    db_type="sra",
    sra_cache_dir="sra_cache",
    fasterq_dump_bin=None,
    prefetch_bin=None,
    pigz_bin=None,
    gzip_bin=None,
    threads=4,
    max_spots=None,
):
    failed_ids = asyncio.run(
        download_all_files(
            file_ids=file_ids,
            output_dir=output_dir,
            max_concurrent=max_concurrent,
            db_type=db_type,
            sra_cache_dir=sra_cache_dir,
            fasterq_dump_bin=fasterq_dump_bin,
            prefetch_bin=prefetch_bin,
            pigz_bin=pigz_bin,
            gzip_bin=gzip_bin,
            threads=threads,
            max_spots=max_spots,
        )
    )
    return failed_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GEO or SRA data")
    parser.add_argument("-i", "--input_path", type=str, required=True, help="Path to input file with IDs")
    parser.add_argument("-d", "--database_name", choices=["GEO", "SRA", "geo", "sra"], required=True, help="Database name")
    parser.add_argument("-o", "--output_dir", type=str, default="Downloads", help="Output directory")
    parser.add_argument("--max_concurrent", type=int, default=4, help="Maximum concurrent downloads")
    parser.add_argument("--sra_cache_dir", type=str, default="sra_cache", help="SRA cache directory")
    parser.add_argument("--threads", type=int, default=4, help="Threads for fasterq-dump")
    parser.add_argument("--max_spots", type=int, default=None, help="Optional read limit for testing")
    parser.add_argument("--fasterq_dump_bin", type=str, default=None, help="Path to fasterq-dump")
    parser.add_argument("--prefetch_bin", type=str, default=None, help="Path to prefetch")
    parser.add_argument("--pigz_bin", type=str, default=None, help="Path to pigz")
    parser.add_argument("--gzip_bin", type=str, default=None, help="Path to gzip")

    args = parser.parse_args()

    with open(args.input_path, "r") as ifile:
        ncbi_ids = ifile.read().strip().splitlines()

    gse_ids, srr_ids = [], []
    for Id in ncbi_ids:
        if Id.startswith("GSE"):
            gse_ids.append(Id)
        elif Id.startswith("SRR"):
            srr_ids.append(Id)
        else:
            print(f"Invalid ID: {Id}", flush=True)

    failed_ids = []
    if args.database_name.lower() == "geo":
        print("GSE IDs:", ", ".join(gse_ids), flush=True)
        failed_ids = main(
            gse_ids,
            output_dir=args.output_dir,
            max_concurrent=args.max_concurrent,
            db_type="geo",
        )
    else:
        print("SRR IDs:", ", ".join(srr_ids), flush=True)
        failed_ids = main(
            srr_ids,
            output_dir=args.output_dir,
            max_concurrent=args.max_concurrent,
            db_type="sra",
            sra_cache_dir=args.sra_cache_dir,
            fasterq_dump_bin=args.fasterq_dump_bin,
            prefetch_bin=args.prefetch_bin,
            pigz_bin=args.pigz_bin,
            gzip_bin=args.gzip_bin,
            threads=args.threads,
            max_spots=args.max_spots,
        )

    if failed_ids:
        print("\nFailed IDs:", flush=True)
        for failed_id in sorted(set(failed_ids)):
            print(failed_id, flush=True)
    else:
        print("\nAll IDs downloaded successfully.", flush=True)
