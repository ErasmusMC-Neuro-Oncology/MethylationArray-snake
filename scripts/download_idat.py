#!/usr/bin/env python3

import csv
import gzip
import html
import os
import shutil
import sys
import time
import urllib.request

samplesheet_path = snakemake.input[0]
out_dir = snakemake.output[0]
out_samplesheet = os.path.join(out_dir, "samplesheet.csv")

os.makedirs(out_dir, exist_ok=True)

with open(samplesheet_path, newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    fieldnames = [fn for fn in reader.fieldnames if fn]  # Filter None fieldnames
    # Clean rows: remove None keys
    rows = [{k: v for k, v in row.items() if k} for row in rows]


def download_and_extract(url, dest_path, max_retries=3):
    url = html.unescape(url)
    is_gzipped = url.endswith(".gz")
    tmp_path = dest_path + ".gz" if is_gzipped else dest_path + ".tmp"

    for attempt in range(max_retries):
        try:
            urllib.request.urlretrieve(url, tmp_path)
            break
        except (urllib.error.URLError, urllib.error.ContentTooShortError) as e:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            if attempt < max_retries - 1:
                wait_time = 5 * (2 ** attempt)
                print(f"Download failed, retrying in {wait_time}s (attempt {attempt + 1}/{max_retries})", flush=True)
                time.sleep(wait_time)
            else:
                raise

    time.sleep(1)

    if is_gzipped:
        with gzip.open(tmp_path, "rb") as f_in, open(dest_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(tmp_path)
    else:
        shutil.move(tmp_path, dest_path)


out_fieldnames = list(fieldnames)
if "sample" not in out_fieldnames:
    out_fieldnames.append("sample")
if "idat_grn" not in out_fieldnames:
    out_fieldnames.append("idat_grn")

completed_rows = []
failed_samples = []

for i, row in enumerate(rows):
    sample_id = row.get('idat_basename', f'sample_{i}')
    try:
        red_filename = os.path.basename(row["idat_red"]).replace(".gz", "")

        # Download both channels
        for channel in ["Red", "Grn"]:
            filename = red_filename if channel == "Red" else red_filename.replace("_Red.idat", "_Grn.idat")
            local_path = os.path.join(out_dir, filename)

            if os.path.exists(local_path):
                print(f"Skipping (exists): {filename}", flush=True)
            else:
                print(f"Downloading: {filename}", flush=True)
                download_and_extract(row[f"url_{channel.lower()}"], local_path)

            row[f"idat_{channel.lower()}"] = local_path

        # Only mark as complete if both channels are present
        if os.path.exists(row["idat_red"]) and os.path.exists(row["idat_grn"]):
            if "sample" not in fieldnames:
                row["sample"] = row["idat_basename"]
            completed_rows.append(row)
        else:
            raise RuntimeError(f"Missing idat files for sample {sample_id}")

    except Exception as e:
        print(f"WARNING: Failed to download sample {sample_id} (row {i}): {e}", file=sys.stderr, flush=True)
        failed_samples.append(sample_id)

if failed_samples:
    print(f"Failed to download {len(failed_samples)} sample(s): {', '.join(failed_samples)}", file=sys.stderr, flush=True)
    print(f"Successfully downloaded {len(completed_rows)} sample(s), continuing with available data", file=sys.stderr, flush=True)

with open(out_samplesheet, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_fieldnames)
    writer.writeheader()
    writer.writerows(completed_rows)

# Write marker file to indicate successful completion (even with partial downloads)
with open(os.path.join(out_dir, ".download_complete"), "w") as f:
    f.write(f"Completed: {len(completed_rows)} samples\n")
    if failed_samples:
        f.write(f"Failed: {len(failed_samples)} samples: {', '.join(failed_samples)}\n")
