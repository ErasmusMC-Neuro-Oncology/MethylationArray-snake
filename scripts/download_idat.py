#!/usr/bin/env python3

import csv
import gzip
import html
import os
import shutil
import time
import urllib.request

samplesheet_path = snakemake.input[0]
out_dir = snakemake.output[0]
out_samplesheet = os.path.join(out_dir, "samplesheet.csv")

os.makedirs(out_dir, exist_ok=True)

with open(samplesheet_path, newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    fieldnames = reader.fieldnames


def download_and_extract(url, dest_path):
    url = html.unescape(url)
    tmp_path = dest_path + ".gz"
    try:
        urllib.request.urlretrieve(url, tmp_path)
    except (urllib.error.URLError, urllib.error.ContentTooShortError):
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    time.sleep(1)
    with gzip.open(tmp_path, "rb") as f_in, open(dest_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_path)


out_fieldnames = list(fieldnames)
if "sample" not in out_fieldnames:
    out_fieldnames.append("sample")

for row in rows:
    red_filename = os.path.basename(row["idat_red"]).replace(".gz", "")
    red_local = os.path.join(out_dir, red_filename)
    if not os.path.exists(red_local):
        print(f"Downloading: {red_filename}", flush=True)
        download_and_extract(row["url_red"], red_local)
    row["idat_red"] = red_local

    grn_filename = red_filename.replace("_Red.idat", "_Grn.idat")
    grn_local = os.path.join(out_dir, grn_filename)
    if not os.path.exists(grn_local):
        print(f"Downloading: {grn_filename}", flush=True)
        download_and_extract(row["url_grn"], grn_local)

    if "sample" not in fieldnames:
        row["sample"] = row["idat_basename"]

with open(out_samplesheet, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_fieldnames)
    writer.writeheader()
    writer.writerows(rows)
