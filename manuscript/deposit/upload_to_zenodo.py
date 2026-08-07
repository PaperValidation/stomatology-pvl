#!/usr/bin/env python3
"""
Publish the manuscript deposit package to Zenodo via the REST API.

Usage:
    export ZENODO_TOKEN=<your_personal_access_token>
    python3 upload_to_zenodo.py --sandbox     # dry-run against Zenodo sandbox
    python3 upload_to_zenodo.py               # publish to production (real DOI)

The script creates a new deposition, uploads every file in this directory
(except itself and hidden files), sets metadata, and publishes.
It prints the resulting DOI / record URL on success.
"""
import os
import sys
import json
import argparse
import urllib.request
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
SANDBOX_API = "https://sandbox.zenodo.org/api"
PROD_API = "https://zenodo.org/api"


def api_request(url, method, token, data=None, binary=None):
    headers = {"Authorization": f"Bearer {token}"}
    if binary is None:
        headers["Content-Type"] = "application/json"
        payload = json.dumps(data).encode() if data is not None else None
    else:
        payload = binary
    req = urllib.request.Request(url, data=payload, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as r:
            body = r.read().decode()
            return r.status, (json.loads(body) if body else {})
    except urllib.error.HTTPError as e:
        sys.stderr.write(f"HTTP {e.code}: {e.read().decode()}\n")
        raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use Zenodo sandbox (dry-run)")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ERROR: set ZENODO_TOKEN environment variable (personal access token).")

    base = SANDBOX_API if args.sandbox else PROD_API

    # 1) create deposition
    status, dep = api_request(f"{base}/deposit/depositions", "POST", token, data={})
    if status >= 400:
        sys.exit(f"Failed to create deposition (HTTP {status}).")
    dep_id = dep["id"]
    print(f"Created deposition id={dep_id}")

    # 2) upload files
    for fn in sorted(os.listdir(HERE)):
        fp = os.path.join(HERE, fn)
        if not os.path.isfile(fp):
            continue
        if fn.startswith(".") or fn == "upload_to_zenodo.py":
            continue
        with open(fp, "rb") as f:
            data = f.read()
        bucket = dep.get("links", {}).get("bucket")
        if not bucket:
            # fallback to files API
            url = f"{base}/deposit/depositions/{dep_id}/files"
            status, _ = api_request(url, "POST", token,
                                    data={"filename": fn},
                                    binary=data)
        else:
            url = f"{bucket}/{fn}"
            status, _ = api_request(url, "PUT", token, binary=data)
        print(f"  uploaded {fn}: HTTP {status}")

    # 3) set metadata
    metadata = {
        "metadata": {
            "title": ("A validated transcriptomic signature of malignant progression in "
                      "oral leukoplakia and proliferative verrucous leukoplakia: "
                      "data and analysis code"),
            "upload_type": "dataset",
            "description": ("Data objects and analysis code reproducing the 101-gene "
                            "progression signature, per-sample direction-weighted scores "
                            "for validation cohorts GSE227919 and GSE85195, and single-cell "
                            "cluster statistics. Re-analysis of public GEO datasets "
                            "GSE26549, GSE227919, GSE196296, GSE85195."),
            "creators": [{"name": "[Author Name]", "affiliation": "[Institution]"}],
            "access_right": "open",
            "license": "cc-by-4.0",
            "keywords": ["oral leukoplakia", "proliferative verrucous leukoplakia",
                         "transcriptomic signature", "single-cell RNA-seq", "GEO"],
            "communities": [{"identifier": "oral-oncology"}],
        }
    }
    status, dep = api_request(f"{base}/deposit/depositions/{dep_id}",
                              "PUT", token, data=metadata)
    print(f"Metadata set: HTTP {status}")

    # 4) publish
    status, dep = api_request(f"{base}/deposit/depositions/{dep_id}/actions/publish",
                              "POST", token, data={})
    if status >= 400:
        sys.exit(f"Publish failed (HTTP {status}). Record draft: {dep.get('links', {}).get('html')}")
    doi = dep.get("doi", "n/a")
    print(f"PUBLISHED. DOI: https://doi.org/{doi}")
    print(f"Record: {dep.get('links', {}).get('html')}")


if __name__ == "__main__":
    main()
