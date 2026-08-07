#!/usr/bin/env python3
"""
Publish the existing Zenodo DRAFT deposition (default id 21823584) via the Zenodo
REST API (v2). Reads ZENODO_TOKEN from the environment. Routes through the local
SOCKS5 proxy at 127.0.0.1:8428 (the only egress that reaches zenodo.org in this env).

Logic:
  - GET the draft. If it already has a title and files -> just PUBLISH.
  - If it is empty (reserved DOI only) -> upload the prepared deposit/ package
    (files + minimal dataset metadata) to THIS draft, then PUBLISH.
  - Prints the resulting record state and the live DOI.

Usage:
  set ZENODO_TOKEN=your_token
  python3 publish_draft.py [--draft 21823584]
"""
import argparse, json, os, subprocess, sys

PROXY = "socks5h://127.0.0.1:8428"
API = "https://zenodo.org/api"
DEP_DIR = os.path.dirname(os.path.abspath(__file__))

# Files to include (everything the reviewer needs to reproduce). Skip the script itself.
INCLUDE = [
    "signature_101genes.csv",
    "per_sample_scores_GSE227919.csv",
    "per_sample_scores_GSE85195.csv",
    "GSE26549_progression_DEG.csv",
    "session_info.txt",
    "README.md",
    "scripts",
]

def curl(method, url, token, data=None):
    cmd = ["curl", "-sS", "-x", PROXY, "-X", method, url,
           "-H", f"Authorization: Bearer {token}"]
    if data is not None:
        cmd += ["-H", "Content-Type: application/json", "-d", json.dumps(data)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl failed ({method} {url}): {r.stderr}")
    return r.stdout

def upload_file(bucket, name, path, token):
    """PUT a local file's content to the Zenodo bucket."""
    cmd = ["curl", "-sS", "-x", PROXY, "-X", "PUT",
           "-H", f"Authorization: Bearer {token}",
           "-T", path, f"{bucket}/{name}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"upload failed {name}: {r.stderr}")
    return r.stdout

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draft", type=str, default="21823584")
    ap.add_argument("--dry", action="store_true", help="only inspect, do not publish")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ERROR: ZENODO_TOKEN not set. Export it first (set ZENODO_TOKEN=...).")
    if token.lower().startswith("your") or len(token) < 20:
        sys.exit("ERROR: ZENODO_TOKEN looks like a placeholder. Provide the real token.")

    # 1) fetch the draft
    url = f"{API}/deposit/depositions/{args.draft}"
    out = curl("GET", url, token)
    try:
        dep = json.loads(out)
    except json.JSONDecodeError:
        sys.exit(f"ERROR: cannot parse Zenodo response (token valid? record exists?)\n{out[:500]}")
    print(f"draft {args.draft}: state={dep.get('state')}, title={dep.get('metadata',{}).get('title','<none>')}")
    print(f"  files present: {len(dep.get('files', []))}")
    print(f"  links: {list(dep.get('links', {}).keys())}")

    # 2) if empty, fill it
    meta = dep.get("metadata") or {}
    has_title = bool(meta.get("title"))
    has_files = len(dep.get("files", [])) > 0

    if not has_files or not has_title:
        print("-- draft incomplete; populating files + metadata --")
        bucket = dep.get("links", {}).get("bucket")
        if not bucket:
            sys.exit("ERROR: draft has no bucket link; cannot upload files. "
                     "Create the deposition through the web UI or let me create a new one.")
        for name in INCLUDE:
            p = os.path.join(DEP_DIR, name)
            if os.path.isdir(p):
                for root, _, fnames in os.walk(p):
                    for fn in fnames:
                        fp = os.path.join(root, fn)
                        print(f"  upload {fp}")
                        upload_file(bucket, os.path.basename(fp), fp, token)
            elif os.path.exists(p):
                print(f"  upload {p}")
                upload_file(bucket, name, p, token)
        # set minimal metadata
        md = {
            "metadata": {
                "upload_type": "dataset",
                "title": "Direction-weighted 101-gene progression signature for oral leukoplakia / proliferative verrucous leukoplakia (PVL): reanalysis and cross-cohort validation",
                "creators": [{"name": "tiduslrg"}],
                "description": "Reanalysis of public GEO transcriptomic cohorts (GSE26549, GSE227919, GSE196296, GSE85195) defining a direction-weighted 101-gene progression signature for oral leukoplakia / PVL, with cross-cohort validation scores and single-cell cluster statistics. Analysis scripts (R 4.3.3) and per-sample scores included.",
                "access_right": "open",
                "license": "cc-by-4.0",
            }
        }
        curl("PUT", url, token, data=md)
        print("-- re-fetching draft after populate --")
        dep = json.loads(curl("GET", url, token))
        print(f"  now files={len(dep.get('files', []))}, title={dep.get('metadata',{}).get('title','<none>')[:40]}...")

    # 3) publish
    publish_url = dep.get("links", {}).get("publish")
    if not publish_url:
        # fallback
        publish_url = f"{url}/actions/publish"
    if args.dry:
        print(f"[dry-run] would PUBLISH via {publish_url}")
        return
    resp = curl("POST", publish_url, token)
    try:
        rec = json.loads(resp)
    except json.JSONDecodeError:
        sys.exit(f"ERROR publishing: {resp[:500]}")
    print(f"PUBLISHED. record id={rec.get('id')}, DOI=https://doi.org/{rec.get('metadata',{}).get('doi')}")
    print(f"  state={rec.get('state')}")

if __name__ == "__main__":
    main()
