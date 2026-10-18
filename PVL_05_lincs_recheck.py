#!/usr/bin/env python
"""Reproduce the single-direction LINCS query of PVL_05_lincs.R and audit the
claim of "727 nominally significant compounds".

Mirrors the R logic:
  prot_down <- enrichr(up_genes,   "LINCS_L1000_Chem_Pert_down")
  prot_up   <- enrichr(down_genes, "LINCS_L1000_Chem_Pert_up")
  promo_up  <- enrichr(up_genes,   "LINCS_L1000_Chem_Pert_up")
  promo_down<- enrichr(down_genes, "LINCS_L1000_Chem_Pert_down")
Significance threshold used downstream in R: Adjusted.P.value < 0.05.

Writes full (untruncated) result tables, one per query, plus a summary of
unique-compound counts under several definitions of "significant".
No existing PVL5_*.csv file is modified.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

import pandas as pd

RES = r"F:\Data Temp\PVL"
BASE = "https://maayanlab.cloud/Enrichr"
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def post_multipart(path, fields, files):
    """Enrichr /addList only accepts multipart/form-data with a file part."""
    boundary = "----PVLboundary7MA4YWxkTrZu0gW"
    body = b""
    for k, v in fields.items():
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n"
                 f"{v}\r\n").encode()
    for k, (fname, content) in files.items():
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"; "
                 f"filename=\"{fname}\"\r\nContent-Type: text/plain\r\n\r\n").encode()
        body += content.encode() + b"\r\n"
    body += f"--{boundary}--\r\n".encode()
    req = urllib.request.Request(BASE + path, data=body, headers={
        "Content-Type": f"multipart/form-data; boundary={boundary}"})
    for attempt in range(4):
        try:
            with OPENER.open(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"  retry {attempt+1} after {e}", flush=True)
            time.sleep(3 * (attempt + 1))


def add_list(genes, desc):
    res = post_multipart("/addList", {"description": desc},
                         {"list": ("genes.txt", "\n".join(genes) + "\n")})
    if "userListId" not in res:
        raise RuntimeError(f"addList failed: {res}")
    return res["userListId"]


def enrich(uid, lib):
    q = urllib.parse.urlencode({"userListId": uid, "backgroundType": lib})
    req = urllib.request.Request(f"{BASE}/enrich?{q}")
    for attempt in range(4):
        try:
            with OPENER.open(req, timeout=120) as r:
                return json.load(r).get(lib, [])
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            print(f"  retry {attempt+1} after {e}", flush=True)
            time.sleep(3 * (attempt + 1))


def to_frame(rows, lib, query, query_size):
    cols = ["Term", "Overlap", "Overlap.N", "Query.Size", "P.value", "Adjusted.P.value",
            "Odds.Ratio", "Combined.Score", "Genes", "direction", "library"]
    if not rows:
        return pd.DataFrame(columns=cols)
    recs = []
    for r in rows:
        # Enrichr /enrich row: [rank, term, p, zscore, combined, genes(list), adjp, oldp, oldadj]
        genes = r[5] if isinstance(r[5], list) else str(r[5]).split(";")
        genes = [g for g in genes if g]
        recs.append({
            "Term": r[1],
            "Overlap": f"{len(genes)}/{query_size}",
            "Overlap.N": len(genes),
            "Query.Size": query_size,
            "P.value": r[2],
            "Adjusted.P.value": r[6],
            "Odds.Ratio": "",  # not returned by Enrichr API v3
            "Combined.Score": r[4],
            "Genes": ";".join(genes),
            "direction": query,
            "library": lib,
        })
    return pd.DataFrame(recs, columns=cols)


DRUG_RE = re.compile(r"-([^-]+)-[^-]*$")


def drug_of(term):
    """Same intent as the R helper: take the compound token from the LINCS term."""
    parts = term.split("-")
    if len(parts) >= 3:
        return "-".join(parts[1:-1])
    return term


def main():
    sig = pd.read_csv(os.path.join(RES, "PVL5_signature_signed.csv"))
    up = [g for g in sig.loc[sig.direction == "up_in_progression", "genes"].iloc[0].split(";") if g]
    down = [g for g in sig.loc[sig.direction == "down_in_progression", "genes"].iloc[0].split(";") if g]
    print(f"signature: up={len(up)} down={len(down)}")

    uid_up = add_list(up, "PVL_progression_up_26")
    uid_down = add_list(down, "PVL_progression_down_60")
    print(f"userListId up={uid_up} down={uid_down}")

    queries = OrderedDict([
        ("up_geneset__Chem_Pert_down", (uid_up, "LINCS_L1000_Chem_Pert_down", "up_in_progression", len(up))),
        ("down_geneset__Chem_Pert_up", (uid_down, "LINCS_L1000_Chem_Pert_up", "down_in_progression", len(down))),
        ("up_geneset__Chem_Pert_up", (uid_up, "LINCS_L1000_Chem_Pert_up", "up_in_progression", len(up))),
        ("down_geneset__Chem_Pert_down", (uid_down, "LINCS_L1000_Chem_Pert_down", "down_in_progression", len(down))),
    ])

    frames = {}
    for key, (uid, lib, direction, qsize) in queries.items():
        rows = enrich(uid, lib)
        df = to_frame(rows, lib, direction, qsize)
        df["Drug"] = df["Term"].map(drug_of)
        frames[key] = df
        print(f"{key}: {lib} -> {len(df)} terms, "
              f"{df['Drug'].nunique()} unique compounds, "
              f"{df.loc[df['Adjusted.P.value'] < 0.05, 'Drug'].nunique()} with adj.P<0.05")

    f1 = os.path.join(RES, "LINCS_single_direction_up_geneset.csv")
    f2 = os.path.join(RES, "LINCS_single_direction_down_geneset.csv")
    frames["up_geneset__Chem_Pert_down"].to_csv(f1, index=False)
    frames["down_geneset__Chem_Pert_up"].to_csv(f2, index=False)
    print("wrote", f1, len(frames["up_geneset__Chem_Pert_down"]))
    print("wrote", f2, len(frames["down_geneset__Chem_Pert_up"]))

    for key, df in frames.items():
        sfx = key.replace("__", "_")
        df.to_csv(os.path.join(RES, f"LINCS_raw_{sfx}.csv"), index=False)

    # ---- count "nominally significant compounds" under several definitions ----
    a = frames["up_geneset__Chem_Pert_down"]   # protective arm 1
    b = frames["down_geneset__Chem_Pert_up"]   # protective arm 2
    c = frames["up_geneset__Chem_Pert_up"]     # promoting arm 1
    d = frames["down_geneset__Chem_Pert_down"]  # promoting arm 2

    def nuniq(df, col="Adjusted.P.value", thr=0.05):
        return int(df.loc[df[col] < thr, "Drug"].nunique())

    def nterm(df, col="Adjusted.P.value", thr=0.05):
        return int((df[col] < thr).sum())

    stats = {
        "prot_down_adjp_lt_0.05_drugs": nuniq(a),
        "prot_down_adjp_lt_0.05_terms": nterm(a),
        "prot_up_adjp_lt_0.05_drugs": nuniq(b),
        "prot_up_adjp_lt_0.05_terms": nterm(b),
        "prot_union_adjp_lt_0.05_drugs": len(set(a.loc[a["Adjusted.P.value"] < 0.05, "Drug"])
                                              | set(b.loc[b["Adjusted.P.value"] < 0.05, "Drug"])),
        "prot_intersection_adjp_lt_0.05_drugs": len(set(a.loc[a["Adjusted.P.value"] < 0.05, "Drug"])
                                                    & set(b.loc[b["Adjusted.P.value"] < 0.05, "Drug"])),
        "prot_down_nominal_p_lt_0.05_drugs": nuniq(a, "P.value"),
        "prot_up_nominal_p_lt_0.05_drugs": nuniq(b, "P.value"),
        "prot_union_nominal_p_lt_0.05_drugs": len(set(a.loc[a["P.value"] < 0.05, "Drug"])
                                                   | set(b.loc[b["P.value"] < 0.05, "Drug"])),
        "all_four_union_adjp_lt_0.05_drugs": len(set().union(*[
            set(x.loc[x["Adjusted.P.value"] < 0.05, "Drug"]) for x in (a, b, c, d)])),
        "all_four_union_nominal_p_lt_0.05_drugs": len(set().union(*[
            set(x.loc[x["P.value"] < 0.05, "Drug"]) for x in (a, b, c, d)])),
        "prot_down_any_term_drugs": int(a["Drug"].nunique()),
        "prot_up_any_term_drugs": int(b["Drug"].nunique()),
    }
    print("\n=== counts ===")
    for k, v in stats.items():
        print(f"{k} = {v}")
    with open(os.path.join(RES, "LINCS_727_recheck_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"signature": {"up": len(up), "down": len(down)},
                   "libs": ["LINCS_L1000_Chem_Pert_up", "LINCS_L1000_Chem_Pert_down"],
                   "userListIds": {"up": uid_up, "down": uid_down},
                   "counts": stats}, fh, indent=2)

    # reversal intersection vs the archived candidates file
    cand_path = os.path.join(RES, "PVL5_lincs_reversal_candidates.csv")
    cand = pd.read_csv(cand_path)
    print(f"\narchived PVL5_lincs_reversal_candidates.csv: {len(cand)} rows, "
          f"consensus_confirmed TRUE = {int(cand['consensus_confirmed'].sum())}")
    return stats


if __name__ == "__main__":
    main()
