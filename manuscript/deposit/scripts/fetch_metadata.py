import urllib.request, gzip, os, csv

# --- resolve DATASET_PATH from local .env (no external deps) ---
import re as _re
def _load_dataset_path():
    _here = os.path.dirname(os.path.abspath(__file__))
    _ds = None
    _envp = os.path.join(_here, ".env")
    if os.path.exists(_envp):
        with open(_envp, encoding="utf-8") as _f:
            for _line in _f:
                _s = _line.strip()
                if not _s or _s.startswith("#"):
                    continue
                _m = _re.match(r"^DATASET_PATH\s*=\s*(.+)$", _s)
                if _m:
                    _ds = _m.group(1).strip().strip("'\"").strip()
                    break
    if not _ds:
        _ds = r"G:\CMCC Sync\Medical Data\Stomatology\PVL"
    return _ds
DATASET_PATH = _load_dataset_path()
RAW = os.path.join(DATASET_PATH, "raw")
MET = os.path.join(DATASET_PATH, "metadata")

OUTDIR = MET
os.makedirs(OUTDIR, exist_ok=True)

PROXIES = [None, 8428, 8429]
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

# correct GEO FTP prefix = first (len-3) digits of the series number
GSES = {
    "GSE26549":  "GSE26nnn/GSE26549/matrix/GSE26549_series_matrix.txt.gz",
    "GSE227919": "GSE227nnn/GSE227919/matrix/GSE227919_series_matrix.txt.gz",
    "GSE196296": "GSE196nnn/GSE196296/matrix/GSE196296_series_matrix.txt.gz",
}

def open_url(url):
    last_err = None
    for p in PROXIES:
        try:
            if p:
                import socks, socket
                socks.set_default_proxy(socks.SOCKS5, "127.0.0.1", p)
                socket.socket = socks.socksocket
            else:
                import socket
                socket.socket = socket.socket
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            return urllib.request.urlopen(req, timeout=90)
        except Exception as e:
            last_err = e
    raise last_err

def fetch_metadata(gse, rel):
    url = BASE + rel
    print(f"\n>>> {gse}: {url}")
    resp = open_url(url)
    sample_rows = []   # (field, [vals]) for Sample_* only
    series_info = {}   # field -> single value
    decomp = gzip.GzipFile(fileobj=resp)
    for raw in decomp:
        line = raw.decode("utf-8", "replace").rstrip("\n")
        if not line.startswith("!"):
            continue
        if line.strip() == "!series_matrix_table_begin":
            break
        parts = line[1:].split("\t")
        field = parts[0]
        vals = parts[1:]
        if field.startswith("Sample_"):
            sample_rows.append((field, vals))
        elif field.startswith("Series_"):
            if vals:
                series_info[field] = vals[0]
    resp.close()
    n = max((len(v) for _, v in sample_rows), default=0)
    samples = []
    for i in range(n):
        rec = {}
        for field, vals in sample_rows:
            rec[field] = vals[i] if i < len(vals) else ""
        samples.append(rec)
    print(f"    samples={n}")
    # dump full CSV (sample-level only)
    if samples:
        cols = list(samples[0].keys())
        path = os.path.join(OUTDIR, f"{gse}_sampletable.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for s in samples:
                w.writerow(s)
        print(f"    -> {path}")
    return series_info, samples

def concise(gse, samples):
    if not samples:
        return
    acc = next((k for k in samples[0] if k == "Sample_geo_accession"), "Sample_geo_accession")
    src = "Sample_source_name_ch1"
    cha = "Sample_characteristics_ch1"
    title = "Sample_title"
    print(f"    === {gse} sample list (n={len(samples)}) ===")
    for s in samples:
        c = s.get(cha, "").replace("\n", " ")
        t = s.get(title, "")
        print(f"      {s.get(acc,''):<12} | {s.get(src,'')[:22]:<22} | {t[:30]:<30} | {c[:60]}")

def group_scan(gse, samples):
    # count keyword occurrences to locate outcome/group columns
    print(f"    === {gse} outcome/group token scan ===")
    toks = {}
    for s in samples:
        for k, v in s.items():
            if "characteristics" in k or "source_name" in k or "title" in k or "group" in k or "description" in k:
                low = (v or "").lower()
                for kw in ["progress","non-progress","transform","malignant","cancer","carcinoma",
                           "normal","low risk","high risk","low-risk","high-risk","risk","dysplasia",
                           "opl","opmd","pvl","leukoplak","male","female","smok","age","group","case","control","involved","uninvolved"]:
                    if kw in low:
                        toks.setdefault(kw,0); toks[kw]+=1
    for kw,c in sorted(toks.items(), key=lambda x:-x[1]):
        print(f"      '{kw}': {c}")
    # show unique values of likely group columns
    for k in samples[0]:
        if k in ("Sample_source_name_ch1",) or "group" in k.lower():
            uniq = {}
            for s in samples:
                uniq[s.get(k,"")] = uniq.get(s.get(k,""),0)+1
            if len(uniq) <= 12:
                print(f"    unique [{k}]:")
                for val,cnt in sorted(uniq.items(), key=lambda x:-x[1]):
                    print(f"        {cnt:>3}  {val[:70]}")

for gse, rel in GSES.items():
    try:
        sinfo, samples = fetch_metadata(gse, rel)
        concise(gse, samples)
        group_scan(gse, samples)
    except Exception as e:
        print(f"    ERROR: {type(e).__name__}: {e}")
print("\nDONE.")
