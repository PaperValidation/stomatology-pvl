#!/usr/bin/env python3
# Download GSE196296 10x single-cell counts (LesionA/B/C) to raw dir.
import socks, socket, urllib.request, os, gzip, io

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

PROXY = ("127.0.0.1", 8428)
OUT = RAW
os.makedirs(OUT, exist_ok=True)

samples = {
    "GSM5870060": "LesionA",
    "GSM5870062": "LesionB",
    "GSM5870063": "LesionC",
}
base = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5870nnn/"

def fetch(url, binary=True):
    socks.set_default_proxy(socks.SOCKS5, *PROXY)
    socket.socket = socks.socksocket
    req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()

ok = 0; fail = 0
for gsm, lesion in samples.items():
    for suf in ["barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"]:
        fname = f"{gsm}_{lesion}_{suf}"
        url = base + gsm + "/suppl/" + fname
        outp = os.path.join(OUT, fname)
        if os.path.exists(outp) and os.path.getsize(outp) > 0:
            print("skip", fname); ok += 1; continue
        try:
            data = fetch(url)
            with open(outp, "wb") as f:
                f.write(data)
            print(f"OK  {fname}  {len(data)//1024} KB")
            ok += 1
        except Exception as e:
            print(f"ERR {fname}: {type(e).__name__} {e}")
            fail += 1
print(f"\nDONE ok={ok} fail={fail}")
