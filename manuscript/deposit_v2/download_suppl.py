import gzip, os, socket, socks, urllib.request

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

OUT = RAW
os.makedirs(OUT, exist_ok=True)
socks.set_default_proxy(socks.SOCKS5, "127.0.0.1", 8428)
socket.socket = socks.socksocket

base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227919/suppl/"
files = [
    "GSE227919_2023_20_03_expr_mat_raw.txt.gz",
    "GSE227919_2023_20_03_expr_mat_microbe.txt.gz",
    "GSE227919_2023_21_03_pml_annotations.csv.gz",
]
for fn in files:
    url = base + fn
    dst = os.path.join(OUT, fn)
    if os.path.exists(dst):
        print(f"[skip] {fn}")
        continue
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as r, open(dst, "wb") as w:
            w.write(r.read())
        print(f"[ok]   {fn} -> {os.path.getsize(dst)} bytes")
    except Exception as e:
        print(f"[FAIL] {fn}: {type(e).__name__}: {e}")
print("SUPPL DONE")
