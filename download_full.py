import urllib.request, os, gzip, shutil

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

OUTDIR = RAW
os.makedirs(OUTDIR, exist_ok=True)
PROXIES = [None, 8428, 8429]
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

TARGETS = {
    "GSE26549": "GSE26nnn/GSE26549/matrix/GSE26549_series_matrix.txt.gz",
    "GSE227919": "GSE227nnn/GSE227919/matrix/GSE227919_series_matrix.txt.gz",
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
            return urllib.request.urlopen(req, timeout=300)
        except Exception as e:
            last_err = e
    raise last_err

for gse, rel in TARGETS.items():
    url = BASE + rel
    dest = os.path.join(OUTDIR, f"{gse}_series_matrix.txt.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"[skip] {gse} already exists ({os.path.getsize(dest)} bytes)")
        continue
    print(f"[down] {gse} <- {url}")
    try:
        resp = open_url(url)
        with open(dest, "wb") as f:
            shutil.copyfileobj(resp, f, length=1024*1024)
        resp.close()
        print(f"[ok]   {gse} -> {dest} ({os.path.getsize(dest)} bytes)")
    except Exception as e:
        print(f"[ERR]  {gse}: {type(e).__name__}: {e}")
print("FULL-DOWNLOAD DONE")
