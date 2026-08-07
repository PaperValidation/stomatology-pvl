import urllib.request, os, sys, subprocess

PY = "C:/Users/tidus/.workbuddy/binaries/python/versions/3.13.12/python.exe"
OUT = "C:/Users/tidus/PycharmProjects/Stomatology-PVL/results/saintigny_overlap"
os.makedirs(OUT, exist_ok=True)

# ---- 1) install xlrd (aliyun mirror, no proxy) ----
subprocess.run([PY, "-m", "pip", "install", "-q",
                "-i", "https://mirrors.aliyun.com/pypi/simple/",
                "--trusted-host", "mirrors.aliyun.com", "xlrd"], check=False)

# ---- 2) download Saintigny 2011 Supplementary Material 9 (2,182 Cox-significant transcripts) ----
URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/3074595/bin/NIHMS261208-supplement-9.xls"
XLS = os.path.join(OUT, "saintigny_2182.xls")

def try_get(proxy=None):
    import socks, socket
    if proxy:
        host, port = proxy
        socks.set_default_proxy(socks.SOCKS5, host, port)
        socket.socket = socks.socksocket
    else:
        # reset to plain
        socket.socket = socket._orig_socket if hasattr(socket, "_orig_socket") else socket.socket
    return urllib.request.urlopen(URL, timeout=40).read()

data, how = None, "none"
for attempt in [("direct", None), ("socks8428", ("127.0.0.1", 8428)), ("socks1987", ("127.0.0.1", 1987))]:
    label, proxy = attempt
    try:
        data = try_get(proxy)
        how = label
        break
    except Exception as e:
        print(f"[{label}] failed: {e}", file=sys.stderr)

if not data:
    print("DOWNLOAD_FAILED")
    sys.exit(2)
with open(XLS, "wb") as f:
    f.write(data)
print(f"downloaded via {how}: {len(data)} bytes -> {XLS}")

# ---- 3) parse ----
import pandas as pd
xl = pd.ExcelFile(XLS, engine="xlrd")
print("sheets:", xl.sheet_names)
df = xl.parse(xl.sheet_names[0])
print("shape:", df.shape)
print("cols:", df.columns.tolist())
print(df.head(3).to_string())

# try to find a gene-symbol-ish column
col = None
for c in df.columns:
    cl = str(c).lower()
    if "gene" in cl and ("symbol" in cl or "sym" in cl):
        col = c; break
if col is None:
    for c in df.columns:
        cl = str(c).lower()
        if "symbol" in cl or cl == "gene":
            col = c; break
if col is None:
    # fallback: first column
    col = df.columns[0]
print("gene-symbol column:", repr(col))

syms_2182 = set()
for v in df[col].dropna():
    s = str(v).strip()
    if s and s != "NA" and s.upper() != "NAN":
        # take first if multiple separated
        s = s.split("//")[0].split(";")[0].split(",")[0].strip()
        if s:
            syms_2182.add(s)
print("2182-list unique symbols parsed:", len(syms_2182))

# ---- 4) our 101 genes ----
our = set()
with open(os.path.join(OUT, "our_101_genes.txt")) as f:
    for line in f:
        g = line.strip()
        if g:
            our.add(g)
print("our 101 genes:", len(our))

# ---- 5) Saintigny 29-transcript model (Table 2 Model 1) ----
model29_raw = [
    ("8095441","CSN1S2A"),("8023314","CCDC110"),("7986442","ENST00000391004"),
    ("8062842","ENST00000387867"),("8084002","KCNMB2"),("7915846","MKNK1"),
    ("8165709","NC_001807"),("8122200","ENST00000385892"),("8046408","PDK1"),
    ("8153223","PTK2"),("8172119","MED14"),("8061092","NA"),
    ("7927106","ENST00000387096"),("7948894","RNU2-1"),("8083939","AK128090"),
    ("3939865","OR4B1"),("916777","hsa-mir-101-1"),("7964360","STAT6"),
    ("8101762","SNCA"),("7962489","PLEKHA9"),("7901361","ENST00000387793"),
    ("8044682","SNRPA1"),("8028950","CYP2G1P"),("7977480","ENST00000386651"),
    ("8067983","ENST00000387011"),("8097743","ENST00000410285"),("8093957","CNO"),
    ("8086536","ENST00000365398"),("8121943","ENST00000384255"),
]
# mappable protein/pseudo gene symbols (drop NA, ENST*, NC_*, hsa-mir-*, AK*, RNU2-1)
def is_mappable(s):
    if s in ("NA",): return False
    if s.startswith("ENST"): return False
    if s.startswith("NC_"): return False
    if s.startswith("hsa-mir"): return False
    if s.startswith("AK"): return False
    if s.startswith("RNU"): return False
    return True
model29_syms = set(s for _, s in model29_raw if is_mappable(s))
print("model29 mappable symbols:", len(model29_syms), sorted(model29_syms))

# ---- 6) overlaps / Jaccard ----
def jac(a, b):
    inter = a & b
    union = a | b
    return inter, union, (len(inter)/len(union) if union else 0.0)

inter29, union29, j29 = jac(our, model29_syms)
inter2182, union2182, j2182 = jac(our, syms_2182)

print("\n=== RESULTS ===")
print(f"Our 101 genes: {len(our)}")
print(f"Saintigny 29-model (mappable syms): {len(model29_syms)}")
print(f"  intersection: {len(inter29)} -> {sorted(inter29)}")
print(f"  Jaccard(101, 29model) = {j29:.4f}")
print(f"Saintigny 2,182 Cox-significant: {len(syms_2182)}")
print(f"  intersection: {len(inter2182)}")
print(f"  Jaccard(101, 2182) = {j2182:.4f}")
print(f"  our 101 within 2182? {len(inter2182)}/{len(our)} = {len(inter2182)/len(our):.3f}")

# ---- 7) write result csv ----
import csv
res = os.path.join(OUT, "saintigny_overlap_results.csv")
with open(res, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["our_101_gene_count", len(our)])
    w.writerow(["saintigny_29model_mappable_sym_count", len(model29_syms)])
    w.writerow(["overlap_101_vs_29model", len(inter29)])
    w.writerow(["jaccard_101_vs_29model", round(j29,4)])
    w.writerow(["saintigny_2182_count", len(syms_2182)])
    w.writerow(["overlap_101_vs_2182", len(inter2182)])
    w.writerow(["jaccard_101_vs_2182", round(j2182,4)])
    w.writerow(["frac_our101_in_2182", round(len(inter2182)/len(our),3)])
print("wrote", res)
