import gzip, csv
from collections import Counter

# --- resolve DATASET_PATH from local .env (no external deps) ---
import os as _os, re as _re
def _load_dataset_path():
    _here = _os.path.dirname(_os.path.abspath(__file__))
    _ds = None
    _envp = _os.path.join(_here, ".env")
    if _os.path.exists(_envp):
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
RAW = _os.path.join(DATASET_PATH, "raw")
MET = _os.path.join(DATASET_PATH, "metadata")

p = _os.path.join(RAW, "GSE26549_series_matrix.txt.gz")
# store each metadata line as list of cell values; keep multi-line characteristics in order
meta_single = {}          # key -> value string (for non-repeated keys)
char_blocks = []          # list of value-lists for each !Sample_characteristics_ch1 line
sample_order = None
in_table = False
with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_table = True; continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            sample_order = [x.strip().strip('"') for x in line.rstrip('\n').split('\t')]
            break
        if line.startswith('!'):
            parts = line.rstrip('\n').split('\t', 1)
            if len(parts) != 2: continue
            key, val = parts
            if key == '!Sample_characteristics_ch1':
                char_blocks.append(val.split('\t'))
            else:
                meta_single[key] = val

assert sample_order is not None
sample_order = sample_order[1:]   # drop ID_REF
n = len(sample_order)
print("N samples:", n, "| first3:", sample_order[:3], "| n_char_blocks:", len(char_blocks))

def stripq(x): return x.strip().strip('"')
samples = []
for i in range(n):
    rec = {'GSM': sample_order[i]}
    for block in char_blocks:
        if i < len(block):
            cell = stripq(block[i])
            if ':' in cell:
                label, v = cell.split(':', 1)
                rec[label.strip()] = v.strip()
            else:
                rec.setdefault('_other_%d' % i, cell)
    samples.append(rec)

keys = set()
for s in samples: keys.update(s.keys())
print("Characteristic labels:", sorted(keys))

def norm_outcome(o):
    if o is None: return 'NA'
    o = o.lower()
    if 'no oral cancer development' in o: return 'non_progressing'
    if 'oral cancer development' in o: return 'progressing'
    return o

out_path = _os.path.join(MET, "GSE26549_clinical.csv")
fields = ['GSM','Sex','race','age','smoking habits','alcohol habits','treatment arm',
          'histology at baseline (hyperplasia versus dysplasia)',
          'histology at baseline (breakdown)','p63 expression','podoplanin expression',
          'oral cancer-free survival time (years)','outcome','time of biopsy']
with open(out_path, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    for s in samples:
        row = dict(s); row['outcome'] = norm_outcome(s.get('outcome')); w.writerow(row)

print("\n=== Outcome (all) ===", dict(Counter(s['outcome'] for s in samples)))
print("=== Time of biopsy ===", dict(Counter(s.get('time of biopsy','?') for s in samples)))
base = [s for s in samples if 'baseline' in (s.get('time of biopsy') or '')]
print("=== Baseline-only (n=%d) outcome ===" % len(base), dict(Counter(s['outcome'] for s in base)))
print("CSV ->", out_path)
