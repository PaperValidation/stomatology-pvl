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

p = _os.path.join(RAW, "GSE227919_2023_21_03_pml_annotations.csv.gz")
out = _os.path.join(MET, "GSE227919_clinical.csv")

with gzip.open(p, 'rt', encoding='utf-8-sig', errors='replace') as f:
    r = csv.reader(f)
    header = next(r)
    idx = {h: i for i, h in enumerate(header)}
    rows = list(r)

def g(row, name):
    i = idx.get(name)
    return row[i].strip() if (i is not None and i < len(row)) else ''

out_rows = []
for row in rows:
    sid = g(row, 'Sample_ID')
    if not sid:
        continue
    out_rows.append({
        'Sample_ID': sid,
        'Class': g(row, 'Class'),
        'Progression_status': g(row, 'Progression_status'),
        'Sex': g(row, 'Sex'),
        'Age': g(row, 'Age'),
        'Smoking_status': g(row, 'Smoking_status'),
        'Histopathology': g(row, 'Histopathology').replace('"', ''),
    })

with open(out, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=['Sample_ID','Class','Progression_status','Sex','Age','Smoking_status','Histopathology'])
    w.writeheader()
    w.writerows(out_rows)

print("N samples:", len(out_rows))
print("Class counts:", dict(Counter(r['Class'] for r in out_rows)))
print("Progression_status counts:", dict(Counter(r['Progression_status'] for r in out_rows)))
print("CSV ->", out)
print("Sample rows:")
for r in out_rows[:6]:
    print(" ", r['Sample_ID'], r['Class'], '|', r['Progression_status'], '|', r['Sex'], r['Age'])
