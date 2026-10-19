from __future__ import annotations
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np

class ModelError(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise ModelError(message)

def read_table(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_table(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def save_json(path, value):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x))
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=convert, allow_nan=False), encoding='utf-8')

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def text_open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def normalize(x):
    x = np.asarray(x, float)
    require(np.all(np.isfinite(x)) and np.all(x >= 0) and x.sum() > 0, 'Invalid probability weights')
    return x / x.sum()

class Fasta:
    """Read indexed FASTA without loading a human genome into memory."""
    def __init__(self,path):
        self.path=Path(path)
        require(Path(str(path)+'.fai').is_file(),'FASTA .fai index required')
        self.index={}
        for line in Path(str(path)+'.fai').read_text().splitlines():
            a=line.split('\t'); self.index[a[0]]=tuple(map(int,a[1:5]))
        self.handle=self.path.open('rb')

    def fetch(self,chrom,start,end):
        require(chrom in self.index,f'Contig absent in FASTA: {chrom}')
        length,offset,bases,width=self.index[chrom]
        require(0<=start<=end<=length,f'Invalid reference interval: {chrom}:{start}-{end}')
        pieces=[]
        while start<end:
            take=min(end-start,bases-start%bases)
            self.handle.seek(offset+(start//bases)*width+start%bases)
            pieces.append(self.handle.read(take).decode('ascii').upper()); start+=take
        return ''.join(pieces)

    def close(self): self.handle.close()
