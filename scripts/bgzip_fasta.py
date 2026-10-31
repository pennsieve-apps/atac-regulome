"""Build-time only: a gzipped FASTA -> BGZF-compressed FASTA of the chromosomes in chrom_sizes.tsv,
indexed (.fai, .gzi) with pyfaidx, so the motif step reads sequences by random access without
holding the genome. Usage: python bgzip_fasta.py <in.fa.gz> <chrom_sizes.tsv> <out.fa.bgz>"""
import gzip
import sys

from Bio import bgzf
from pyfaidx import Fasta

src, sizes, out = sys.argv[1:4]
keep = {line.split("\t")[0] for line in open(sizes).read().splitlines()}
written, on = set(), False
with gzip.open(src, "rb") as fin, bgzf.BgzfWriter(out, "wb") as fout:
    for line in fin:
        if line.startswith(b">"):
            name = line[1:].split()[0].decode()
            on = name in keep
            if on:
                written.add(name)
                line = f">{name}\n".encode()
        if on:
            fout.write(line)
missing = keep - written
if missing:
    sys.exit(f"chromosomes missing from {src}: {sorted(missing)}")
fa = Fasta(out)
bad = {c: len(fa[c]) for c in keep if len(fa[c]) != int(dict(l.split("\t") for l in open(sizes).read().splitlines())[c])}
fa.close()
if bad:
    sys.exit(f"chromosome lengths differ from chrom_sizes.tsv: {bad}")
