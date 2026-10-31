"""Build-time only: a gzipped FASTA -> genome.fa.gz holding only the chromosomes in chrom_sizes.tsv,
plus genome.fa.fai (the samtools/pyfaidx index of its uncompressed form). At run time the motif step
decompresses it to /tmp and puts the index next to it: pyfaidx reads BGZF records from their start,
which makes random access into a compressed chromosome far too slow.
Usage: python primary_fasta.py <in.fa.gz> <chrom_sizes.tsv> <out dir>"""
import gzip
import shutil
import sys
from pathlib import Path

from pyfaidx import Fasta

src, sizes_file, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
sizes = dict(line.split("\t") for line in open(sizes_file).read().splitlines())
plain = out / "genome.fa"
written, on = set(), False
with gzip.open(src, "rb") as fin, open(plain, "wb") as fout:
    for line in fin:
        if line.startswith(b">"):
            name = line[1:].split()[0].decode()
            on = name in sizes
            if on:
                written.add(name)
                line = f">{name}\n".encode()
        if on:
            fout.write(line)
if set(sizes) - written:
    sys.exit(f"chromosomes missing from {src}: {sorted(set(sizes) - written)}")
fa = Fasta(str(plain))                    # writes genome.fa.fai
bad = {c: len(fa[c]) for c in sizes if len(fa[c]) != int(sizes[c])}
fa.close()
if bad:
    sys.exit(f"chromosome lengths differ from chrom_sizes.tsv: {bad}")
with open(plain, "rb") as fin, gzip.open(out / "genome.fa.gz", "wb", compresslevel=6) as fout:
    shutil.copyfileobj(fin, fout, 1 << 24)
plain.unlink()
