# Case marking for the cross-assembly alignments: the source row is 300 bp flank + copy + 300 bp flank (unaligned input
# aln/<id>.in.fa on KIT, first record). Columns holding the source copy are written UPPERCASE in every row; all other
# columns (flanks) lowercase.
import sys, glob, csv, os
from Bio import SeqIO
SHOW = 300
L = {r[0]: int(r[4]) - int(r[3]) + 1 for r in csv.reader(open(sys.argv[1], encoding='utf-8'), delimiter='\t')}
for fn in sorted(glob.glob('traced/*.aln.fa')):
    i = os.path.basename(fn).split('.')[0]
    recs = list(SeqIO.parse(fn, 'fasta'))
    src = str(recs[0].seq)
    k = 0; copycols = set()
    for c, ch in enumerate(src):
        if ch == '-': continue
        if SHOW <= k < SHOW + L[i]: copycols.add(c)
        k += 1
    with open(fn, 'w', encoding='utf-8', newline='\n') as f:
        for r in recs:
            s = ''.join(ch.upper() if c in copycols else ch.lower() for c, ch in enumerate(str(r.seq)))
            f.write(f'>{r.description}\n{s}\n')
    print(i, 'copy columns', len(copycols), 'of', len(src))
