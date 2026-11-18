# Stacked monomer view of a solo copy: cut the copy-vs-tiled-consensus alignment (aln/<id>.fa) at consensus position 1;
# one row per monomer, columns = consensus positions 1..L plus the insertion columns any monomer needs. No realignment.
import sys, csv
from Bio import SeqIO
CONS = {r.id: str(r.seq).upper().replace('-', '') for r in SeqIO.parse(sys.argv[1], 'fasta')}
def stack(fn):
    recs = list(SeqIO.parse(fn, 'fasta'))
    a, b = str(recs[0].seq), str(recs[1].seq)
    name = recs[0].id.replace('_consensus_tiled', ''); c = CONS[name]; L = len(c)
    t = a.replace('-', '').upper(); off = (c + c).find(t[:40]); assert off >= 0
    tiles = []; j = 0; p = None; cur = None
    for x, y in zip(a, b):
        if x != '-':
            p = (off + j) % L; j += 1
            if cur is None or p == 0: cur = {'pos': {}, 'ins': {}}; tiles.append(cur)
            cur['pos'][p] = y.upper() if y.isupper() or y == '-' else '-'
        elif y.isupper() and cur is not None:
            cur['ins'].setdefault(p, []).append(y)
    width = {q: max((len(tl['ins'].get(q, [])) for tl in tiles), default=0) for q in range(L)}
    def row(tl):
        out = []
        for q in range(L):
            out.append(tl['pos'].get(q, '-')); ins = tl['ins'].get(q, []); out.append(''.join(ins) + '-' * (width[q] - len(ins)))
        return ''.join(out)
    cons = ''.join(c[q] + '-' * width[q] for q in range(L))
    rows = [(f'{name}_consensus', cons)] + [(f'm{k + 1:02d}', row(tl)) for k, tl in enumerate(tiles)]
    rows = [(n, s) for n, s in rows if n.endswith('consensus') or s.strip('-')]
    return name, rows, recs[1].id
if __name__ == '__main__':
    ids = [l.split('\t')[0] for l in open(sys.argv[2], encoding='utf-8') if l.strip()]
    for i in ids:
        name, rows, hdr = stack(f'aln/{i}.fa')
        with open(f'stacked/{i}.stacked.fa', 'w', encoding='utf-8', newline='\n') as f:
            for n, s in rows: f.write(f'>{n}{"" if n.endswith("consensus") else "|" + hdr.split("|")[1]}\n{s}\n')
        print(i, name, len(rows) - 1, 'monomer rows')
