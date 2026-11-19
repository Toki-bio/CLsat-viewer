# Stacked monomer view of a solo copy: cut the copy-vs-tiled-consensus alignment (aln/<id>.fa) at consensus position 1;
# one row per monomer (UPPERCASE), columns = consensus positions 1..L plus the insertion columns any monomer needs;
# flank rows above and below in lowercase, continued in the consensus frame. No realignment.
import sys, csv
from Bio import SeqIO
CONS = {r.id: str(r.seq).upper().replace('-', '') for r in SeqIO.parse(sys.argv[1], 'fasta')}
def stack(fn):
    recs = list(SeqIO.parse(fn, 'fasta'))
    a, b = str(recs[0].seq), str(recs[1].seq)
    name = recs[0].id.replace('_consensus_tiled', ''); c = CONS[name]; L = len(c)
    t = a.replace('-', '').upper(); off = (c + c).find(t[:40]); assert off >= 0
    tiles = []; j = 0; p = None; cur = None; lead = []
    for x, y in zip(a, b):
        if x == '-' and y.isupper() and cur is None: lead.append(y); continue   # copy bases before the first consensus base
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
    # flanks (lowercase in aln/<id>.fa) continued in the consensus frame without alignment: the base just before the
    # copy sits in the column before the copy's first base; two rows on each side, written lowercase
    core = [k for k, y in enumerate(b) if y.isupper()]
    fL = ''.join(y for y in b[:core[0]] if y != '-').lower(); fR = ''.join(y for y in b[core[-1] + 1:] if y != '-').lower()
    pos_core = [((off + sum(1 for x in a[:k] if x != '-')) % L) for k in (core[0], core[-1])]
    p0, pe = pos_core
    def frow(d):
        return ''.join(d.get(q, '-') + '-' * width[q] for q in range(L))
    left = {}; j = 0
    for q_abs in range(p0 - 1 - len(lead), p0 - 1 - len(lead) - 2 * L, -1):   # left flank ends before the unaligned start, if any
        if j >= len(fL): break
        left.setdefault(q_abs // L, {})[q_abs % L] = fL[-1 - j]; j += 1
    right = {}; j = 0
    for q_abs in range(pe + 1, pe + 1 + 2 * L):
        if j >= len(fR): break
        right.setdefault(q_abs // L, {})[q_abs % L] = fR[j]; j += 1
    lrows = [(f'left_flank_{-k}', frow(left[k])) for k in sorted(left)]
    rrows = [(f'right_flank_{k + 1}', frow(right[k])) for k in sorted(right)]
    lead_rows = []
    if lead:   # not aligned: placed like the flanks, its last base in the column before the copy's first aligned base
        d = {}
        for j2, ch in enumerate(reversed(lead)): d[(p0 - 1 - j2) % L] = ch
        lead_rows = [(f'm00_start_not_aligned_{len(lead)}bp', frow(d))]
    rows = [(f'{name}_consensus', cons)] + lrows + lead_rows + [(f'm{k + 1:02d}', row(tl)) for k, tl in enumerate(tiles)] + rrows
    rows = [(n, s) for n, s in rows if n.endswith('consensus') or s.strip('-')]
    return name, rows, recs[1].id
if __name__ == '__main__':
    ids = [l.split('\t')[0] for l in open(sys.argv[2], encoding='utf-8') if l.strip()]
    for i in ids:
        name, rows, hdr = stack(f'aln/{i}.fa')
        with open(f'stacked/{i}.stacked.fa', 'w', encoding='utf-8', newline='\n') as f:
            for n, s in rows: f.write(f'>{n}{"" if n.endswith("consensus") else "|" + hdr.split("|")[1]}\n{s}\n')
        print(i, name, len(rows) - 1, 'monomer rows')
