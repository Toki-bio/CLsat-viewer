# BED tracks (chrom "alignment" = alignment columns, 0-based, end exclusive) for the solo alignments in MSA-viewer:
#   aln/<id>.bed     copy vs tiled consensus: left flank, monomers m01.. (cut at consensus position 1), right flank,
#                    insertion from empty sites where tested
#   traced/<id>.bed  copy across assemblies: the same features placed on the source row
# Monomers: each core base of the copy-vs-consensus alignment belongs to the consensus repeat it is aligned in; a new
# monomer starts where the consensus position returns to 1. Insertion: boundaries/insertion_bounds.tsv (empty sites).
import csv, os, sys, glob
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
CONS = {r.id: str(r.seq).upper().replace('-', '') for r in SeqIO.parse(sys.argv[1], 'fasta')}
INS = {}
for r in csv.DictReader(open('boundaries/insertion_bounds.tsv', encoding='utf-8'), delimiter='\t'):
    INS.setdefault(r['id'], []).append((int(r['insertion_start']), int(r['insertion_end'])))
INS = {k: Counter(v).most_common(1)[0][0] for k, v in INS.items()}     # targets agree; take the common value
COL = {'flank': '160,160,160', 'm1': '30,136,229', 'm2': '0,150,136', 'part': '255,179,0', 'ins': '211,47,47'}
def monomers(fn):
    recs = list(SeqIO.parse(fn, 'fasta')); a, b = str(recs[0].seq), str(recs[1].seq)
    name = recs[0].id.replace('_consensus_tiled', ''); c = CONS[name]; L = len(c)
    off = (c + c).find(a.replace('-', '').upper()[:40])
    res_of_col = {}; r = 0
    for k, y in enumerate(b):
        if y != '-': res_of_col[k] = r; r += 1
    N = r
    tiles = []; j = 0; cur = None; p = None; seen = []; lead = []
    for k, (x, y) in enumerate(zip(a, b)):
        if x == '-' and y.isupper() and cur is None: lead.append(res_of_col[k]); continue   # copy bases before the first consensus base
        if x != '-':
            p = (off + j) % L; j += 1
            if cur is None or p == 0: cur = []; tiles.append(cur); seen.append(set())
            if y.isupper(): cur.append(res_of_col[k]); seen[-1].add(p)
        elif y.isupper() and cur is not None: cur.append(res_of_col[k])
    core = [res_of_col[k] for k, y in enumerate(b) if y.isupper()]
    mons = [(min(t), max(t), len(s) / L) for t, s in zip(tiles, seen) if t]
    return name, L, N, core[0], core[-1], mons, lead
def feats(i, fn, strand_rc):
    name, L, N, c0, c1, mons, lead = monomers(fn)
    F = [('left flank', 0, c0 - 1, COL['flank']), ('right flank', c1 + 1, N - 1, COL['flank'])]
    if lead: F.append((f'start not aligned to {name} ({len(lead)} bp)', min(lead), max(lead), COL['part']))
    for k, (s, e, cov) in enumerate(mons):
        col = COL['part'] if cov < 0.9 else (COL['m1'] if k % 2 == 0 else COL['m2'])
        F.append((f'm{k + 1:02d} {name} {round(cov * L)}/{L}', s, e, col))
    if i in INS:
        s, e = INS[i]                                    # genomic-orientation region coordinates
        if strand_rc: s, e = N - 1 - e, N - 1 - s       # to the copy-vs-consensus orientation
        F.append(('insertion (empty sites)', s, e, COL['ins']))
    return F, N
def write_bed(path, F, col_of_res, title):
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(f'track name="{title}"\n')
        for nm, s, e, col in sorted(F, key=lambda x: x[1]):
            if e < s: continue
            cs, ce = col_of_res[s], col_of_res[e]
            f.write(f'alignment\t{cs}\t{ce + 1}\t{nm}\t0\t.\t{cs}\t{ce + 1}\t{col}\n')
n_aln = n_tr = 0
for fn in sorted(glob.glob('aln/*.fa')):
    i = os.path.basename(fn)[:-3]
    b = str(list(SeqIO.parse(fn, 'fasta'))[1].seq)
    col_of_res = [k for k, y in enumerate(b) if y != '-']
    tr = f'traced/{i}.aln.fa'; rc = False
    if os.path.exists(tr):
        t = str(list(SeqIO.parse(tr, 'fasta'))[0].seq)
        rc = b.replace('-', '').upper() != t.replace('-', '').upper()
    F, N = feats(i, fn, rc)
    write_bed(f'aln/{i}.bed', F, col_of_res, f'{i} monomers'); n_aln += 1
    if os.path.exists(tr):
        recs = list(SeqIO.parse(tr, 'fasta')); src = str(recs[0].seq)
        tcol = [k for k, y in enumerate(src) if y != '-']
        G = []
        for nm, s, e, col in F:                          # map copy-vs-consensus residues to the source row
            if rc:
                s, e = N - 1 - e, N - 1 - s
                nm = {'left flank': 'right flank', 'right flank': 'left flank'}.get(nm, nm)   # flanks named in genomic order
            if i in INS:                                 # corrected copy: flanks from the insertion, features clipped to it
                if nm in ('left flank', 'right flank'): continue
                s, e = max(s, INS[i][0]), min(e, INS[i][1])
                if e < s: continue
            G.append((nm, s, e, col))
        if i in INS:
            G += [('left flank', 0, INS[i][0] - 1, COL['flank']), ('right flank', INS[i][1] + 1, len(tcol) - 1, COL['flank'])]
        write_bed(f'traced/{i}.bed', G, tcol, f'{i} monomers')
        if i in INS:                                     # uppercase = insertion from the empty sites
            s, e = INS[i]; up = set(range(tcol[s], tcol[e] + 1))   # whole span, also columns where only other rows have bases
            with open(tr, 'w', encoding='utf-8', newline='\n') as f:
                for r in recs:
                    q = ''.join(ch.upper() if k in up else ch.lower() for k, ch in enumerate(str(r.seq)))
                    f.write(f'>{r.description}\n{q}\n')
        n_tr += 1
print('aln beds', n_aln, 'traced beds', n_tr, 'insertion-corrected', sorted(INS))
