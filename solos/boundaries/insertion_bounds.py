# Insertion boundaries of solo copies from the places in other assemblies where the copy is absent.
# Evidence = empty sites only (a place found by both flanks where the copy is absent). Global alignment (affine, free end
# gaps) of the source region (300 bp + copy + 300 bp, genomic orientation of the source) with the empty site region; the
# insertion = the longest run of columns where the empty site has a gap, runs separated by <= 30 bp of the empty site's
# own sequence (filler) merged. Output: per copy and target, insertion start/end in source-region coordinates (0-based) and the shift
# against the catalog copy (300 .. 300+len-1). Copies with no empty site cannot be tested this way.
import csv, sys
from Bio import SeqIO, Align
SHOW = 300; MERGE = 30
al = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-3, open_gap_score=-6, extend_gap_score=-1,
                           end_gap_score=0)
LEN = {r[0]: int(r[4]) - int(r[3]) + 1 for r in csv.reader(open(sys.argv[1], encoding='utf-8'), delimiter='\t')}
ST = {}
for r in csv.reader(open('traced/status.tsv', encoding='utf-8'), delimiter='\t'):
    ST[(r[0], r[4])] = r[7]
def blocks(q, t):
    a = al.align(q, t)[0]
    rows = list(a)
    s, tt = str(rows[0]), str(rows[1])
    n = len(s)
    # runs of columns where the empty site has a gap and the source a base
    runs = []; k = 0
    while k < n:
        if tt[k] == '-' and s[k] != '-':
            j = k
            while j < n and tt[j] == '-': j += 1
            runs.append([k, j - 1]); k = j
        else: k += 1
    if not runs: return None, None, 0
    # merge runs separated by <= MERGE target bases (the empty site's own short filler)
    merged = [runs[0]]
    for r in runs[1:]:
        between = sum(1 for k in range(merged[-1][1] + 1, r[0]) if tt[k] != '-')
        if between <= MERGE: merged[-1][1] = r[1]
        else: merged.append(r)
    best = max(merged, key=lambda r: sum(1 for k in range(r[0], r[1] + 1) if s[k] != '-'))
    srcpos = [-1] * n; p = 0
    for k in range(n):
        if s[k] != '-': srcpos[k] = p; p += 1
    ins_s = next(srcpos[k] for k in range(best[0], n) if srcpos[k] >= 0)
    ins_e = next(srcpos[k] for k in range(best[1], -1, -1) if srcpos[k] >= 0)
    filler = sum(1 for k in range(best[0], best[1] + 1) if tt[k] != '-')
    return ins_s, ins_e, filler
out = open('boundaries/insertion_bounds.tsv', 'w', encoding='utf-8', newline='\n')
out.write('id\ttarget\tstatus\tcatalog_start\tcatalog_end\tinsertion_start\tinsertion_end\tshift_start\tshift_end\tempty_site_filler_bp\n')
import os
ids = sorted({k[0] for k in ST if os.path.exists(f'traced/{k[0]}.aln.fa')})
for i in ids:
    recs = {r.description.split('|')[1]: str(r.seq).replace('-', '') for r in SeqIO.parse(f'traced/{i}.aln.fa', 'fasta')}
    src = list(SeqIO.parse(f'traced/{i}.aln.fa', 'fasta'))[0].description.split('|')[1]
    q = recs[src]
    for t, seq in recs.items():
        if t == src or ST.get((i, t)) != 'empty': continue
        s0, s1, tb = blocks(q.upper(), seq.upper())
        c0, c1 = SHOW, SHOW + LEN[i] - 1
        out.write(f'{i}\t{t}\t{ST[(i, t)]}\t{c0}\t{c1}\t{s0}\t{s1}\t{s0 - c0}\t{s1 - c1}\t{tb}\n')
        print(f'{i:11s} {t:5s} {ST[(i, t)]:15s} catalog {c0}-{c1}  insertion {s0}-{s1}  shift {s0 - c0:+d}/{s1 - c1:+d}  filler in the empty site {tb} bp')
