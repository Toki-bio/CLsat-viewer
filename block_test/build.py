import json, os, csv
H = os.path.dirname(os.path.abspath(__file__))
# all called blocks: half (unm = unv, unp = unn), contig, start, nai_diff, val_diff, call
B = []
for r in csv.reader(open(f'{H}/B_blocks_called.tsv', encoding='utf-8'), delimiter='\t'):
    if r[0] == 'half': continue
    B.append([0 if r[0] == 'unm' else 1, round(float(r[4]), 5), round(float(r[3]), 5), {'valentini fits better': 'v', 'nairensis fits better': 'n'}.get(r[5], 'c')])
# switch contig: reads per block (MAPQ>=30) and difference rate, both animals
C = 'JAWDKN010000070.1'
T = {}
for an, f in (('nai', 'B_nai_17_78_on_unp.tsv'), ('val', 'B_val_5_38_on_unp.tsv')):
    for r in csv.reader(open(f'{H}/{f}', encoding='utf-8'), delimiter='\t'):
        if r[2] != C: continue
        s = int(r[3]); T.setdefault(s, {})[an] = [int(r[4]), int(r[5]), int(r[6])]
calls = {int(r[2]): r[5] for r in csv.reader(open(f'{H}/B_blocks_called.tsv', encoding='utf-8'), delimiter='\t') if r[1] == C}
track = []
for s in range(0, 23112228, 100000):
    d = T.get(s, {}); n = d.get('nai', [0, 0, 0]); v = d.get('val', [0, 0, 0])
    track.append([s, n[0], round(n[2] / n[1], 5) if n[1] else None, v[0], round(v[2] / v[1], 5) if v[1] else None,
                  {'valentini fits better': 'v', 'nairensis fits better': 'n', 'close': 'c'}.get(calls.get(s), '')])
reads = None
if os.path.exists(f'{H}/reads.tsv'):
    WIN = {'nai_part': dict(s=18499250, e=18500750, title='Window 1: 1.5 kb in the D. raddei nairensis part (18.50 Mb)'),
           'val_part': dict(s=19499250, e=19500750, title='Window 2: 1.5 kb in the D. valentini part after the switch (19.50 Mb)'),
           'L17_start': dict(s=22776000, e=22780000, arr=22777927, title='Window 3, why arrays get no call: 4 kb at the start of the CLsat1.1 array L17 (shaded)')}
    for k in WIN: WIN[k]['reads'] = []
    ALL = {}
    FULL = {'nai_part': '10 kb', 'val_part': '10 kb', 'L17_start': '20 kb'}
    for r in csv.reader(open(f'{H}/reads.tsv', encoding='utf-8'), delimiter='\t'):
        w, a, p, e, q, m = r[0], r[1], int(r[2]), int(r[3]), int(r[4]), (r[5] if len(r) > 5 else '')
        rec = dict(a=a, p=p, e=e, q=q, m=[int(x) for x in m.split(',') if x])
        ALL.setdefault(w, []).append(rec)                 # whole 10/20-kb window: for the numbers
        if e < WIN[w]['s'] or p > WIN[w]['e']: continue   # zoomed window: for the drawing
        WIN[w]['reads'].append(rec)
    for k, w in WIN.items():
        st = {}
        for a in ('nai', 'val'):
            rs = [x for x in ALL[k] if x['a'] == a and x['q'] >= 30]
            bases = sum(x['e'] - x['p'] + 1 for x in rs); dif = sum(len(x['m']) for x in rs)
            st[a] = (len(rs), len([x for x in ALL[k] if x['a'] == a]), dif / bases if bases else 0)
        w['cap'] = (f"In the surrounding {FULL[k]}: <i>D. raddei nairensis</i> reads: {st['nai'][1]} ({st['nai'][0]} with mapping quality ≥ 30), "
                    f"{100*st['nai'][2]:.1f}% of their bases differ; <i>D. valentini</i> reads: {st['val'][1]} ({st['val'][0]} with "
                    f"mapping quality ≥ 30), {100*st['val'][2]:.1f}% differ.")
        if k == 'L17_start': w['cap'] = 'Not a comparison of the two animals. Every read was mapped to the whole assembly and went to its best place, so the reads shown here are those that chose this place; a read of a CLsat type found elsewhere in the assembly goes there instead. What the window shows: the 400 bp before the array are themselves a masked repeat, and no read of either animal starts in unique flank sequence and continues into the array. Reads made of array alone stack in columns one monomer apart, mostly with mapping quality 0, because they fit every monomer equally. So the block test cannot reach the array, and its parental side is judged from the called blocks next to it (22.6 Mb: <i>D. valentini</i>).'
        w['title'] = w['title'].replace('D. raddei nairensis', '<i>D. raddei nairensis</i>').replace('D. valentini', '<i>D. valentini</i>')
        print(k, st)
    reads = dict(windows=[WIN[k] for k in ('nai_part', 'val_part', 'L17_start')],
                 note='Each bar is one read where it maps on the unn contig, upper group <i>D. raddei nairensis</i>, lower group '
                      '<i>D. valentini</i>; dark ticks mark bases that differ from the assembly (mismatches, insertion and deletion '
                      'sites). Faded bars are reads with mapping quality below 30, which the block test does not count. A read that '
                      'lies wholly inside the array cannot be placed: it fits equally well at many monomers.')
ARM = []   # D. armeniaca: KIT clsat_map/armsplit/full/called_*.tsv; col4 = D. mixta reads, col5 = D. valentini reads
for k, f in enumerate(('called_val_38_75_mix_DM2.tsv', 'called_val_33_75_mix_DM3.tsv')):
    for r in csv.reader(open(f'{H}/arm/{f}', encoding='utf-8'), delimiter='\t'):
        ARM.append([0 if r[0] == 'armf' else 1, k, round(float(r[4]), 5), round(float(r[3]), 5), {'valentini': 'v', 'mixta': 'm'}.get(r[5], 'c')])
# how common are switches: blocks called for the other parent, contigs with a change of call, CLsat loci in them
import collections
LM = 'C:/work/CLsat/review/results_revision/run_gb_v2/regions/locus_master.tsv'
LOCI = collections.defaultdict(list)
for r in csv.reader(open(LM, encoding='utf-8'), delimiter='\t'): LOCI[r[4]].append((r[1], int(r[2]), int(r[3]), r[0]))
def _calls(f, names):
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    for r in csv.reader(open(f, encoding='utf-8'), delimiter='\t'):
        if r[0] == 'half': continue
        c = 'v' if 'valentini' in r[5] else 'n' if 'nairensis' in r[5] else 'm' if 'mixta' in r[5] else 'c'
        out[names[r[0]]][r[1]][int(r[2])] = c
    return out
EXP = {'unv': 'v', 'unn': 'n', 'armf': 'v', 'armm': 'm'}
SW = []
for f, names in ((f'{H}/B_blocks_called.tsv', {'unm': 'unv', 'unp': 'unn'}), (f'{H}/arm/called_val_38_75_mix_DM2.tsv', {'armf': 'armf', 'armm': 'armm'})):
    for asm, ctg in _calls(f, names).items():
        e = EXP[asm]; wrong = []; nctg = 0; tot = 0
        for c, bl in ctg.items():
            tot += len(bl); w = [s0 for s0, x in bl.items() if x not in (e, 'c')]
            if w: nctg += 1; wrong += [(c, s0) for s0 in w]
        loci = sorted({l[3] for l in LOCI[asm] for c, s0 in wrong if l[0] == c and l[2] > s0 - 100000 and l[1] < s0 + 200000})
        SW.append(dict(asm=asm, tot=tot, wrong=len(wrong), ctg=nctg, loci=loci))
data = dict(sw=SW, arm=ARM, blocks=B, track=track, contig=C, contig_len=23112228, array=[22777927, 23112228], reads=reads)
tpl = open(f'{H}/page.tpl.html', encoding='utf-8').read()
open(f'{H}/short_read_blocks.html', 'w', encoding='utf-8').write(tpl.replace('/*DATA*/null', json.dumps(data, separators=(',', ':'))))
print(len(B), 'blocks', len(track), 'track rows', 'reads' if reads else 'no reads yet')
