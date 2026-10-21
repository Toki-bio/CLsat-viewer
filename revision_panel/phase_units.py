"""Phase cluster units to their family consensus before aligning: strand and start from the best ssearch36 hit
against the doubled consensus; the unit position aligned to consensus position 1 is read
from the alignment itself (BTOP), so gaps between the hit start and that position do not shift the start; every unit is turned to the consensus strand and rotated to begin at consensus position 1.
Usage: phase_units.py cluster.fa consensus.fa out.fa"""
import subprocess, sys, os, tempfile
fa, cons_fa, out = sys.argv[1:4]
cons = ''.join(l.strip() for l in open(cons_fa) if not l.startswith('>')).upper(); L = len(cons)
tmp = tempfile.mkdtemp(); x2 = os.path.join(tmp, 'x2.fa')
open(x2, 'w').write('>c\n' + cons + cons + '\n')
seqs = {}; order = []; n = None
for l in open(fa):
    if l.startswith('>'): n = l[1:].strip(); seqs[n] = ''; order.append(n)
    else: seqs[n] += l.strip().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(comp)[::-1]
import re
def best_hit(s, extra):
    q = os.path.join(tmp, 'q.fa'); open(q, 'w').write('>q\n' + s + '\n')
    r = subprocess.run(['ssearch36', '-n', '-m', '8CB', '-E', '10', '-b', '1', '-d', '1'] + extra + [q, x2],
                       capture_output=True, text=True).stdout
    return next((l.split('\t') for l in r.splitlines() if l and not l.startswith('#')), None)
def start_at_cons1(qs, ss, btop):
    """walk the alignment (BTOP) from the hit start; return the 0-based unit position aligned to consensus position 1"""
    qp, sp = qs, ss
    for tok in re.findall(r'\d+|[A-Z*-]{2}', btop.upper()):
        cols = [('x', 'x')] * int(tok) if tok.isdigit() else [(tok[0], tok[1])]
        for a, b in cols:
            if b != '-' and (sp - 1) % L == 0: return qp - 1
            qp += a != '-'; sp += b != '-'
    return None
o = []; walked = 0
for k in order:
    s = seqs[k]
    hit = best_hit(s, [])
    if not hit: o.append((k + '|unphased', s)); continue
    if int(hit[6]) > int(hit[7]):     # minus strand: turn the unit, then align again on the plus strand only
        s = rc(s); hit = best_hit(s, ['-3'])
    qs, ss = int(hit[6]), int(hit[8])
    start = start_at_cons1(qs, ss, hit[12])
    if start is None:                 # the hit never reaches consensus position 1: count back from the hit start (no gaps assumed)
        start = (qs - 1 - (ss - 1) % L) % len(s)
    else: walked += 1
    o.append((k, s[start % len(s):] + s[:start % len(s)]))
print('%d units, %d phased by the alignment walk, %d by counting back' % (len(o), walked, len(o) - walked), file=sys.stderr)
with open(out, 'w') as f:
    f.write('>consensus\n' + cons + '\n')
    for k, s in o: f.write('>' + k + '\n' + s + '\n')
