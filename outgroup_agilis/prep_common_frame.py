# Common frame for CLsat consensuses, the L. agilis 147-bp unit, Agi160 clones and L. agilis Agi160 array units.
from Bio import Align, SeqIO
from Bio.Seq import Seq
rc=lambda s:str(Seq(s).reverse_complement())
C={r.id:str(r.seq).upper().replace('-','') for r in SeqIO.parse('C:/work/CLsat/manuscript/subsubfamily_consensuses.fasta','fasta')}
A={r.id:str(r.seq).upper() for r in SeqIO.parse('C:/work/CLsat/review/nikitin_agi160_2026-11-10/agi160.fa','fasta')}
T=[l.rstrip('\n').split('\t') for l in open('trf_hits.tsv')]
al=Align.PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,open_gap_score=-5,extend_gap_score=-2)
def frame(q,anchor):
    """strand and rotation of q that put it in anchor's frame: returns (rotated q, strand, score, identity, alnlen)"""
    best=None
    for strand,s in (('+',q),('-',rc(q))):
        d=s+s
        a=al.align(anchor,d)[0]
        (a0,a1),(q0,q1)=a.aligned[0][0],a.aligned[1][0]
        rot=(q0-a0)%len(s)
        L=sum(1 for _ in zip(*a)); m=sum(1 for x,y in zip(*a) if x==y and x!='-')
        if best is None or a.score>best[2]: best=(s[rot:]+s[:rot],strand,a.score,m/L,L)
    return best
anc=C['CLsat3.2']
u147=[r for r in T if r[0]=='147']
U,st,sc,idn,L=frame(u147[0][6],anc); print('147 unit -> CLsat3.2',st,sc,round(idn,3),L)
out=[('CLsat'+k[5:] if False else k, v) for k,v in C.items()]
out.append(('Lagilis_147bp_unit_NC_046330.1', U))
agi_units=[]
AGA=frame(A['AY184829.1'],U)[0]          # Agi160 anchor, itself in the unit's frame
print('Agi160 anchor AY184829 framed on the 147-bp unit')
for k,v in A.items():
    s,st2,sc2,i2,L2=frame(v,AGA); agi_units.append((f'Agi160_clone_{k}_{len(v)}bp',s)); print(k,len(v),st2,sc2,round(i2,3),L2)
seen=set()
for r in T:
    if r[0]=='agi' and 140<=int(r[4])<=200:
        key=(r[1],r[2])
        s,st2,sc2,i2,L2=frame(r[6],AGA); print(r[1],r[2],r[4],sc2,round(i2,3),L2)
        agi_units.append((f'Lagilis_Agi160array_{r[1]}_{r[2]}_{r[3]}_p{r[4]}',s))
out+=agi_units
with open('units_framed.fa','w') as f:
    for k,v in out: f.write(f'>{k}\n{v}\n')
print('units',len(out))
# monomers of the two 147-bp arrays, cut in the unit's frame
with open('arrays147_monomers.fa','w') as f:
    f.write(f'>CLsat3.2_consensus\n{anc}\n>Lagilis_147bp_unit\n{U}\n')
    for r in u147:
        arr=r[7] if st=='+' else rc(r[7])
        a=al.align(U,arr[:400])[0]; (a0,_),(q0,_)=a.aligned[0][0],a.aligned[1][0]
        off=(q0-a0)%147
        if off>0: f.write(f'>{r[1]}:{r[2]}-{r[3]}_part0\n{arr[:off]}\n')
        k=1
        for i in range(off,len(arr),147):
            f.write(f'>{r[1]}:{r[2]}-{r[3]}_m{k:02d}\n{arr[i:i+147]}\n'); k+=1
print('monomer file written')
