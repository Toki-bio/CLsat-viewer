#!/usr/bin/env python3
"""
Automated test suite for CLsat-viewer MSA analysis.
Port of the core JS functions to Python, fetches sequences from NCBI,
runs the full analysis + rendering pipeline, and validates for common bugs.

Usage: python test_msa.py
"""

import sys, time, urllib.request, json, math, re
import numpy as np
from typing import List, Dict, Any, Tuple, Optional

# ============================================================================
# CONSENSUS SEQUENCES (same as index.html)
# ============================================================================
CONSENSUS_SEQS = {
    "CLsat1": "AAAGAAGCTAGTTNGGATGTATTTGGCTTTTCTAACGTTCAGTTTGGCTTACTTGCGTGATTTTGCATTAACAGCTCAAATACAGCTAACTTTGGAAATGAACACGATGGAAACTTGTTGGTGTGTTTTCTATGCATTTCGACCTG",
    "CLsat2": "AAAGAAGCTCATTGGGATGTAGTTGTGTTTCTAAGCTTCATTTTAGCTTATTTGCGTGATTTCGCTTTTAAGGCTCAAATACAGCTATCTTTTGAAAGCAACACGTTAGAAACTTCTTGGTGTGTTTTTCATGCATTTGGACCCA",
    "CLsat3": "AAAGAAGCTCGTTCGGATGTCTTTCTGTTTTCTAAGCTTCATTTGAGCTTCTTTGCGTCAATTTGCACTCCAGACTCAAAAACGGCTCTCTTTTGAAAGGAACACCAAAGAAGCTTGTTCGTGTGTTTCCCATGCATTTCCAGCCC",
    "CLsat4": "AAAGAAGCTCGTTTGGATGTATTTGCGTGACTTTGCATTCACAACTCAAAAACACATATCTTTTGAAATGAAGGTGATGGAAGCTTGTTCGTGTGTTTTCTATGCATTTCGAGCCT",
}

# ============================================================================
# CORE ANALYSIS FUNCTIONS (faithful port from index.html JS)
# ============================================================================

def rev_comp(seq: str) -> str:
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N',
            'a':'t','t':'a','g':'c','c':'g','n':'n'}
    return ''.join(comp.get(c, c) for c in reversed(seq))

def smith_waterman(query: str, subject: str, match=2, mismatch=-2, gap_open=-13, gap_ext=-3) -> list:
    qLen = len(query)
    sLen = len(subject)
    qU = query.upper()
    sU = subject.upper()

    # Allocate numpy matrices
    H = np.zeros((qLen + 1, sLen + 1), dtype=np.float32)
    E = np.zeros((qLen + 1, sLen + 1), dtype=np.float32)
    F = np.zeros((qLen + 1, sLen + 1), dtype=np.float32)
    traceH = np.zeros((qLen + 1, sLen + 1), dtype=np.int8)

    # Pre-compute score row for each query position
    q_arr = np.array([ord(c) for c in qU], dtype=np.int32)
    s_arr = np.array([ord(c) for c in sU], dtype=np.int32)
    N_ord = ord('N')

    # Fill matrices row by row (vectorized over columns)
    for i in range(1, qLen + 1):
        # Score vector: match if same char and not N
        sc = np.where((s_arr == q_arr[i-1]) & (q_arr[i-1] != N_ord), match, mismatch)

        E[i, 1:] = np.maximum(H[i, :-1] + gap_open + gap_ext, E[i, :-1] + gap_ext)
        # F needs sequential dependency on rows, but we can compute it
        F[i, 1:] = np.maximum(H[i-1, 1:] + gap_open + gap_ext, F[i-1, 1:] + gap_ext)

        diag = H[i-1, :-1] + sc
        best = np.maximum(0, np.maximum(diag, np.maximum(E[i, 1:], F[i, 1:])))
        H[i, 1:] = best

        # Trace: 0=stop, 1=diag, 2=up(F), 3=left(E)
        tr = np.zeros(sLen, dtype=np.int8)
        tr[best == diag] = 1
        tr[best == F[i, 1:]] = 2  # may overwrite diag ties
        tr[best == E[i, 1:]] = 3  # may overwrite F ties
        tr[best == 0] = 0
        # Priority: diag > F > E (same as JS)
        tr = np.where(best == 0, 0,
             np.where(best == diag, 1,
             np.where(best == F[i, 1:], 2, 3)))
        traceH[i, 1:] = tr

    # Note: E has col-wise dependency that numpy handles correctly row-by-row
    # But E[i,j] depends on E[i,j-1], so the vectorized version above is wrong for E.
    # We need to do E column by column within each row. Let's fix this with a loop.
    # Actually, let's redo the fill more carefully.
    H[:] = 0; E[:] = 0; F[:] = 0; traceH[:] = 0

    for i in range(1, qLen + 1):
        qi = q_arr[i-1]
        for j in range(1, sLen + 1):
            sc = match if (s_arr[j-1] == qi and qi != N_ord) else mismatch
            e_val = max(H[i, j-1] + gap_open + gap_ext, E[i, j-1] + gap_ext)
            f_val = max(H[i-1, j] + gap_open + gap_ext, F[i-1, j] + gap_ext)
            E[i, j] = e_val
            F[i, j] = f_val
            diag = H[i-1, j-1] + sc
            best = max(0.0, diag, e_val, f_val)
            H[i, j] = best
            if best == 0: traceH[i, j] = 0
            elif best == diag: traceH[i, j] = 1
            elif best == f_val: traceH[i, j] = 2
            else: traceH[i, j] = 3

    threshold = max(qLen * match * 0.35, 20)
    # Create a masked copy for multi-hit extraction
    H_masked = H.copy()
    hits = []

    for _pass in range(5000):
        # Find max in masked matrix (numpy is fast for this)
        flat_idx = np.argmax(H_masked[1:, 1:])
        max_i = flat_idx // sLen + 1
        max_j = flat_idx % sLen + 1
        max_val = H_masked[max_i, max_j]
        if max_val < threshold:
            break

        i, j = int(max_i), int(max_j)
        align_len, matches = 0, 0
        end_j = j
        while i > 0 and j > 0 and H[i, j] > 0:
            t = traceH[i, j]
            if t == 1:
                if qU[i-1] == sU[j-1] and qU[i-1] != 'N': matches += 1
                align_len += 1; i -= 1; j -= 1
            elif t == 2: align_len += 1; i -= 1
            elif t == 3: align_len += 1; j -= 1
            else: break
        start_j = j
        identity = matches / align_len if align_len > 0 else 0
        hits.append({
            'start': start_j, 'end': end_j - 1,
            'score': float(max_val), 'identity': identity,
            'length': end_j - start_j, 'alignLen': align_len, 'matches': matches,
        })
        # Mask used columns for all rows
        H_masked[1:, start_j+1:end_j+1] = 0

    hits.sort(key=lambda h: h['start'])
    return hits

def pairwise_align(query: str, subject: str):
    match, mismatch, gap_open, gap_ext = 2, -2, -13, -3
    qLen, sLen = len(query), len(subject)
    qU, sU = query.upper(), subject.upper()
    H = np.zeros((qLen+1, sLen+1), dtype=np.float32)
    E = np.zeros((qLen+1, sLen+1), dtype=np.float32)
    F = np.zeros((qLen+1, sLen+1), dtype=np.float32)
    tr = np.zeros((qLen+1, sLen+1), dtype=np.int8)
    q_arr = np.array([ord(c) for c in qU], dtype=np.int32)
    s_arr = np.array([ord(c) for c in sU], dtype=np.int32)
    N_ord = ord('N')
    max_val, max_i, max_j = 0.0, 0, 0
    for i in range(1, qLen+1):
        qi = q_arr[i-1]
        for j in range(1, sLen+1):
            sc = match if (s_arr[j-1] == qi and qi != N_ord) else mismatch
            e_val = max(H[i, j-1] + gap_open + gap_ext, E[i, j-1] + gap_ext)
            f_val = max(H[i-1, j] + gap_open + gap_ext, F[i-1, j] + gap_ext)
            E[i, j] = e_val
            F[i, j] = f_val
            diag = H[i-1, j-1] + sc
            best = max(0.0, diag, e_val, f_val)
            H[i, j] = best
            if best == 0: tr[i, j] = 0
            elif best == diag: tr[i, j] = 1
            elif best == f_val: tr[i, j] = 2
            else: tr[i, j] = 3
            if best > max_val: max_val, max_i, max_j = best, i, j
    qA, sA = [], []
    ci, cj, matches, align_len = max_i, max_j, 0, 0
    while ci>0 and cj>0 and H[ci, cj]>0:
        t = tr[ci, cj]
        if t==1:
            qA.append(qU[ci-1]); sA.append(sU[cj-1])
            if qU[ci-1]==sU[cj-1] and qU[ci-1]!='N': matches+=1
            align_len+=1; ci-=1; cj-=1
        elif t==2: qA.append(qU[ci-1]); sA.append('-'); align_len+=1; ci-=1
        elif t==3: qA.append('-'); sA.append(sU[cj-1]); align_len+=1; cj-=1
        else: break
    qA.reverse(); sA.reverse()
    return {
        'qAln':''.join(qA), 'sAln':''.join(sA),
        'score':float(max_val), 'identity':matches/align_len if align_len>0 else 0,
        'matches':matches, 'alignLen':align_len,
        'sStart':cj, 'sEnd':max_j, 'qStart':ci, 'qEnd':max_i
    }

def split_wide_hits(hits, fallback_consensus, raw_seq, max_hit_len):
    if not hits: return hits
    out = []
    for m in hits:
        span = m['end'] - m['start'] + 1
        if span <= max_hit_len:
            out.append(m); continue
        cons = CONSENSUS_SEQS.get(m.get('clsat', ''), fallback_consensus) or fallback_consensus
        cons_len = len(cons)
        cons_u = cons.upper()
        fwd_slice = raw_seq[m['start']:m['end']+1]
        read_seq = rev_comp(fwd_slice) if m.get('strand') == '-' else fwd_slice
        off = 0; split_ok = False
        for si in range(20):
            if off >= len(read_seq): break
            rem = read_seq[off:]
            if len(rem) < cons_len * 0.25: break
            wa = pairwise_align(cons_u, rem)
            if wa['alignLen'] < cons_len * 0.25 or wa['identity'] < 0.35: break
            read_start = off + wa['sStart']
            read_end = off + wa['sEnd'] - 1
            if m.get('strand') == '-':
                raw_start = m['start'] + (span - 1 - read_end)
                raw_end = m['start'] + (span - 1 - read_start)
            else:
                raw_start = m['start'] + read_start
                raw_end = m['start'] + read_end
            out.append({**m, 'start':raw_start, 'end':raw_end, 'length':raw_end-raw_start+1,
                        'identity':wa['identity'], 'matches':wa['matches'], 'alignLen':wa['alignLen']})
            split_ok = True
            off += wa['sEnd']
        if not split_ok: out.append(m)
    out.sort(key=lambda h: h['start'])
    return out

def merge_hits(hits, max_merge_len=float('inf')):
    if not hits: return hits
    s = sorted(hits, key=lambda h: h['start'])
    merged = [dict(s[0])]
    for i in range(1, len(s)):
        prev = merged[-1]; cur = s[i]
        same_type = prev.get('clsat','?') == cur.get('clsat','?')
        if same_type and cur['start'] <= prev['end'] + 1:
            new_end = max(prev['end'], cur['end'])
            new_len = new_end - prev['start'] + 1
            if new_len <= max_merge_len:
                prev['end'] = new_end; prev['length'] = new_len
                prev['score'] = max(prev['score'], cur['score'])
                prev['identity'] = max(prev['identity'], cur['identity'])
            else: merged.append(dict(cur))
        else: merged.append(dict(cur))
    return merged

def split_into_arrays(hits, min_tandems=3):
    if not hits: return []
    s = sorted(hits, key=lambda h: h['start'])
    arrays = [[s[0]]]
    for i in range(1, len(s)):
        gap = s[i]['start'] - s[i-1]['end'] - 1
        if gap > 100: arrays.append([s[i]])
        else: arrays[-1].append(s[i])
    return [a for a in arrays if len(a) >= min_tandems]

def split_arrays_at_transitions(arrays):
    result = []
    for arr in arrays:
        hits = arr['hits']
        if len(hits) <= 1: result.append(arr); continue
        types = [h.get('clsat','?') for h in hits]
        if all(t == types[0] for t in types): result.append(arr); continue
        sub = [{'strand':arr['strand'], 'hits':[hits[0]]}]
        for i in range(1, len(hits)):
            if hits[i].get('clsat','?') != hits[i-1].get('clsat','?'):
                sub.append({'strand':arr['strand'], 'hits':[hits[i]]})
            else: sub[-1]['hits'].append(hits[i])
        for sa in sub:
            if len(sa['hits']) >= 1: result.append(sa)
    return result

def scan_one_consensus(seq, cons_name, min_identity):
    consensus = CONSENSUS_SEQS.get(cons_name)
    if not consensus: return []
    seq_len = len(seq)
    cons_len = len(consensus)
    Slf = 0.9
    min_len = int(Slf * cons_len)
    fwd = smith_waterman(consensus, seq)
    for h in fwd: h['strand'] = '+'; h['clsat'] = cons_name
    rc_seq = rev_comp(seq)
    rc_raw = smith_waterman(consensus, rc_seq)
    rc_hits = []
    for h in rc_raw:
        ns = seq_len - h['end'] - 1
        ne = seq_len - h['start'] - 1
        rc_hits.append({**h, 'start':ns, 'end':ne, 'length':ne-ns+1, 'strand':'-', 'clsat':cons_name})
    combined = [h for h in fwd + rc_hits if h['identity'] >= min_identity and h['length'] >= min_len]
    max_span = int(cons_len * 1.3)
    combined = split_wide_hits(combined, consensus, seq, max_span)
    return combined

def resolve_overlaps(hits):
    if not hits: return hits
    hits.sort(key=lambda h: (h['start'], -h['score']))
    kept = [hits[0]]
    for i in range(1, len(hits)):
        prev = kept[-1]; cur = hits[i]
        ov_start = max(prev['start'], cur['start'])
        ov_end = min(prev['end'], cur['end'])
        ov = max(0, ov_end - ov_start + 1)
        shorter = min(prev['length'], cur['length'])
        if ov > shorter * 0.5:
            if cur['score'] > prev['score']: kept[-1] = cur
        else: kept.append(cur)
    return kept

def analyze_monomers_sw(sequence, consensus_name, min_identity=0.75):
    seq = re.sub(r'^>.*\n', '', sequence, flags=re.MULTILINE).replace('\n','').replace(' ','').upper()
    seq_len = len(seq)
    if seq_len > 200000:
        return {'success': False, 'message': 'Too large'}
    scan_names = list(CONSENSUS_SEQS.keys())
    all_hits = []
    for cn in scan_names:
        hits = scan_one_consensus(seq, cn, min_identity)
        all_hits.extend(hits)
    all_hits = resolve_overlaps(all_hits)
    all_hits.sort(key=lambda h: h['start'])
    if not all_hits:
        return {'success':False, 'message':'No monomers', 'monomer_count':0, 'monomers':[], 'arrays':[]}
    max_cons_len = max(len(CONSENSUS_SEQS[cn]) for cn in scan_names)
    max_merge_len = int(max_cons_len * 1.3)
    plus_hits = merge_hits([h for h in all_hits if h['strand']=='+'], max_merge_len)
    minus_hits = merge_hits([h for h in all_hits if h['strand']=='-'], max_merge_len)
    primary_cons = CONSENSUS_SEQS[scan_names[0]]
    plus_hits = split_wide_hits(plus_hits, primary_cons, seq, max_merge_len)
    minus_hits = split_wide_hits(minus_hits, primary_cons, seq, max_merge_len)
    plus_arrays = split_into_arrays(plus_hits, 3)
    minus_arrays = split_into_arrays(minus_hits, 3)
    all_arrays = sorted(
        [{'strand':'+','hits':a} for a in plus_arrays] + [{'strand':'-','hits':a} for a in minus_arrays],
        key=lambda x: x['hits'][0]['start']
    )
    all_arrays = split_arrays_at_transitions(all_arrays)
    all_monomers = []
    for a in all_arrays:
        for h in a['hits']:
            all_monomers.append({**h, 'strand':a['strand']})
    avg_id = sum(m['identity'] for m in all_monomers)/len(all_monomers) if all_monomers else 0
    used_clsats = list(set(m['clsat'] for m in all_monomers))
    return {
        'success': True,
        'monomer_count': len(all_monomers),
        'array_count': len(all_arrays),
        'consensus_name': '+'.join(used_clsats) if len(used_clsats)>1 else used_clsats[0],
        'consensus_names': used_clsats,
        'consensus_length': len(CONSENSUS_SEQS.get(used_clsats[0],'')) if used_clsats else 146,
        'avg_identity': avg_id,
        'monomers': all_monomers,
        'arrays': all_arrays,
    }


# ============================================================================
# SIMULATED renderMonomerMSA — returns structured data instead of HTML
# ============================================================================
def simulate_render_msa(result, raw_seq, genomic_start):
    if not result or not result.get('success') or not result.get('arrays'):
        return []
    g_off = genomic_start or 0
    FLANK_SIZE = 300
    global_idx = 0
    array_results = []

    for ai, arr in enumerate(result['arrays']):
        strand = arr['strand']
        first_start = arr['hits'][0]['start']
        last_end = arr['hits'][-1]['end']
        clsat_counts = {}
        for m in arr['hits']:
            c = m.get('clsat','?')
            clsat_counts[c] = clsat_counts.get(c,0)+1
        dominant_clsat = max(clsat_counts, key=clsat_counts.get)
        consensus = CONSENSUS_SEQS.get(dominant_clsat)
        if not consensus: continue
        cons_u = consensus.upper()
        cons_len = len(cons_u)

        left_flank_start = max(0, first_start - FLANK_SIZE)
        left_flank = raw_seq[left_flank_start:first_start]
        right_flank_end = min(len(raw_seq), last_end + 1 + FLANK_SIZE)
        right_flank = raw_seq[last_end+1:right_flank_end]

        if strand == '-':
            flank_before = rev_comp(right_flank)
            flank_after = rev_comp(left_flank)
        else:
            flank_before = left_flank
            flank_after = right_flank
        n_before_orig = len(flank_before)
        n_after_orig = len(flank_after)

        max_hit_len = int(cons_len * 1.3)
        display_hits = split_wide_hits(list(arr['hits']), consensus, raw_seq, max_hit_len)
        # Key fix: reverse AFTER splitWideHits for minus strand
        if strand == '-':
            display_hits.reverse()

        proj_rows = []
        for mi, m in enumerate(display_hits):
            global_idx += 1
            mon_seq = raw_seq[m['start']:m['end']+1]
            if strand == '-': mon_seq = rev_comp(mon_seq)

            ext_left = ''
            ext_right = ''
            if mi == 0 and n_before_orig > 0: ext_left = flank_before
            if mi == len(display_hits)-1 and n_after_orig > 0: ext_right = flank_after

            ext_left_len = len(ext_left)
            mon_seq_len = len(mon_seq)
            mon_end_in_ext = ext_left_len + mon_seq_len
            ext_seq = ext_left + mon_seq + ext_right
            aln = pairwise_align(consensus, ext_seq)

            proj = ['-'] * cons_len
            q_pos = aln['qStart']
            for k in range(len(aln['qAln'])):
                if aln['qAln'][k] != '-':
                    if q_pos < cons_len: proj[q_pos] = aln['sAln'][k]
                    q_pos += 1

            flank_left_str = ''; mon_over_left = ''; mon_over_right = ''; flank_right_str = ''
            if aln['sStart'] > 0:
                if ext_left_len > 0:
                    flank_end = min(aln['sStart'], ext_left_len)
                    flank_left_str = ext_seq[:flank_end]
                    if aln['sStart'] > ext_left_len:
                        mon_over_left = ext_seq[ext_left_len:aln['sStart']]
                else:
                    mon_over_left = ext_seq[:aln['sStart']]
            if aln['sEnd'] < len(ext_seq):
                if ext_right:
                    if aln['sEnd'] < mon_end_in_ext:
                        mon_over_right = ext_seq[aln['sEnd']:mon_end_in_ext]
                    flank_right_str = ext_seq[mon_end_in_ext:]
                else:
                    mon_over_right = ext_seq[aln['sEnd']:]

            # Fill leading
            if mon_over_left and aln['qStart'] > 0:
                fill = min(len(mon_over_left), aln['qStart'])
                for k in range(fill):
                    proj[aln['qStart']-fill+k] = mon_over_left[len(mon_over_left)-fill+k]
                mon_over_left = mon_over_left[:len(mon_over_left)-fill]
            # Fill trailing
            if mon_over_right and aln['qEnd'] < cons_len:
                fill = min(len(mon_over_right), cons_len - aln['qEnd'])
                for k in range(fill):
                    proj[aln['qEnd']+k] = mon_over_right[k]
                mon_over_right = mon_over_right[fill:]

            # Spacer
            spacer_seq = ''; spacer_len = 0
            if mi > 0:
                prev_m = display_hits[mi-1]
                gap_start_g = min(prev_m['end'], m['end']) + 1
                gap_end_g = max(prev_m['start'], m['start']) - 1
                if gap_end_g >= gap_start_g:
                    spacer_len = gap_end_g - gap_start_g + 1
                    spacer_seq = raw_seq[gap_start_g:gap_end_g+1]
                    if strand == '-': spacer_seq = rev_comp(spacer_seq)

            g_start = m['start'] + g_off
            g_end = m['end'] + g_off
            strand_ch = '+' if strand == '+' else '-'
            label = f'#{global_idx} {g_start}-{g_end}({strand_ch})'

            proj_rows.append({
                'name': label,
                'bases': proj, 'identity': aln['identity'],
                'flankLeftStr': flank_left_str, 'monOverLeft': mon_over_left,
                'monOverRight': mon_over_right, 'flankRightStr': flank_right_str,
                'isFirst': mi == 0, 'isLast': mi == len(display_hits)-1,
                'clsat': m.get('clsat', dominant_clsat),
                'spacerSeq': spacer_seq, 'spacerLen': spacer_len,
                'gStart': g_start, 'gEnd': g_end,
                'hitSpan': m['end'] - m['start'] + 1,
                'isPartial': False,
            })

        # REFINE STEP
        for i in range(len(proj_rows)-1):
            ov = proj_rows[i]['monOverRight']
            if not ov: continue
            if proj_rows[i+1]['spacerLen'] > 0: continue
            lead_gaps = 0
            for k in range(cons_len):
                if proj_rows[i+1]['bases'][k] == '-': lead_gaps += 1
                else: break
            to_move = min(len(ov), lead_gaps)
            if to_move > 0:
                for k in range(to_move):
                    proj_rows[i+1]['bases'][k] = ov[k]
                proj_rows[i]['monOverRight'] = ov[to_move:]

        # ABSORB SHORT SPACERS
        for i in range(len(proj_rows)):
            if proj_rows[i]['spacerLen'] == 0 or proj_rows[i]['spacerLen'] > 10: continue
            sp = proj_rows[i]['spacerSeq'].upper()
            if i > 0:
                trail_gaps = 0
                for k in range(cons_len-1, -1, -1):
                    if proj_rows[i-1]['bases'][k] == '-': trail_gaps += 1
                    else: break
                fill_trail = min(len(sp), trail_gaps)
                if fill_trail > 0:
                    for k in range(fill_trail):
                        proj_rows[i-1]['bases'][cons_len-trail_gaps+k] = sp[k]
                    sp = sp[fill_trail:]
            if sp:
                lead_gaps = 0
                for k in range(cons_len):
                    if proj_rows[i]['bases'][k] == '-': lead_gaps += 1
                    else: break
                fill_lead = min(len(sp), lead_gaps)
                if fill_lead > 0:
                    for k in range(fill_lead):
                        proj_rows[i]['bases'][lead_gaps-fill_lead+k] = sp[len(sp)-fill_lead+k]
                    sp = sp[:len(sp)-fill_lead]
            proj_rows[i]['spacerSeq'] = sp
            proj_rows[i]['spacerLen'] = len(sp)

        # BOUNDARY PARTIAL DETECTION - reading-last
        for _iter in range(10):
            last_row = proj_rows[-1]
            if not last_row['isLast']: break
            boundary_seq = last_row['monOverRight'] + last_row['flankRightStr']
            if len(boundary_seq) <= 10: break
            b_aln = pairwise_align(consensus, boundary_seq)
            b_accept = (b_aln['alignLen']>10 and b_aln['identity']>0.8) or (b_aln['alignLen']>20 and b_aln['identity']>0.5)
            if not b_accept: break
            dominated = False
            for other_name in CONSENSUS_SEQS:
                if other_name == dominant_clsat: continue
                other_aln = pairwise_align(CONSENSUS_SEQS[other_name], boundary_seq)
                if other_aln['score'] > b_aln['score'] and other_aln['identity'] > b_aln['identity']:
                    dominated = True; break
            if dominated: break
            global_idx += 1
            part_proj = ['-'] * cons_len
            pq = b_aln['qStart']
            for k in range(len(b_aln['qAln'])):
                if b_aln['qAln'][k] != '-':
                    if pq < cons_len: part_proj[pq] = b_aln['sAln'][k]
                    pq += 1
            pre_bases = boundary_seq[:b_aln['sStart']]
            pre_used = 0
            if pre_bases and b_aln['qStart'] > 0:
                pre_used = min(len(pre_bases), b_aln['qStart'])
                for k in range(pre_used):
                    part_proj[b_aln['qStart']-pre_used+k] = pre_bases[len(pre_bases)-pre_used+k]
            pre_remain = pre_bases[:len(pre_bases)-pre_used]
            last_row['monOverRight'] = pre_remain if pre_remain else ''
            true_flank = boundary_seq[b_aln['sEnd']:]
            if true_flank and b_aln['qEnd'] < cons_len:
                fill = min(len(true_flank), cons_len - b_aln['qEnd'])
                for k in range(fill): part_proj[b_aln['qEnd']+k] = true_flank[k]
                true_flank = true_flank[fill:]
            proj_rows.append({
                'name': f'#{global_idx} partial({last_row["clsat"][-1]})',
                'bases': part_proj, 'identity': b_aln['identity'],
                'flankLeftStr': '', 'monOverLeft': '',
                'monOverRight': '', 'flankRightStr': true_flank,
                'isFirst': False, 'isLast': True,
                'clsat': last_row['clsat'], 'spacerLen':0, 'spacerSeq':'',
                'gStart': last_row['gEnd'], 'gEnd': last_row['gEnd'],
                'isPartial': True,
            })
            last_row['flankRightStr'] = ''
            last_row['isLast'] = False

        # BOUNDARY PARTIAL DETECTION - reading-first
        for _iter in range(10):
            first_row = proj_rows[0]
            if not first_row['isFirst']: break
            boundary_seq = first_row['flankLeftStr'] + first_row['monOverLeft']
            if len(boundary_seq) <= 10: break
            b_aln = pairwise_align(consensus, boundary_seq)
            b_accept = (b_aln['alignLen']>10 and b_aln['identity']>0.8) or (b_aln['alignLen']>20 and b_aln['identity']>0.5)
            if not b_accept: break
            dominated = False
            for other_name in CONSENSUS_SEQS:
                if other_name == dominant_clsat: continue
                other_aln = pairwise_align(CONSENSUS_SEQS[other_name], boundary_seq)
                if other_aln['score'] > b_aln['score'] and other_aln['identity'] > b_aln['identity']:
                    dominated = True; break
            if dominated: break
            global_idx += 1
            part_proj = ['-'] * cons_len
            pq = b_aln['qStart']
            for k in range(len(b_aln['qAln'])):
                if b_aln['qAln'][k] != '-':
                    if pq < cons_len: part_proj[pq] = b_aln['sAln'][k]
                    pq += 1
            post_bases = boundary_seq[b_aln['sEnd']:]
            post_used = 0
            if post_bases and b_aln['qEnd'] < cons_len:
                post_used = min(len(post_bases), cons_len - b_aln['qEnd'])
                for k in range(post_used): part_proj[b_aln['qEnd']+k] = post_bases[k]
            post_remain = post_bases[post_used:]
            first_row['monOverLeft'] = post_remain if post_remain else ''
            true_flank = boundary_seq[:b_aln['sStart']]
            pre_used = 0
            if true_flank and b_aln['qStart'] > 0:
                pre_used = min(len(true_flank), b_aln['qStart'])
                for k in range(pre_used):
                    part_proj[b_aln['qStart']-pre_used+k] = true_flank[len(true_flank)-pre_used+k]
                true_flank = true_flank[:len(true_flank)-pre_used]
            proj_rows.insert(0, {
                'name': f'#0 partial({first_row["clsat"][-1]})',
                'bases': part_proj, 'identity': b_aln['identity'],
                'flankLeftStr': true_flank, 'monOverLeft': '',
                'monOverRight': '', 'flankRightStr': '',
                'isFirst': True, 'isLast': False,
                'clsat': first_row['clsat'], 'spacerLen':0, 'spacerSeq':'',
                'gStart': first_row['gStart'], 'gEnd': first_row['gStart'],
                'isPartial': True,
            })
            first_row['flankLeftStr'] = ''
            first_row['monOverLeft'] = ''
            first_row['isFirst'] = False

        array_results.append({
            'arrayIndex': ai+1,
            'strand': strand,
            'dominantClsat': dominant_clsat,
            'consLen': cons_len,
            'nMonomers': len(arr['hits']),
            'rows': proj_rows,
        })

    return array_results


# ============================================================================
# VALIDATION CHECKS
# ============================================================================
def validate_results(tc, result, msa_arrays, raw_seq):
    errors = []
    warnings = []

    if not result.get('success'):
        errors.append(f"Analysis failed: {result.get('message','?')}")
        return errors, warnings

    if not result['arrays']:
        errors.append('No arrays found')
        return errors, warnings

    # CHECK 3: No double-width monomers
    for arr in result['arrays']:
        cons_len = len(CONSENSUS_SEQS.get(arr['hits'][0].get('clsat',''), '')) or 146
        for m in arr['hits']:
            span = m['end'] - m['start'] + 1
            if span > cons_len * 1.5:
                errors.append(f"Double-width monomer: {m['start']}-{m['end']} span={span}bp ({span/cons_len:.1f}x consensus {cons_len}bp)")

    # CHECK 4: No mixed CLsat types in array
    for arr in result['arrays']:
        types = list(set(h.get('clsat','?') for h in arr['hits']))
        if len(types) > 1:
            errors.append(f"Mixed CLsat types in array: {types} ({len(arr['hits'])} monomers)")

    for msa_arr in msa_arrays:
        cons_len = msa_arr['consLen']

        for row in msa_arr['rows']:
            # 5a: bases length
            if len(row['bases']) != cons_len:
                errors.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: bases length {len(row['bases'])} != {cons_len}")

            # 5b: identity
            if not row['isPartial'] and row['identity'] < 0.50:
                errors.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: very low identity {row['identity']*100:.1f}%")

            # 5c: excessive gap ratio
            lead_gaps = 0
            for k in range(cons_len):
                if row['bases'][k] == '-': lead_gaps += 1
                else: break
            trail_gaps = 0
            for k in range(cons_len-1, -1, -1):
                if row['bases'][k] == '-': trail_gaps += 1
                else: break
            if not row['isPartial'] and (lead_gaps + trail_gaps) > cons_len * 0.5:
                warnings.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: {lead_gaps} leading + {trail_gaps} trailing gaps (>50%)")

            # 5d: flank on wrong side
            if not row['isFirst'] and row['flankLeftStr']:
                errors.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: has left flank ({len(row['flankLeftStr'])}bp) but NOT first row")
            if not row['isLast'] and row['flankRightStr']:
                errors.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: has right flank ({len(row['flankRightStr'])}bp) but NOT last row")

        # CHECK 6: ordering
        non_partial = [r for r in msa_arr['rows'] if not r['isPartial']]
        if msa_arr['strand'] == '-' and len(non_partial) >= 2:
            for i in range(1, len(non_partial)):
                if non_partial[i]['gStart'] >= non_partial[i-1]['gStart']:
                    errors.append(f"[Array {msa_arr['arrayIndex']}] Minus-strand ordering error: {non_partial[i-1]['name']} ({non_partial[i-1]['gStart']}) then {non_partial[i]['name']} ({non_partial[i]['gStart']})")
        if msa_arr['strand'] == '+' and len(non_partial) >= 2:
            for i in range(1, len(non_partial)):
                if non_partial[i]['gStart'] <= non_partial[i-1]['gStart']:
                    errors.append(f"[Array {msa_arr['arrayIndex']}] Plus-strand ordering error: {non_partial[i-1]['name']} ({non_partial[i-1]['gStart']}) then {non_partial[i]['name']} ({non_partial[i]['gStart']})")

        # CHECK 8: isFirst/isLast consistency
        first_rows = [r for r in msa_arr['rows'] if r['isFirst']]
        last_rows = [r for r in msa_arr['rows'] if r['isLast']]
        if len(first_rows) != 1:
            errors.append(f"[Array {msa_arr['arrayIndex']}] Expected 1 isFirst row, got {len(first_rows)}")
        if len(last_rows) != 1:
            errors.append(f"[Array {msa_arr['arrayIndex']}] Expected 1 isLast row, got {len(last_rows)}")

        # CHECK 9: cross-type partials
        for row in msa_arr['rows']:
            if not row['isPartial']: continue
            bases = ''.join(b for b in row['bases'] if b != '-')
            if len(bases) < 20: continue
            aln_self = pairwise_align(CONSENSUS_SEQS[msa_arr['dominantClsat']], bases)
            for other_name in CONSENSUS_SEQS:
                if other_name == msa_arr['dominantClsat']: continue
                aln_other = pairwise_align(CONSENSUS_SEQS[other_name], bases)
                if aln_other['score'] > aln_self['score'] and aln_other['identity'] > aln_self['identity']:
                    errors.append(f"[Array {msa_arr['arrayIndex']}] {row['name']}: cross-type partial — matches {other_name} better ({aln_other['identity']*100:.0f}% vs {aln_self['identity']*100:.0f}%)")

        # CHECK 10: minus-strand first/last coords
        if msa_arr['strand'] == '-' and len(msa_arr['rows']) >= 2:
            first = next((r for r in msa_arr['rows'] if r['isFirst'] and not r['isPartial']), None)
            last = next((r for r in msa_arr['rows'] if r['isLast'] and not r['isPartial']), None)
            if first and last and first['gStart'] < last['gStart']:
                errors.append(f"[Array {msa_arr['arrayIndex']}] Minus-strand first/last flipped: first.gStart={first['gStart']} < last.gStart={last['gStart']}")

    # CHECK 11: overlapping monomers
    for arr in result['arrays']:
        s = sorted(arr['hits'], key=lambda h: h['start'])
        for i in range(1, len(s)):
            ov = s[i-1]['end'] - s[i]['start'] + 1
            if ov > 0:
                shorter = min(s[i-1]['length'], s[i]['length'])
                if ov > shorter * 0.5:
                    errors.append(f"Monomers overlap >50%: {s[i-1]['start']}-{s[i-1]['end']} and {s[i]['start']}-{s[i]['end']} (overlap={ov}bp)")

    return errors, warnings


# ============================================================================
# NCBI FETCH
# ============================================================================
def fetch_ncbi(accession, start, end):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text&seq_start={start}&seq_stop={end}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = resp.read().decode('utf-8')
    if not data.startswith('>'):
        raise ValueError(f"Not FASTA: {data[:200]}")
    seq = re.sub(r'^>.*\n', '', data, count=1).replace('\n','').replace(' ','').upper()
    return seq


# ============================================================================
# TEST CASES
# ============================================================================
TEST_CASES = [
    {'id':'L5',   'clsat':'CLsat3','species':'uni','scf':'JAWDKO010000017.1','start':8410089,'end':8413382,'size':3294},
    {'id':'L8',   'clsat':'CLsat4','species':'uni','scf':'JAWDKO010000018.1','start':148914, 'end':155000, 'size':6087},
    {'id':'L33',  'clsat':'CLsat1','species':'uni','scf':'JAWDKO010000060.1','start':10849164,'end':10851647,'size':2484},
    {'id':'L51',  'clsat':'CLsat3','species':'arm','scf':'JAWWNF010000057.1','start':533041, 'end':535120, 'size':2080},
    {'id':'L97',  'clsat':'CLsat1','species':'arm','scf':'JAWWNF010000837.1','start':59347,  'end':65999,  'size':6653},
    {'id':'L139', 'clsat':'CLsat3','species':'val','scf':'JAWWNG010000030.1','start':357599, 'end':360577, 'size':2979},
    {'id':'L142', 'clsat':'CLsat1','species':'val','scf':'JAWWNG010000038.1','start':294908, 'end':301407, 'size':6500},
    {'id':'L227', 'clsat':'CLsat1','species':'mix','scf':'JAWWNH010000002.1','start':100181, 'end':102255, 'size':2075},
    {'id':'L230', 'clsat':'CLsat2','species':'mix','scf':'JAWWNH010000040.1','start':2976815,'end':2979088,'size':2274},
    {'id':'L356', 'clsat':'CLsat4','species':'nai','scf':'JAWWNI010000088.1','start':18099,  'end':27610,  'size':9512},
    {'id':'L360', 'clsat':'CLsat3','species':'nai','scf':'JAWWNI010000158.1','start':1320321,'end':1321775,'size':1455},
    {'id':'L396', 'clsat':'CLsat1','species':'nai','scf':'JAWWNI010000715.1','start':283003, 'end':284557, 'size':1555},
]


# ============================================================================
# MAIN
# ============================================================================
def main():
    print('=' * 60)
    print('CLsat-viewer MSA Test Suite (Python)')
    print('=' * 60)
    print(f'Testing {len(TEST_CASES)} loci across 5 species, 4 CLsat types\n')

    total_errors = 0
    total_warnings = 0
    passed = 0
    failed = 0

    for ti, tc in enumerate(TEST_CASES):
        label = f"[{ti+1}/{len(TEST_CASES)}] {tc['id']} {tc['clsat']} {tc['species']} {tc['scf']}:{tc['start']}-{tc['end']} ({tc['size']}bp)"
        print(f"\n{label}")

        try:
            # Fetch
            print('  Fetching FASTA...', end='', flush=True)
            raw_seq = fetch_ncbi(tc['scf'], tc['start'], tc['end'])
            print(f' {len(raw_seq)}bp')

            if len(raw_seq) < 100:
                print('  SKIP: sequence too short')
                continue

            # Analyze
            print('  Analyzing...', end='', flush=True)
            t0 = time.time()
            result = analyze_monomers_sw(raw_seq, tc['clsat'])
            dt = time.time() - t0
            mc = result.get('monomer_count', 0)
            ac = result.get('array_count', 0)
            print(f' {dt:.1f}s, {mc} monomers, {ac} arrays')

            if not result.get('success'):
                print(f"  SKIP: {result.get('message','?')}")
                continue

            # Simulate MSA rendering
            print('  Rendering MSA...', end='', flush=True)
            msa_arrays = simulate_render_msa(result, raw_seq, tc['start'])
            print(f' {len(msa_arrays)} MSA blocks')

            for ma in msa_arrays:
                print(f"    Array {ma['arrayIndex']}: {ma['nMonomers']} monomers ({ma['strand']} strand, {ma['dominantClsat']}) -> {len(ma['rows'])} rows")

            # Validate
            errs, warns = validate_results(tc, result, msa_arrays, raw_seq)

            if not errs:
                print(f'  PASS ({len(warns)} warnings)')
                passed += 1
            else:
                print(f'  FAIL -- {len(errs)} errors:')
                for e in errs: print(f'    ERROR: {e}')
                failed += 1
            for w in warns:
                print(f'    WARN: {w}')
            total_errors += len(errs)
            total_warnings += len(warns)

        except Exception as ex:
            import traceback
            print(f'  EXCEPTION: {ex}')
            traceback.print_exc()
            failed += 1
            total_errors += 1

        # Rate limit
        if ti < len(TEST_CASES) - 1:
            time.sleep(0.5)

    print('\n' + '=' * 60)
    print(f'RESULTS: {passed} passed, {failed} failed')
    print(f'         {total_errors} errors, {total_warnings} warnings')
    print('=' * 60)
    sys.exit(1 if failed > 0 else 0)

if __name__ == '__main__':
    main()
