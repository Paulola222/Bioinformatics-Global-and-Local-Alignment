#!/usr/bin/env python3
class _State:
    __slots__ = ('length','link','next','first_pos')
    def __init__(self):
        self.length = 0
        self.link = -1
        self.next = {}
        self.first_pos = -1

def _build_automaton(s):
    sa = [_State()]
    last = 0
    for i, c in enumerate(s):
        cur = len(sa)
        sa.append(_State())
        sa[cur].length = sa[last].length + 1
        sa[cur].first_pos = i
        p = last
        while p >= 0 and c not in sa[p].next:
            sa[p].next[c] = cur
            p = sa[p].link
        if p == -1:
            sa[cur].link = 0
        else:
            q = sa[p].next[c]
            if sa[p].length + 1 == sa[q].length:
                sa[cur].link = q
            else:
                clone = len(sa)
                sa.append(_State())
                sa[clone].length = sa[p].length + 1
                sa[clone].next = sa[q].next.copy()
                sa[clone].link = sa[q].link
                sa[clone].first_pos = sa[q].first_pos
                while p >= 0 and sa[p].next[c] == q:
                    sa[p].next[c] = clone
                    p = sa[p].link
                sa[q].link = sa[cur].link = clone
        last = cur
    return sa

def find_lcs(s1, s2):
    sa = _build_automaton(s1)
    v = 0
    l = 0
    best = 0
    best_end2 = 0
    best_state = 0
    for i, c in enumerate(s2):
        if c in sa[v].next:
            v = sa[v].next[c]
            l += 1
        else:
            while v != -1 and c not in sa[v].next:
                v = sa[v].link
            if v == -1:
                v = 0
                l = 0
            else:
                l = sa[v].length + 1
                v = sa[v].next[c]
        if l > best:
            best = l
            best_end2 = i
            best_state = v
    if best == 0:
        return 0, (0, 0)
    end1 = sa[best_state].first_pos
    start1 = end1 - best + 1
    end2 = best_end2
    start2 = end2 - best + 1
    return best, (start1, start2)

def global_align(s1, s2, match=1, mismatch=-2, gap_open=-5, gap_extend=-1):
    m, n = len(s1), len(s2)
    if m == 0 or n == 0:
        return 0, 0
    INF = -10**9
    prevS = [INF]*(n+1)
    prevD = [INF]*(n+1)
    prevI = [INF]*(n+1)
    prevS[0] = 0
    prevD[0] = prevI[0] = gap_open
    for j in range(1, n+1):
        prevI[j] = gap_open + (j-1)*gap_extend
    max_score = INF
    mi = mj = 0
    state = 'S'
    for i in range(1, m+1):
        currS = [INF]*(n+1)
        currD = [INF]*(n+1)
        currI = [INF]*(n+1)
        currD[0] = gap_open + (i-1)*gap_extend
        for j in range(1, n+1):
            sub = match if s1[i-1] == s2[j-1] else mismatch
            s_val = max(prevS[j-1], prevD[j-1], prevI[j-1]) + sub
            d_val = max(prevS[j] + gap_open, prevD[j] + gap_extend)
            i_val = max(currS[j-1] + gap_open, currI[j-1] + gap_extend)
            currS[j], currD[j], currI[j] = s_val, d_val, i_val
            for sc, st in ((s_val,'S'), (d_val,'D'), (i_val,'I')):
                if sc > max_score:
                    max_score, mi, mj, state = sc, i, j, st
        prevS, prevD, prevI = currS, currD, currI
    matches = 0
    i, j = mi, mj
    while i > 0 and j > 0:
        if state == 'S':
            if s1[i-1] == s2[j-1]:
                matches += 1
            i -= 1; j -= 1; state = 'S'
        elif state == 'D':
            i -= 1; state = 'S'
        else:
            j -= 1; state = 'S'
    return max_score, matches
