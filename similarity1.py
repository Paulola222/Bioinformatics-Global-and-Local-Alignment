import sys
import time
from collections import defaultdict, deque


class Node:
    def __init__(self, start=None, end=None):
        self.children = {}
        self.suffix_link = None
        self.start = start
        self.end = end
        self.index = -1

class SuffixTree:
    def __init__(self, text, alphabet):
        self.text = text
        self.alphabet = alphabet
        self.root = Node()
        self.build()

    def build(self):
        s = self.text
        n = len(s)
        root = self.root
        root.suffix_link = root
        active_node = root
        active_edge = 0
        active_length = 0
        remainder = 0
        last_new_node = None

    
        def edge_length(node):
            return (node.end or n-1) - node.start + 1

        for pos in range(n):
            remainder += 1
            last_new_node = None
            while remainder > 0:
                if active_length == 0:
                    active_edge = pos
                edge_char = s[active_edge]
                if edge_char not in active_node.children:
                    
                    leaf = Node(pos, None)
                    active_node.children[edge_char] = leaf
                    leaf.index = pos - remainder + 1
                    if last_new_node:
                        last_new_node.suffix_link = active_node
                        last_new_node = None
                else:
                    next_node = active_node.children[edge_char]
                    if active_length >= edge_length(next_node):
                        active_edge += edge_length(next_node)
                        active_length -= edge_length(next_node)
                        active_node = next_node
                        continue
                    if s[next_node.start + active_length] == s[pos]:
                        active_length += 1
                        if last_new_node:
                            last_new_node.suffix_link = active_node
                            last_new_node = None
                        break
                    
                    split_end = next_node.start + active_length - 1
                    split = Node(next_node.start, split_end)
                    active_node.children[edge_char] = split
                    leaf = Node(pos, None)
                    split.children[s[pos]] = leaf
                    leaf.index = pos - remainder + 1
                    next_node.start += active_length
                    split.children[s[next_node.start]] = next_node
                    if last_new_node:
                        last_new_node.suffix_link = split
                    last_new_node = split
                    split.suffix_link = root
                remainder -= 1
                if active_node == root and active_length > 0:
                    active_length -= 1
                    active_edge = pos - remainder + 1
                else:
                    active_node = active_node.suffix_link or root
        
        def set_leaf_ends(node, cur_end):
            if not node.children:
                if node.end is None:
                    node.end = cur_end
            else:
                for c in node.children.values():
                    set_leaf_ends(c, cur_end)
        set_leaf_ends(root, n-1)

    def longest_common_substring(self, s1, s2):
        
        best = {'length':0, 'pos':(0,0)}
        def dfs(node, depth):
            mask = 0
            if not node.children:
                idx = node.index
                mask = 1 if idx < len(s1) else 2
                return mask
            for child in node.children.values():
                m = dfs(child, depth + (child.end - child.start + 1))
                mask |= m
            if mask == 3 and depth > best['length']:
                best['length'] = depth
               
                best['pos'] = (find_leaf(node, True), find_leaf(node, False))
            return mask

        def find_leaf(node, in_s1):
            stack = [node]
            while stack:
                v = stack.pop()
                if not v.children:
                    idx = v.index
                    if in_s1 and idx < len(s1): return idx
                    if not in_s1 and idx >= len(s1): return idx - len(s1)
                else:
                    for c in v.children.values(): stack.append(c)
            return 0

        dfs(self.root, 0)
        return best['pos'][0], best['pos'][1], best['length']


MATCH = 1
MISMATCH = -1
GAP_OPEN = -2
GAP_EXTEND = -1

def global_align(s1, s2):
    m, n = len(s1), len(s2)
    
    S = [[0]*(n+1) for _ in range(m+1)]
    D = [[0]*(n+1) for _ in range(m+1)]
    I = [[0]*(n+1) for _ in range(m+1)]
    
    S[0][0] = 0
    for i in range(1, m+1):
        S[i][0] = float('-inf'); D[i][0] = GAP_OPEN + (i-1)*GAP_EXTEND; I[i][0] = float('-inf')
    for j in range(1, n+1):
        S[0][j] = float('-inf'); I[0][j] = GAP_OPEN + (j-1)*GAP_EXTEND; D[0][j] = float('-inf')
    
    max_score = float('-inf')
    max_i = max_j = 0
    for i in range(1,m+1):
        for j in range(1,n+1):
            sub = MATCH if s1[i-1]==s2[j-1] else MISMATCH
            
            S[i][j] = max(S[i-1][j-1], D[i-1][j-1], I[i-1][j-1]) + sub
            
            D[i][j] = max(S[i-1][j] + GAP_OPEN, D[i-1][j] + GAP_EXTEND)
            
            I[i][j] = max(S[i][j-1] + GAP_OPEN, I[i][j-1] + GAP_EXTEND)
            
            for score in (S[i][j], D[i][j], I[i][j]):
                if score > max_score:
                    max_score, max_i, max_j = score, i, j
    
    i,j = max_i, max_j; matches=0; state='S'
    while i>0 and j>0:
        if state=='S':
            if s1[i-1]==s2[j-1]: matches+=1
            prev = max((S[i-1][j-1], 'S'), (D[i-1][j-1], 'D'), (I[i-1][j-1], 'I'))[1]
            i-=1; j-=1; state=prev
        elif state=='D':
            prev = 'S' if S[i-1][j]+GAP_OPEN>=D[i-1][j]+GAP_EXTEND else 'D'
            i-=1; state=prev
        else:
            prev = 'S' if S[i][j-1]+GAP_OPEN>=I[i][j-1]+GAP_EXTEND else 'I'
            j-=1; state=prev
    return matches


def compute_similarity_matrix(seqs, alphabet):
    n = len(seqs)
    D = [[(0,0)]*n for _ in range(n)]
    for i in range(n):
        D[i][i] = (len(seqs[i]) - 1, len(seqs[i]) - 1)
        for j in range(i+1, n):
            concat = seqs[i] + seqs[j]
            t0 = time.time()
            st = SuffixTree(concat, alphabet)
            st_time = time.time() - t0
            x1,x2,b = st.longest_common_substring(seqs[i], seqs[j])
            
            a = global_align(seqs[i][:x1][::-1], seqs[j][:x2][::-1])
            
            y1 = x1+b; y2=x2+b
            c = global_align(seqs[i][y1:], seqs[j][y2:])
            score = a + b + c
            D[i][j] = D[j][i] = (score, b)
    return D


def load_fasta(fname):
    seqs = []
    with open(fname) as f:
        header=None; seq=[]
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if header: seqs.append(''.join(seq)+'$')
                header=line; seq=[]
            else:
                seq.append(line)
        if header: seqs.append(''.join(seq)+'$')
    return seqs

if __name__=='__main__':
    fasta_file, alpha_file = sys.argv[1], sys.argv[2]
    seqs = load_fasta(fasta_file)
    alphabet = open(alpha_file).read().strip()
    matrix = compute_similarity_matrix(seqs, alphabet)
    
    for row in matrix:
        print("\t".join(f"{cell[0]}({cell[1]})" for cell in row))
