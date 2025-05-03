#!/usr/bin/env python3
import sys
import csv
from suffix_tree import find_lcs, global_align

def read_fasta(path):
    seqs = []
    with open(path) as f:
        name = None
        buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs.append((name, ''.join(buf)))
                name = line[1:]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs.append((name, ''.join(buf)))
    return seqs

def main():
    if len(sys.argv) != 3:
        print("Usage: python similarity.py all_sequences.fasta DNA_alphabet.txt")
        sys.exit(1)
    fasta, _ = sys.argv[1], sys.argv[2]
    seqs = read_fasta(fasta)
    k = len(seqs)
    names = [n for n,_ in seqs]

    lcs_mat = [[0]*k for _ in range(k)]
    sim_mat = [[0]*k for _ in range(k)]

    for i in range(k):
        Li = len(seqs[i][1])
        lcs_mat[i][i] = Li
        sim_mat[i][i] = Li
        for j in range(i+1, k):
            print(f"Comparing {names[i]} vs {names[j]}", flush=True)
            s1 = seqs[i][1]
            s2 = seqs[j][1]
            b, (x1, x2) = find_lcs(s1, s2)
            pre1 = s1[:x1][::-1]
            pre2 = s2[:x2][::-1]
            _, a = global_align(pre1, pre2)
            y1, y2 = x1 + b, x2 + b
            suf1 = s1[y1:]
            suf2 = s2[y2:]
            _, c = global_align(suf1, suf2)
            score = a + b + c
            lcs_mat[i][j] = lcs_mat[j][i] = b
            sim_mat[i][j] = sim_mat[j][i] = score

    with open('lcs_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([''] + names)
        for name,row in zip(names, lcs_mat):
            w.writerow([name] + row)

    with open('similarity_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([''] + names)
        for name,row in zip(names, sim_mat):
            w.writerow([name] + row)

    print("Wrote lcs_results.csv and similarity_results.csv")

if __name__ == '__main__':
    main()
