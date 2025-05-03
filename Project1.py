import sys

def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

def create_matrix(rows, cols):
    return [[0] * cols for _ in range(rows)]

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend):
    n, m = len(seq1), len(seq2)

    
    subtitue = create_matrix(n+1, m+1)
    deletion = create_matrix(n+1, m+1)
    insertion = create_matrix(n+1, m+1)
    traceback = create_matrix(n+1, m+1)

   
    subtitue[0][0] = 0
    for i in range(1, n+1):
        deletion[i][0] = gap_open + (i-1) * gap_extend
        subtitue[i][0] = float('-inf')
        insertion[i][0] = float('-inf')
    for j in range(1, m+1):
        insertion[0][j] = gap_open + (j-1) * gap_extend
        subtitue[0][j] = float('-inf')
        deletion[0][j] = float('-inf')

   
    for i in range(1, n+1):
        for j in range(1, m+1):
            
            prev_max = max(subtitue[i-1][j-1], deletion[i-1][j-1], insertion[i-1][j-1])
            match = prev_max + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            deletion[i][j] = max(subtitue[i-1][j] + gap_open + gap_extend,
                                 deletion[i-1][j] + gap_extend)
            insertion[i][j] = max(subtitue[i][j-1] + gap_open + gap_extend,
                                  insertion[i][j-1] + gap_extend)
            subtitue[i][j] = max(match, deletion[i][j], insertion[i][j])
            
            
            if subtitue[i][j] == match:
                traceback[i][j] = 'diag'
            elif subtitue[i][j] == deletion[i][j]:
                traceback[i][j] = 'up'
            else:
                traceback[i][j] = 'left'

    
    score = subtitue[n][m]
    
    
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 'diag':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or traceback[i][j] == 'up'):
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return score, align1, align2

def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend):
    n, m = len(seq1), len(seq2)

    subtitue = create_matrix(n+1, m+1)
    deletion = create_matrix(n+1, m+1)
    insertion = create_matrix(n+1, m+1)
    traceback = create_matrix(n+1, m+1)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n+1):
        for j in range(1, m+1):
            prev_max = max(subtitue[i-1][j-1], deletion[i-1][j-1], insertion[i-1][j-1])
            match = prev_max + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            deletion[i][j] = max(subtitue[i-1][j] + gap_open + gap_extend, deletion[i-1][j] + gap_extend)
            insertion[i][j] = max(subtitue[i][j-1] + gap_open + gap_extend, insertion[i][j-1] + gap_extend)
            subtitue[i][j] = max(0, match, deletion[i][j], insertion[i][j])

            if subtitue[i][j] == match:
                traceback[i][j] = 'diag'
            elif subtitue[i][j] == deletion[i][j]:
                traceback[i][j] = 'up'
            elif subtitue[i][j] == insertion[i][j]:
                traceback[i][j] = 'left'

            if subtitue[i][j] > max_score:
                max_score = subtitue[i][j]
                max_pos = (i, j)

    align1, align2 = '', ''
    i, j = max_pos
    while i > 0 and j > 0 and subtitue[i][j] > 0:
        if traceback[i][j] == 'diag':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == 'up':
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return max_score, align1, align2

def print_alignment(align1, align2, seq1_name, seq2_name, line_length=60):
    for i in range(0, len(align1), line_length):
        align1_segment = align1[i:i+line_length]
        align2_segment = align2[i:i+line_length]
        match_line = ''.join('|' if a != '-' and b != '-' and a == b else ' ' for a, b in zip(align1_segment, align2_segment))
        print(f"{seq1_name} {i+1:4} {align1_segment} {i+len(align1_segment):4}")
        print(f"          {match_line}")
        print(f"{seq2_name} {i+1:4} {align2_segment} {i+len(align2_segment):4}")
        print()

def calculate_statistics(align1, align2):
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-' and b != '-')
    mismatches = sum(1 for a, b in zip(align1, align2) if a != b and a != '-' and b != '-')
    gaps = sum(1 for a, b in zip(align1, align2) if a == '-' or b == '-')
    
    gap_openings = 0
    if align1[0] == '-' or align2[0] == '-':
        gap_openings += 1
    for i in range(1, len(align1)):
        if (align1[i] == '-' and align1[i-1] != '-') or (align2[i] == '-' and align2[i-1] != '-'):
            gap_openings += 1

    gap_extensions = gaps

    identity = matches / len(align1) * 100
    return matches, mismatches, gap_openings, gap_extensions, identity

def main():
    if len(sys.argv) != 6:
        print("Usage: python alignment.py <seq_file> <match_score> <mismatch_penalty> <gap_open> <gap_extend>")
        sys.exit(1)

    seq_file = sys.argv[1]
    match_score = int(sys.argv[2])
    mismatch_penalty = int(sys.argv[3])
    gap_open = int(sys.argv[4])
    gap_extend = int(sys.argv[5])

    sequences = read_fasta(seq_file)
    if len(sequences) < 2:
        print("Error: FASTA file should contain at least two sequences.")
        sys.exit(1)
    
    seq1, seq2 = sequences[0], sequences[1]

    print("Global Alignment (Needleman-Wunsch):")
    global_score, global_align1, global_align2 = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend)
    print_alignment(global_align1, global_align2, "Seq1", "Seq2")
    matches, mismatches, gap_openings, gap_extensions, identity = calculate_statistics(global_align1, global_align2)
    print(f"Global optimal score: {global_score}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Opening gaps: {gap_openings}")
    print(f"Gap extensions: {gap_extensions}")
    print(f"Identity percentage: {identity:.2f}%")
    print()

    print("Local Alignment (Smith-Waterman):")
    local_score, local_align1, local_align2 = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend)
    print_alignment(local_align1, local_align2, "Seq1", "Seq2")
    matches, mismatches, gap_openings, gap_extensions, identity = calculate_statistics(local_align1, local_align2)
    print(f"Local optimal score: {local_score}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Opening gaps: {gap_openings}")
    print(f"Gap extensions: {gap_extensions}")
    print(f"Identity percentage: {identity:.2f}%")

if __name__ == "__main__":
    main()
