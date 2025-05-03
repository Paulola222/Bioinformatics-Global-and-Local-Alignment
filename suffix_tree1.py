
import sys, time, os, platform

INTERNAL_NODE_COUNTER = 0



class Node:
    def __init__(self):
        self.id = 0
        self.suff_link = None
        self.parent = None
        self.children = {}  
        self.depth = 0
        self.edge_label = (0, 0)

class LongestRepeat:
    def __init__(self):
        self.length = 0
        self.positions = []
        self.count = 0



def create_node():
    return Node()

def generate_id(is_leaf, suff_index, str_len):
    global INTERNAL_NODE_COUNTER
    if is_leaf:
        return suff_index
    else:
        node_id = str_len + INTERNAL_NODE_COUNTER
        INTERNAL_NODE_COUNTER += 1
        return node_id

def get_char_child(c, alphabet):
    return c if c in alphabet else None



def find_path(root, sequence_string, suff_index, start_pos, alphabet):
    v = root
    last_internal = root
    curr_pos = start_pos
    str_len = len(sequence_string)
    leaf_inserted = False
    while not leaf_inserted and curr_pos < str_len:
        branch_char = sequence_string[curr_pos]
        child = v.children.get(branch_char)
        if child is None:
            new_leaf = create_node()
            new_leaf.id = generate_id(True, suff_index, str_len)
            new_leaf.edge_label = (curr_pos, str_len - 1)
            new_leaf.parent = v
            new_leaf.depth = v.depth + (new_leaf.edge_label[1] - new_leaf.edge_label[0] + 1)
            v.children[branch_char] = new_leaf
            leaf_inserted = True
        else:
            edge_start, edge_end = child.edge_label
            edge_pos = edge_start
            while edge_pos <= edge_end and curr_pos < str_len and sequence_string[edge_pos] == sequence_string[curr_pos]:
                edge_pos += 1
                curr_pos += 1
            if edge_pos > edge_end:
                v = child
            else:
                new_internal = create_node()
                new_internal.id = generate_id(False, suff_index, str_len)
                new_internal.edge_label = (edge_start, edge_pos - 1)
                new_internal.parent = v
                new_internal.depth = v.depth + (edge_pos - edge_start)
                child.edge_label = (edge_pos, child.edge_label[1])
                child.parent = new_internal
                new_internal.children[sequence_string[edge_pos]] = child
                v.children[sequence_string[edge_start]] = new_internal
                new_leaf = create_node()
                new_leaf.id = generate_id(True, suff_index, str_len)
                new_leaf.edge_label = (curr_pos, str_len - 1)
                new_leaf.parent = new_internal
                new_leaf.depth = new_internal.depth + (str_len - curr_pos)
                new_internal.children[sequence_string[curr_pos]] = new_leaf
                leaf_inserted = True
    return last_internal

def node_hops(v_prime, sequence_string, suff_index, alphabet, beta_len, beta_start):
    if v_prime is None:
        sys.exit("Error: v_prime is None in node_hops")
    v = v_prime
    str_len = len(sequence_string)
    beta_counter = 0
    str_pos = beta_start
    if beta_start < 0 or beta_start >= str_len:
        sys.exit(f"Error: Invalid beta_start position {beta_start}")
    while beta_counter < beta_len and str_pos < str_len:
        current_char = sequence_string[str_pos]
        next_node = v.children.get(current_char)
        if next_node is None:
            sys.exit(f"Error in node_hops: No child for '{current_char}' at pos {str_pos}")
        edge_start, edge_end = next_node.edge_label
        edge_len = edge_end - edge_start + 1
        remaining_beta = beta_len - beta_counter
        if edge_len > remaining_beta:
            new_internal = create_node()
            new_internal.id = generate_id(False, suff_index, str_len)
            new_internal.edge_label = (next_node.edge_label[0], next_node.edge_label[0] + remaining_beta - 1)
            new_internal.parent = v
            new_internal.depth = v.depth + remaining_beta
            next_node.edge_label = (next_node.edge_label[0] + remaining_beta, next_node.edge_label[1])
            next_node.parent = new_internal
            new_internal.children[sequence_string[next_node.edge_label[0]]] = next_node
            v.children[current_char] = new_internal
            v = new_internal
            break
        else:
            beta_counter += edge_len
            str_pos += edge_len
            v = next_node
    return v

def suff_link_known(u, sequence_string, suff_index, alphabet):
    v = u.suff_link
    k = v.depth
    if suff_index + k <= len(sequence_string):
        return find_path(v, sequence_string, suff_index, suff_index + k, alphabet)
    return v

def suff_link_unknown_internal(u, sequence_string, suff_index, alphabet):
    u_prime = u.parent
    v_prime = u_prime.suff_link
    u_start_edge = u.edge_label[0]
    beta_len = u.edge_label[1] - u.edge_label[0] + 1
    v = node_hops(v_prime, sequence_string, suff_index, alphabet, beta_len, u_start_edge)
    u.suff_link = v
    alpha = v.depth
    return find_path(v, sequence_string, suff_index, suff_index + alpha, alphabet)

def suff_link_unknown_root(u, sequence_string, suff_index, alphabet):
    str_len = len(sequence_string)
    u_prime = u.parent
    beta_len = u.depth - 1
    beta_start = u.edge_label[0] + 1
    v = node_hops(u_prime, sequence_string, suff_index, alphabet, beta_len, beta_start)
    u.suff_link = v
    alpha = v.depth
    return find_path(v, sequence_string, suff_index, suff_index + alpha, alphabet)

def build_suffix_tree(sequence_string, alphabet, is_naive=False):
    str_len = len(sequence_string)
    root = create_node()
    root.suff_link = root
    root.parent = root
    root.id = generate_id(False, 0, str_len)
    if is_naive:
        for suff_ind in range(str_len):
            find_path(root, sequence_string, suff_ind, suff_ind, alphabet)
    else:
        last_internal = None
        for suff_ind in range(str_len):
            u = last_internal if last_internal is not None else root
            if u.suff_link is not None:
                last_internal = suff_link_known(u, sequence_string, suff_ind, alphabet)
            elif u != root:
                last_internal = suff_link_unknown_internal(u, sequence_string, suff_ind, alphabet)
            else:
                last_internal = suff_link_unknown_root(u, sequence_string, suff_ind, alphabet)
    return root



def is_leaf(node, str_len):
    return node.id < str_len

def is_root(node):
    return node.suff_link == node

def compute_tree_stats(root, seq_str, alphabet):
    stats = {"total_nodes": 0, "internal_nodes": 0, "leaves": 0, "internal_depth_sum": 0, "deepest_internal": 0}
    str_len = len(seq_str)
    def traverse(n):
        stats["total_nodes"] += 1
        if is_leaf(n, str_len):
            stats["leaves"] += 1
        else:
            stats["internal_nodes"] += 1
            stats["internal_depth_sum"] += n.depth
            if n.depth > stats["deepest_internal"]:
                stats["deepest_internal"] = n.depth
        for letter in alphabet:
            child = n.children.get(letter)
            if child:
                traverse(child)
    traverse(root)
    avg_depth = stats["internal_depth_sum"] / stats["internal_nodes"] if stats["internal_nodes"] > 0 else 0
    return stats, avg_depth

def write_dfs_to_file(root, seq_str, alphabet, filename):
    with open(filename, "w") as f:
        def dfs(n):
            if n.parent == n:
                line = f"[DFS Root id={n.id}, depth={n.depth}]"
            else:
                edge = seq_str[n.edge_label[0]: n.edge_label[1] + 1]
                line = f"[DFS Node id={n.id}, depth={n.depth}, edge='{edge}']"
            f.write(line + "\n")
            for letter in alphabet:
                child = n.children.get(letter)
                if child:
                    dfs(child)
        dfs(root)
    print(f"DFS enumeration written to: {filename}")

def write_postorder_to_file(root, seq_str, alphabet, filename):
    with open(filename, "w") as f:
        def postorder(n):
            for letter in alphabet:
                child = n.children.get(letter)
                if child:
                    postorder(child)
            if n.parent == n:
                line = f"[PostOrder Root id={n.id}, depth={n.depth}]"
            else:
                edge = seq_str[n.edge_label[0]: n.edge_label[1] + 1]
                line = f"[PostOrder Node id={n.id}, depth={n.depth}, edge='{edge}']"
            f.write(line + "\n")
        postorder(root)
    print(f"Post-order enumeration written to: {filename}")

def print_suffix_tree(node, seq_str, alphabet, str_len, depth=0):
    indent = "  " * depth
    if is_root(node):
        print(f"{indent}[Root id={node.id}]")
    else:
        edge = seq_str[node.edge_label[0]: node.edge_label[1] + 1]
        if is_leaf(node, str_len):
            print(f"{indent}[Leaf id={node.id}, suffix={node.id}, edge='{edge}']")
        else:
            print(f"{indent}[Internal id={node.id}, edge='{edge}']")
    if node.suff_link is not None:
        print(f"{indent} --> [id={node.suff_link.id}]")
    for letter in alphabet:
        child = node.children.get(letter)
        if child:
            print_suffix_tree(child, seq_str, alphabet, str_len, depth + 1)

def print_tree(root, seq_str, alphabet):
    str_len = len(seq_str)
    print(f"\nSuffix Tree for: '{seq_str}' (length={str_len})")
    print(f"Alphabet: '{alphabet}'")
    print("Tree structure:")
    print_suffix_tree(root, seq_str, alphabet, str_len)
    print()



def compute_bwt_index(root, sequence_file, seq_str, alphabet):
    n = len(seq_str)
    BWT = [''] * n
    bwt_index = 0
    stack = [root]
    while stack:
        curr = stack.pop()
        if is_leaf(curr, n):
            suffix_id = curr.id
            bwt_pos = n - 1 if suffix_id == 0 else suffix_id - 1
            BWT[bwt_index] = seq_str[bwt_pos]
            bwt_index += 1
        for letter in reversed(alphabet):
            child = curr.children.get(letter)
            if child:
                stack.append(child)
    base_name = os.path.splitext(sequence_file)[0]
    output_filename = f"{base_name}_bwt.txt"
    with open(output_filename, "w") as f:
        for ch in BWT:
            f.write(f"{ch}\n")
    print(f"BWT index written to: {output_filename}")



def report_space_usage(root, seq_str, alphabet):
    input_bytes = len(seq_str)
    import sys
    node_size = sys.getsizeof(Node()) + (len(alphabet) * 8)
    estimated_memory = (input_bytes * 2) * node_size
    space_constant = estimated_memory / input_bytes
    print("Space Usage Estimate:")
    print(f"Input size: {input_bytes} bytes")
    print(f"Estimated tree memory: {estimated_memory} bytes (~{estimated_memory/(1024*1024):.2f} MB)")
    print(f"Space constant: ~{space_constant:.1f} bytes per input byte")



def collect_leaf_positions(node, seq_str, alphabet, result):
    if node is None:
        return
    if is_leaf(node, len(seq_str)):
        result.positions.append(node.id)
        result.count += 1
    else:
        for letter in alphabet:
            child = node.children.get(letter)
            if child:
                collect_leaf_positions(child, seq_str, alphabet, result)

def find_longest_repeat(node, seq_str, alphabet, result, current_depth):
    if node is None:
        return
    edge_length = (node.edge_label[1] - node.edge_label[0] + 1) if node.parent != node else 0
    current_depth += edge_length
    child_count = sum(1 for letter in alphabet if letter in node.children)
    if child_count >= 2 and current_depth > result.length:
        result.length = current_depth
        result.positions = []
        result.count = 0
        collect_leaf_positions(node, seq_str, alphabet, result)
    for letter in alphabet:
        child = node.children.get(letter)
        if child:
            find_longest_repeat(child, seq_str, alphabet, result, current_depth)

def find_repeats(root, seq_str, alphabet):
    result = LongestRepeat()
    find_longest_repeat(root, seq_str, alphabet, result, 0)
    return result

def print_repeats(repeat, seq_str):
    if repeat.length == 0:
        print("No repeats found.")
        return
    repeat_str = seq_str[repeat.positions[0]: repeat.positions[0] + repeat.length]
    print(f"Longest exact repeat: '{repeat_str}' (length {repeat.length})")
    print("Positions:", ", ".join(str(pos) for pos in repeat.positions))



def report_system_configuration():
    print("System Configuration:")
    print(f"Processor: {platform.processor()}")
    print(f"Machine: {platform.machine()}")
    print(f"System: {platform.system()} {platform.version()}\n")



def read_string_sequence(filename, num_seq=1):
    sequences = []
    try:
        with open(filename, "r") as f:
            seq_name, seq_data = None, ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq_name is not None:
                        sequences.append({"name": seq_name, "sequence": seq_data})
                        seq_data = ""
                    seq_name = line[1:]
                else:
                    seq_data += line
            if seq_name is not None:
                sequences.append({"name": seq_name, "sequence": seq_data})
    except IOError as e:
        sys.exit(f"Error opening file: {e}")
    if not sequences:
        sys.exit("No sequence found in the file.")
    sequences[0]["sequence"] += "$"
    return sequences

def read_alphabet(filename):
    try:
        with open(filename, "r") as f:
            content = f.read().strip()
    except IOError as e:
        sys.exit(f"Error opening alphabet file: {e}")
    seen = set()
    alphabet = "$"
    for ch in content:
        if ch.isalnum() and ch not in seen:
            seen.add(ch)
            alphabet += ch
    return alphabet



def main():
    if len(sys.argv) < 3:
        print("Usage: python suffix_tree.py <sequence_file> <alphabet_file>")
        sys.exit(1)
    report_system_configuration()
    sequence_file = sys.argv[1]
    alphabet_file = sys.argv[2]
    print(f"Sequence File: {sequence_file}")
    sequences = read_string_sequence(sequence_file)
    seq_name = sequences[0]["name"]
    seq_str = sequences[0]["sequence"]
    print(f"Sequence Name: {seq_name}")
    print(f"Alphabet File: {alphabet_file}")
    alphabet = read_alphabet(alphabet_file)
    print(f"Alphabet: {alphabet}")
    print("**************************************************")
    start_time = time.time()
    root = build_suffix_tree(seq_str, alphabet, is_naive=False)
    end_time = time.time()
    print(f"Suffix Tree Construction Time: {(end_time - start_time):.6f} seconds")
    print("**************************************************")
    report_space_usage(root, seq_str, alphabet)
    print("**************************************************")
    stats, avg_depth = compute_tree_stats(root, seq_str, alphabet)
    print("Suffix Tree Statistics:")
    print(f"  Internal nodes: {stats['internal_nodes']}")
    print(f"  Leaves: {stats['leaves']}")
    print(f"  Total nodes: {stats['total_nodes']}")
    print(f"  Average string-depth of internal nodes: {avg_depth:.2f}")
    print(f"  Deepest internal node string-depth: {stats['deepest_internal']}")
    print("**************************************************")
    compute_bwt_index(root, sequence_file, seq_str, alphabet)
    print("**************************************************")
    repeats = find_repeats(root, seq_str, alphabet)
    print_repeats(repeats, seq_str)
    print("**************************************************")
    
    write_dfs_to_file(root, seq_str, alphabet, "dfs_output.txt")
    write_postorder_to_file(root, seq_str, alphabet, "postorder_output.txt")

if __name__ == "__main__":
    main()
