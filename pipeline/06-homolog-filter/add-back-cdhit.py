def parse_fasta_with_prefix(file_path):
    """Parses a FASTA file and returns mappings for sequences, headers, and prefixes."""
    seqid_to_seq = {}
    seqid_to_header = {}
    seqid_to_prefix = {}
    with open(file_path, 'r') as f:
        seq_id = ''
        header = ''
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    seqid_to_seq[seq_id] = sequence
                    seqid_to_header[seq_id] = header
                    prefix = extract_prefix(seq_id)
                    seqid_to_prefix[seq_id] = prefix
                header = line
                sequence = ''
                seq_id = line[1:].split()[0]
            else:
                sequence += line
        # Add last sequence
        if seq_id:
            seqid_to_seq[seq_id] = sequence
            seqid_to_header[seq_id] = header
            prefix = extract_prefix(seq_id)
            seqid_to_prefix[seq_id] = prefix
    return seqid_to_seq, seqid_to_header, seqid_to_prefix

def extract_prefix(seq_id):
    """Extracts the prefix (e.g., EP00314) from the sequence ID."""
    return seq_id.split('_', 1)[0]

def parse_cluster_file(cluster_file):
    """
    Parses the CD-HIT cluster file and returns a list of clusters.
    Each cluster is a list of tuples: (seq_id, prefix, length, is_star).

    We extract:
    - length from the segment before 'aa,'
    - seq_id from after the '>' and before '...'
    - prefix from seq_id
    - star presence from '*'
    """
    clusters = []
    current_cluster = []
    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                if current_cluster:
                    clusters.append(current_cluster)
                current_cluster = []
            else:
                # Example line:
                # 0   378aa, >EP00742_Colponemidia_sp_Colp-10_P020711... at 100.00%
                # length = 378 (before 'aa,')
                length_str = line.split('aa,')[0].split()[-1]  # last token before 'aa,' is length
                length_val = int(length_str)

                # Extract seq_id
                after_aa = line.split('aa,', 1)[1].strip()
                seq_part = after_aa.split('>')[1]
                seq_id = seq_part.split('...')[0].strip()

                prefix = extract_prefix(seq_id)
                is_star = '*' in line
                current_cluster.append((seq_id, prefix, length_val, is_star))
        if current_cluster:
            clusters.append(current_cluster)
    return clusters

def identify_minority_sequences(clusters):
    """
    WHY? If there are different EP ids are in the cluster, the minority might me removed but by this step it is readded as a representative.
    
    Identifies minority sequences with the following rules:

    1. Determine majority/minority by count and star logic:
       - If one prefix has the highest count, it's majority; others minority.
       - If tie in count:
         * If exactly one star prefix in tie, that star prefix is majority; others minority.
         * If multiple star prefixes tie, non-star tied prefixes are minority.
         * If no star in tie, pick one prefix as majority arbitrarily, others minority.

    2. Star Override for fewer prefixes:
       - If a prefix is fewer-copy and has a star, it is NOT minority (remove it from minority).

    3. Longest Sequence Logic:
       - In a tie scenario where a single star prefix is chosen as majority and multiple prefixes become minority,
         instead of adding all minority sequences, add only the single longest sequence among all those minority prefixes combined.

    We'll implement step 3 after determining minority.

    Returns a set of minority sequence IDs.
    """

    minority_seq_ids = set()

    for cluster in clusters:
        prefix_count = {}
        prefix_has_star = set()
        prefix_to_seqs = {}  # prefix -> list of (seq_id, length, is_star)
        for seq_id, prefix, length, is_star in cluster:
            prefix_count[prefix] = prefix_count.get(prefix, 0) + 1
            if prefix not in prefix_to_seqs:
                prefix_to_seqs[prefix] = []
            prefix_to_seqs[prefix].append((seq_id, length, is_star))
            if is_star:
                prefix_has_star.add(prefix)

        if len(prefix_count) <= 1:
            # Only one prefix, no minority
            continue

        print("Cluster with multiple prefixes: {{prefix_count}}")
        max_count = max(prefix_count.values())
        max_prefixes = [p for p, c in prefix_count.items() if c == max_count]

        # Determine initial minority prefixes
        cluster_minority_prefixes = set()
        tie_scenario = len(max_prefixes) > 1
        single_star_majority = False

        if len(max_prefixes) == 1:
            # One clear majority
            majority_prefix = max_prefixes[0]
            for p, c in prefix_count.items():
                if p != majority_prefix:
                    cluster_minority_prefixes.add(p)
        else:
            # Tie scenario
            star_in_tie = [p for p in max_prefixes if p in prefix_has_star]
            if len(star_in_tie) == 1:
                # One star prefix wins
                star_majority = star_in_tie[0]
                single_star_majority = True
                # Other tied max are minority
                for p in max_prefixes:
                    if p != star_majority:
                        cluster_minority_prefixes.add(p)
                # Also fewer count prefixes minority
                for p, c in prefix_count.items():
                    if c < max_count:
                        cluster_minority_prefixes.add(p)
            elif len(star_in_tie) > 1:
                # Multiple star prefixes tie
                # Non-star tied prefixes are minority
                non_star_in_tie = [p for p in max_prefixes if p not in prefix_has_star]
                for p in non_star_in_tie:
                    cluster_minority_prefixes.add(p)
                # Fewer count also minority
                for p, c in prefix_count.items():
                    if c < max_count:
                        cluster_minority_prefixes.add(p)
            else:
                # Tie, no star prefix
                chosen_majority = max_prefixes[0]
                for p in max_prefixes[1:]:
                    cluster_minority_prefixes.add(p)
                # Fewer count also minority
                for p, c in prefix_count.items():
                    if c < max_count:
                        cluster_minority_prefixes.add(p)

        # Star override for fewer prefixes
        # If a fewer prefix has a star, remove it from minority
        for p, c in prefix_count.items():
            if p in cluster_minority_prefixes and c < max_count and p in prefix_has_star:
                cluster_minority_prefixes.remove(p)

        # Now we have the final minority prefixes based on initial logic
        # If in a tie scenario with a single star majority and multiple minority prefixes:
        # Instead of adding all minority sequences, add only the single longest sequence.
        # Check conditions:
        # - tie scenario
        # - single_star_majority == True
        # - more than one minority prefix or more than one sequence in minority?

        minority_sequences = []
        for p in cluster_minority_prefixes:
            # Add all sequences from these minority prefixes
            for (seq_id, length, is_star) in prefix_to_seqs[p]:
                minority_sequences.append((seq_id, length))

        if tie_scenario and single_star_majority and len(minority_sequences) > 1:
            # Pick only the longest sequence
            longest_seq = max(minority_sequences, key=lambda x: x[1])
            minority_seq_ids.add(longest_seq[0])
        else:
            # Add all minority sequences normally
            for seq_id, length in minority_sequences:
                minority_seq_ids.add(seq_id)

    print("Identified minority sequence IDs: {{minority_seq_ids}}")
    return minority_seq_ids

def write_sequences_to_file(seq_ids, seqid_to_seq, seqid_to_header, output_file):
    """Writes sequences corresponding to seq_ids to the output file."""
    if not seq_ids:
        print("No sequences to re-add.")
        return
    with open(output_file, 'a') as f:
        for seq_id in seq_ids:
            header = seqid_to_header.get(seq_id, '>' + seq_id)
            sequence = seqid_to_seq.get(seq_id, '')
            f.write(f"{header}\n")
            # Write sequence in lines of 60 characters
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
            print(f"Re-added sequence: {seq_id}")

def main():
    import sys
    if len(sys.argv) != 5:
        print("Usage: python script.py <original_fasta> <cdhit_fasta> <cluster_file> <output_fasta>")
        sys.exit(1)
    original_fasta = sys.argv[1]
    cdhit_fasta = sys.argv[2]
    cluster_file = sys.argv[3]
    output_fasta = sys.argv[4]

    # Step 1: Parse original fasta file
    seqid_to_seq, seqid_to_header, seqid_to_prefix = parse_fasta_with_prefix(original_fasta)

    # Step 2: Copy CD-HIT output fasta to new output fasta file
    with open(cdhit_fasta, 'r') as src, open(output_fasta, 'w') as dst:
        dst.write(src.read())

    # Step 3: Parse cluster file
    clusters = parse_cluster_file(cluster_file)

    # Step 4: Identify minority sequences with the extended logic
    minority_seq_ids = identify_minority_sequences(clusters)

    # Step 5: Write minority sequences to output fasta file
    write_sequences_to_file(minority_seq_ids, seqid_to_seq, seqid_to_header, output_fasta)

if __name__ == "__main__":
    main()
