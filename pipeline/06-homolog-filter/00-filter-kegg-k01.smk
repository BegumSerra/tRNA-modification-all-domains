"""
00-kegg-filter: Find IDs that are NOT in KEGG lists (to be filtered out)
Input: homology-analyses-tables/raw_tables/{protein}_table.tsv
       input/kegg/{protein}_kegg.txt
Output: seqs/filtering/ids/k01/{protein}_non_kegg.txt
"""

configfile: "01-config.yaml"
PROTEIN = config["PROTEIN"]
results = "../../output/"

rule all:
    input:
        expand("../../output/seqs/filtering/ids/k01/{protein}_non_kegg.txt", protein=PROTEIN),
        expand(results + "seqs/filtering/fasta/{protein}_k01.fa", protein=PROTEIN)

rule filter_kegg:
    input:
        table = "../../output/homology-analyses-tables/raw_tables/{protein}_table.tsv",
        kegg = "../../input/kegg/{protein}_kegg.txt"
    output:
        "../../output/seqs/filtering/ids/k01/{protein}_non_kegg.txt"
    shell:
        """
        python filter-kegg-ids.py {input.table} {input.kegg} {output}
        """

rule filter_kegg_fasta:
    input:
        fasta = results + "seqs/blast/merged/{protein}.fa",
        kegg_filter = results + "seqs/filtering/ids/k01/{protein}_non_kegg.txt"
    output:
        results + "seqs/filtering/fasta/{protein}_k01.fa"
    shell:
        "python fasta-remove.py {input.fasta} {input.kegg_filter} {output}"
