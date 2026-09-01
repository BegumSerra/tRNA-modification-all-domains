# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

configfile: "config.yaml"
wildcard_constraints:
    protein_arc = "|".join(config["PROTEIN_ARC"]),
    protein_bact = "|".join(config["PROTEIN_BACT"]),
    protein_euk = "|".join(config["PROTEIN_EUK"])

results = "../../output/"
data = "../../input/"
database = "../../input/database/"
logs = "logs/"

PROTEIN_ARC = config["PROTEIN_ARC"]
PROTEIN_BACT = config["PROTEIN_BACT"]
PROTEIN_EUK = config["PROTEIN_EUK"]

rule all:
    input:
        expand(results + "domain/ref_seq/arc/{protein_arc}_ref_arc.tsv", protein_arc=PROTEIN_ARC),
        expand(results + "domain/ref_seq/bact/{protein_bact}_ref_bact.tsv", protein_bact=PROTEIN_BACT),
        expand(results + "domain/ref_seq/euk/{protein_euk}_ref_euk.tsv", protein_euk=PROTEIN_EUK),

        expand(results + "domain/blast/euk-arc/{protein_euk}_filtered_length.tsv", protein_euk=PROTEIN_EUK),
        expand(results + "domain/blast/bact-arc/{protein_bact}_filtered_length.tsv", protein_bact=PROTEIN_BACT),
        expand(results + "domain/blast/arc-arc/{protein_arc}_filtered_length.tsv", protein_arc=PROTEIN_ARC),

        expand(results + "blast/blast/euk-arc/{protein_euk}_unique_filtered.tab", protein_euk=PROTEIN_EUK),
        expand(results + "blast/blast/bact-arc/{protein_bact}_unique_filtered.tab", protein_bact=PROTEIN_BACT),
        expand(results + "blast/blast/arc-arc/{protein_arc}_unique_filtered.tab", protein_arc=PROTEIN_ARC)


rule interpro_ref_arc:
    input:
        data + "ref_seqs/arc/{protein_arc}.fa"
    output:
        results + "domain/ref_seq/arc/{protein_arc}_ref_arc.tsv"
    log:
        logs + "interpro_ref_arc_{protein_arc}.log"
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

rule interpro_ref_bact:
    input:
        data + "ref_seqs/bact/{protein}.fa"
    output:
        results + "domain/ref_seq/bact/{protein}_ref_bact.tsv"
    log:
        logs + "interpro_ref_bact_{protein}.log"
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

rule interpro_ref_euk:
    input:
        data + "ref_seqs/euk/{protein}.fa"
    output:
        results + "domain/ref_seq/euk/{protein}_ref_euk.tsv"
    log:
        logs + "interpro_ref_euk_{protein}.log"
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

### ALL to Archaea database BLAST Search  ###

rule blastp_euk_arc:
    input:
        seqs = data + "ref_seqs/euk/{protein_euk}.fa"
    output:
        results + "blast/blast/euk-arc/{protein_euk}.tab"
    params:
        db = results + "blastdb/archaea/blastdb"
    log:
        logs + "blastp_euk_arc_{protein_euk}.log"
    shell:
        """
        blastp -query {input.seqs} -db {params.db} -evalue 0.001 \
                 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
                 -max_target_seqs 1000 -out {output} 2>&1 | tee {log}
        """

rule filter_duplicates_euk_arc:
    input:
        results + "blast/blast/euk-arc/{protein_euk}.tab"
    output:
        results + "blast/blast/euk-arc/{protein_euk}_unique.tab"
    shell:
        "awk '!seen[$2]++' {input}  > {output}"

rule extract_fasta_euk_arc:
    input:
        blast = results + "blast/blast/euk-arc/{protein_euk}_unique.tab",
        proteome = database + "archaea/archaea_ncbi_proteins_unique.faa"
    output:
        results + "seqs/blast/euk-arc/{protein_euk}_blast.fa"
    log:
        logs + "extract_fasta_euk_arc_{protein_euk}.log"
    shell:
        "python extract-fasta-blast.py {input.blast} {input.proteome} {output} 2>&1 | tee {log}"

rule remove_asteriks_euk_arc:
    input:
        results + "seqs/blast/euk-arc/{protein_euk}_blast.fa"
    output:
        results + "seqs/blast/euk-arc/{protein_euk}_cleaned.fa"
    shell:
        """
        awk '/^>/ {{print; next}} {{gsub(/\\*/, "", $0); print}}' {input} > {output}
        """

rule interpro_euk_arc:  # --resources interpro_limit=1
    input:
        results + "seqs/blast/euk-arc/{protein_euk}_cleaned.fa"
    output:
        results + "domain/blast/euk-arc/{protein_euk}.tsv"
    log:
        logs + "interpro_euk_arc_{protein_euk}.log"
    resources:
        interpro_limit=1
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

rule filter_pfam_euk_arc:
    input:
        ref=results + "domain/ref_seq/euk/{protein_euk}_ref_euk.tsv",
        query=results + "domain/blast/euk-arc/{protein_euk}.tsv"
    output:
        results + "domain/blast/euk-arc/{protein_euk}_filtered.tsv"
    log:
        logs + "filter_pfam_euk_arc_{protein_euk}.log"
    shell:
        "python filter-pfam-EBA.py {input.ref} {input.query} {output} 2>&1 | tee {log}"

rule filter_blast_euk_arc:
    input:
        domain=results + "domain/blast/euk-arc/{protein_euk}_filtered.tsv",
        blast=results + "blast/blast/euk-arc/{protein_euk}_unique.tab"
    output:
        results + "blast/blast/euk-arc/{protein_euk}_unique_filtered.tab",
    log:
        logs + "filter_blast_euk_arc_{protein_euk}.log"
    shell:
        "python filter-blast.py {input.domain} {input.blast} {output} 2>&1 | tee {log}"

rule length_pfam_euk_arc:
    input:
        results + "domain/blast/euk-arc/{protein_euk}_filtered.tsv",
    output:
        results + "domain/blast/euk-arc/{protein_euk}_filtered_length.tsv"
    log:
        logs + "length_pfam_euk_arc_{protein_euk}.log"
    shell:
        "python domain-length-column.py {input} {output} 2>&1 | tee {log}"

rule blastp_bact_arc:
    input:
        seqs=data + "ref_seqs/bact/{protein_bact}.fa"
    output:
        results + "blast/blast/bact-arc/{protein_bact}.tab"
    params:
        db = results + "blastdb/archaea/blastdb"
    log:
        logs + "blastp_bact_arc_{protein_bact}.log"
    shell:
        """
        blastp -query {input.seqs} -db {params.db} -evalue 0.001 \
                 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
                 -max_target_seqs 1000 -out {output} 2>&1 | tee {log}
        """

rule filter_duplicates_bact_arc:
    input:
        results + "blast/blast/bact-arc/{protein_bact}.tab"
    output:
        results + "blast/blast/bact-arc/{protein_bact}_unique.tab"
    shell:
        "awk '!seen[$2]++' {input}  > {output}"

rule extract_fasta_bact_arc:
    input:
        blast = results + "blast/blast/bact-arc/{protein_bact}_unique.tab",
        proteome = database + "archaea/archaea_ncbi_proteins_unique.faa"
    output:
        results + "seqs/blast/bact-arc/{protein_bact}_blast.fa"
    log:
        logs + "extract_fasta_bact_arc_{protein_bact}.log"
    shell:
        "python extract-fasta-blast.py {input.blast} {input.proteome} {output} 2>&1 | tee {log}"

rule remove_asteriks_bact_arc:
    input:
        results + "seqs/blast/bact-arc/{protein_bact}_blast.fa"
    output:
        results + "seqs/blast/bact-arc/{protein_bact}_cleaned.fa"
    shell:
        """
        awk '/^>/ {{print; next}} {{gsub(/\\*/, "", $0); print}}' {input} > {output}
        """

rule interpro_bact_arc:  # --resources interpro_limit=1
    input:
        results + "seqs/blast/bact-arc/{protein_bact}_cleaned.fa"
    output:
        results + "domain/blast/bact-arc/{protein_bact}.tsv"
    log:
        logs + "interpro_bact_arc_{protein_bact}.log"
    resources:
        interpro_limit=1
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

rule filter_pfam_bact_arc:
    input:
        ref=results + "domain/ref_seq/bact/{protein_bact}_ref_bact.tsv",
        query=results + "domain/blast/bact-arc/{protein_bact}.tsv"
    output:
        results + "domain/blast/bact-arc/{protein_bact}_filtered.tsv"
    log:
        logs + "filter_pfam_bact_arc_{protein_bact}.log"
    shell:
        "python filter-pfam-EBA.py {input.ref} {input.query} {output} 2>&1 | tee {log}"

rule filter_blast_bact_arc:
    input:
        domain=results + "domain/blast/bact-arc/{protein_bact}_filtered.tsv",
        blast=results + "blast/blast/bact-arc/{protein_bact}_unique.tab"
    output:
        results + "blast/blast/bact-arc/{protein_bact}_unique_filtered.tab",
    log:
        logs + "filter_blast_bact_arc_{protein_bact}.log"
    shell:
        "python filter-blast.py {input.domain} {input.blast} {output} 2>&1 | tee {log}"

rule length_pfam_bact_arc:
    input:
        results + "domain/blast/bact-arc/{protein_bact}_filtered.tsv",
    output:
        results + "domain/blast/bact-arc/{protein_bact}_filtered_length.tsv"
    log:
        logs + "length_pfam_bact_arc_{protein_bact}.log"
    shell:
        "python domain-length-column.py {input} {output} 2>&1 | tee {log}"

rule blastp_arc_arc:
    input:
        seqs = data + "ref_seqs/arc/{protein_arc}.fa"
    output:
        results + "blast/blast/arc-arc/{protein_arc}.tab"
    params:
        db = results + "blastdb/archaea/blastdb"
    log:
        logs + "blastp_arc_arc_{protein_arc}.log"
    shell:
        """
        blastp -query {input.seqs} -db {params.db} -evalue 0.001 \
                 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
                 -max_target_seqs 1000 -out {output} 2>&1 | tee {log}
        """

rule filter_duplicates_arc_arc:
    input:
        results + "blast/blast/arc-arc/{protein_arc}.tab"
    output:
        results + "blast/blast/arc-arc/{protein_arc}_unique.tab"
    shell:
        "awk '!seen[$2]++' {input}  > {output}"

rule extract_fasta_arc_arc:
    input:
        blast = results + "blast/blast/arc-arc/{protein_arc}_unique.tab",
        proteome = database + "archaea/archaea_ncbi_proteins_unique.faa"
    output:
        results + "seqs/blast/arc-arc/{protein_arc}_blast.fa"
    log:
        logs + "extract_fasta_arc_arc_{protein_arc}.log"
    shell:
        "python extract-fasta-blast.py {input.blast} {input.proteome} {output} 2>&1 | tee {log}"

rule remove_asteriks_arc_arc:
    input:
        results + "seqs/blast/arc-arc/{protein_arc}_blast.fa"
    output:
        results + "seqs/blast/arc-arc/{protein_arc}_cleaned.fa"
    shell:
        """
        awk '/^>/ {{print; next}} {{gsub(/\\*/, "", $0); print}}' {input} > {output}
        """

rule interpro_arc_arc:  # --resources interpro_limit=1
    input:
        results + "seqs/blast/arc-arc/{protein_arc}_cleaned.fa"
    output:
        results + "domain/blast/arc-arc/{protein_arc}.tsv"
    log:
        logs + "interpro_arc_arc_{protein_arc}.log"
    resources:
        interpro_limit=1
    shell:
        "interproscan.sh -i {input} -o {output} -f TSV -cpu 50 -appl Pfam"

rule filter_pfam_arc_arc:
    input:
        ref=results + "domain/ref_seq/arc/{protein_arc}_ref_arc.tsv",
        query=results + "domain/blast/arc-arc/{protein_arc}.tsv"
    output:
        results + "domain/blast/arc-arc/{protein_arc}_filtered.tsv"
    log:
        logs + "filter_pfam_arc_arc_{protein_arc}.log"
    shell:
        "python filter-pfam-EBA.py {input.ref} {input.query} {output} 2>&1 | tee {log}"

rule filter_blast_arc_arc:
    input:
        domain=results + "domain/blast/arc-arc/{protein_arc}_filtered.tsv",
        blast=results + "blast/blast/arc-arc/{protein_arc}_unique.tab"
    output:
        results + "blast/blast/arc-arc/{protein_arc}_unique_filtered.tab",
    log:
        logs + "filter_blast_arc_arc_{protein_arc}.log"
    shell:
        "python filter-blast.py {input.domain} {input.blast} {output} 2>&1 | tee {log}"

rule length_pfam_arc_arc:
    input:
        results + "domain/blast/arc-arc/{protein_arc}_filtered.tsv",
    output:
        results + "domain/blast/arc-arc/{protein_arc}_filtered_length.tsv"
    log:
        logs + "length_pfam_arc_arc_{protein_arc}.log"
    shell:
        "python domain-length-column.py {input} {output} 2>&1 | tee {log}"

