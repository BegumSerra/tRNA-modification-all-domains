# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "02-config.yaml"

wildcard_constraints:
    protein1 = "|".join(config["PROTEIN1"]),
    protein2 = "|".join(config["PROTEIN2"]),
    protein3 = "|".join(config["PROTEIN3"])

results = "../../output/"
data = "../../input/"
database = "../../database/"

PROTEIN1 = config["PROTEIN1"]
PROTEIN2 = config["PROTEIN2"]
PROTEIN3 = config["PROTEIN3"]

rule all:
    input:
        expand(results + "seqs/blast/merged/{protein1}.fa", protein1=PROTEIN1),
        expand(results + "seqs/blast/merged/{protein2}.fa", protein2=PROTEIN2),
        expand(results + "seqs/blast/merged/{protein3}.fa", protein3=PROTEIN3),
        expand(results + "domain/blast/merged/{protein1}.tsv", protein1=PROTEIN1),
        expand(results + "domain/blast/merged/{protein2}.tsv", protein2=PROTEIN2),
        expand(results + "domain/blast/merged/{protein3}.tsv", protein3=PROTEIN3),
        expand(results + "blast/blast/merged/{protein1}.tab", protein1=PROTEIN1),
        expand(results + "blast/blast/merged/{protein2}.tab", protein2=PROTEIN2),
        expand(results + "blast/blast/merged/{protein3}.tab", protein3=PROTEIN3)

rule merge_bact_eba:
    input:
        bactbact = results +"seqs/blast/bact-bact/{protein1}_cleaned.fa",
        bacteuk = results +"seqs/blast/bact-euk/{protein1}_cleaned.fa",
        bactarc = results +"seqs/blast/bact-arc/{protein1}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein1}.fa"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} {input.bactarc} > {output.merged}
        """

rule merge_bact_eba_domain:
    input:
        bactbact = results + "domain/blast/bact-bact/{protein1}.tsv",
        bacteuk = results + "domain/blast/bact-euk/{protein1}.tsv",
        bactarc = results + "domain/blast/bact-arc/{protein1}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein1}.tsv"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} {input.bactarc} > {output.merged}
        """

rule merge_bact_eba_blast:
    input:
        bactbact = results + "blast/blast/bact-bact/{protein1}.tab",
        bacteuk = results + "blast/blast/bact-euk/{protein1}.tab",
        bactarc = results + "blast/blast/bact-arc/{protein1}.tab"
    output:
        merged = results + "blast/blast/merged/{protein1}.tab"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} {input.bactarc} > {output.merged}
        """

rule merge_bact_ba:
    input:
        bactbact = results +"seqs/blast/bact-bact/{protein2}_cleaned.fa",
        bactarc = results +"seqs/blast/bact-arc/{protein2}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein2}.fa"
    shell:
        """
        cat {input.bactbact} {input.bactarc} > {output.merged}
        """

rule merge_bact_ba_domain:
    input:
        bactbact = results + "domain/blast/bact-bact/{protein2}.tsv",
        bactarc = results + "domain/blast/bact-arc/{protein2}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein2}.tsv"
    shell:
        """
        cat {input.bactbact} {input.bactarc} > {output.merged}
        """

rule merge_bact_ba_blast:
    input:
        bactbact = results + "blast/blast/bact-bact/{protein2}.tab",
        bactarc = results + "blast/blast/bact-arc/{protein2}.tab"
    output:
        merged = results + "blast/blast/merged/{protein2}.tab"
    shell:
        """
        cat {input.bactbact} {input.bactarc} > {output.merged}
        """

rule merge_bact_eb:
    input:
        bactbact = results +"seqs/blast/bact-bact/{protein3}_cleaned.fa",
        bacteuk = results +"seqs/blast/bact-euk/{protein3}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein3}.fa"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} > {output.merged}
        """

rule merge_bact_eb_domain:
    input:
        bactbact = results + "domain/blast/bact-bact/{protein3}.tsv",
        bacteuk = results + "domain/blast/bact-euk/{protein3}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein3}.tsv"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} > {output.merged}
        """

rule merge_bact_eb_blast:
    input:
        bactbact = results + "blast/blast/bact-bact/{protein3}.tab",
        bacteuk = results + "blast/blast/bact-euk/{protein3}.tab"
    output:
        merged = results + "blast/blast/merged/{protein3}.tab"
    shell:
        """
        cat {input.bactbact} {input.bacteuk} > {output.merged}
        """