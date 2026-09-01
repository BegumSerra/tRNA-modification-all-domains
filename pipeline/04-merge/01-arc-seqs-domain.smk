# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "01-config.yaml"
wildcard_constraints:
    protein1 = "|".join(config["PROTEIN1"]),
    protein2 = "|".join(config["PROTEIN2"])

results = "../../output/"
data = "../../input/"
database = "../../database/"

PROTEIN1 = config["PROTEIN1"]
PROTEIN2 = config["PROTEIN2"]

rule all:
    input:
        expand(results + "seqs/blast/merged/{protein1}.fa", protein1=PROTEIN1),
        expand(results + "seqs/blast/merged/{protein2}.fa", protein2=PROTEIN2),
        expand(results + "domain/blast/merged/{protein1}.tsv", protein1=PROTEIN1),
        expand(results + "domain/blast/merged/{protein2}.tsv", protein2=PROTEIN2),
        expand(results + "blast/blast/merged/{protein1}.tab", protein1=PROTEIN1),
        expand(results + "blast/blast/merged/{protein2}.tab", protein2=PROTEIN2)

rule merge_arc_eba:
    input:
        arcbact = results +"seqs/blast/arc-bact/{protein1}_cleaned.fa",
        arceuk = results +"seqs/blast/arc-euk/{protein1}_cleaned.fa",
        arcarc = results +"seqs/blast/arc-arc/{protein1}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein1}.fa"
    shell:
        """
        cat {input.arcbact} {input.arceuk} {input.arcarc} > {output.merged}
        """

rule copy_arc_arc:
    input:
        arcarc = results + "seqs/blast/arc-arc/{protein2}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein2}.fa"
    shell:
        """
        cp {input.arcarc} {output.merged}
        """

rule merge_arc_eba_domain:
    input:
        arcbact = results + "domain/blast/arc-bact/{protein1}.tsv",
        arceuk = results + "domain/blast/arc-euk/{protein1}.tsv",
        arcarc = results + "domain/blast/arc-arc/{protein1}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein1}.tsv"
    shell:
        """
        cat {input.arcbact} {input.arceuk} {input.arcarc} > {output.merged}
        """

rule copy_arc_arc_domain:
    input:
        arcarc = results +"domain/blast/arc-arc/{protein2}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein2}.tsv"
    shell:
        """
        cp {input.arcarc} {output.merged}
        """

rule merge_arc_eba_blast:
    input:
        arcbact = results + "blast/blast/arc-bact/{protein1}.tab",
        arceuk = results + "blast/blast/arc-euk/{protein1}.tab",
        arcarc = results + "blast/blast/arc-arc/{protein1}.tab"
    output:
        merged = results + "blast/blast/merged/{protein1}.tab"
    shell:
        """
        cat {input.arcbact} {input.arceuk} {input.arcarc} > {output.merged}
        """

rule copy_arc_arc_blast:
    input:
        arcarc = results +"blast/blast/arc-arc/{protein2}.tab"
    output:
        merged = results + "blast/blast/merged/{protein2}.tab"
    shell:
        """
        cp {input.arcarc} {output.merged}
        """
