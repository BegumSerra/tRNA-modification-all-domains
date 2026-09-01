# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "03-config.yaml"

wildcard_constraints:
    protein1 = "|".join(config["PROTEIN1"]),
    protein2 = "|".join(config["PROTEIN2"]),
    protein3 = "|".join(config["PROTEIN3"]),
    protein4 = "|".join(config["PROTEIN4"])

results = "../../output/"
data = "../../input/"
database = "../../database/"

PROTEIN1 = config["PROTEIN1"]
PROTEIN2 = config["PROTEIN2"]
PROTEIN3 = config["PROTEIN3"]
PROTEIN4 = config["PROTEIN4"]

rule all:
    input:
        expand(results + "seqs/blast/merged/{protein1}.fa", protein1=PROTEIN1),
        expand(results + "seqs/blast/merged/{protein2}.fa", protein2=PROTEIN2),
        expand(results + "seqs/blast/merged/{protein3}.fa", protein3=PROTEIN3),
        expand(results + "seqs/blast/merged/{protein4}.fa", protein4=PROTEIN4),
        expand(results + "domain/blast/merged/{protein1}.tsv", protein1=PROTEIN1),
        expand(results + "domain/blast/merged/{protein2}.tsv", protein2=PROTEIN2),
        expand(results + "domain/blast/merged/{protein3}.tsv", protein3=PROTEIN3),
        expand(results + "domain/blast/merged/{protein4}.tsv", protein4=PROTEIN4),
        expand(results + "blast/blast/merged/{protein1}.tab", protein1=PROTEIN1),
        expand(results + "blast/blast/merged/{protein2}.tab", protein2=PROTEIN2),
        expand(results + "blast/blast/merged/{protein3}.tab", protein3=PROTEIN3),
        expand(results + "blast/blast/merged/{protein4}.tab", protein4=PROTEIN4)

rule merge_euk_eba:
    input:
        eukbact = results + "seqs/blast/euk-bact/{protein1}_cleaned.fa",
        eukeuk = results + "seqs/blast/euk-euk/{protein1}_cleaned.fa",
        eukarc = results + "seqs/blast/euk-arc/{protein1}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein1}.fa"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule merge_euk_eba_domain:
    input:
        eukbact = results + "domain/blast/euk-bact/{protein1}.tsv",
        eukeuk = results + "domain/blast/euk-euk/{protein1}.tsv",
        eukarc = results + "domain/blast/euk-arc/{protein1}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein1}.tsv"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule merge_euk_eba_blast:
    input:
        eukbact = results + "blast/blast/euk-bact/{protein1}.tab",
        eukeuk = results + "blast/blast/euk-euk/{protein1}.tab",
        eukarc = results + "blast/blast/euk-arc/{protein1}.tab"
    output:
        merged = results + "blast/blast/merged/{protein1}.tab"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule merge_euk_eb:
    input:
        eukbact = results +"seqs/blast/euk-bact/{protein2}_cleaned.fa",
        eukeuk = results +"seqs/blast/euk-euk/{protein2}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein2}.fa"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} > {output.merged}
        """

rule merge_euk_eb_domain:
    input:
        eukbact = results + "domain/blast/euk-bact/{protein2}.tsv",
        eukeuk = results + "domain/blast/euk-euk/{protein2}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein2}.tsv"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} > {output.merged}
        """

rule merge_euk_eb_blast:
    input:
        eukbact = results + "blast/blast/euk-bact/{protein2}.tab",
        eukeuk = results + "blast/blast/euk-euk/{protein2}.tab"
    output:
        merged = results + "blast/blast/merged/{protein2}.tab"
    shell:
        """
        cat {input.eukbact} {input.eukeuk} > {output.merged}
        """

rule merge_euk_ea:
    input:
        eukeuk = results + "seqs/blast/euk-euk/{protein3}_cleaned.fa",
        eukarc = results + "seqs/blast/euk-arc/{protein3}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein3}.fa"
    shell:
        """
        cat {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule merge_euk_ea_domain:
    input:
        eukeuk = results + "domain/blast/euk-euk/{protein3}.tsv",
        eukarc = results + "domain/blast/euk-arc/{protein3}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein3}.tsv"
    shell:
        """
        cat {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule merge_euk_ea_blast:
    input:
        eukeuk = results + "blast/blast/euk-euk/{protein3}.tab",
        eukarc = results + "blast/blast/euk-arc/{protein3}.tab"
    output:
        merged = results + "blast/blast/merged/{protein3}.tab"
    shell:
        """
        cat {input.eukeuk} {input.eukarc} > {output.merged}
        """

rule copy_euk_euk:
    input:
        eukeuk = results + "seqs/blast/euk-euk/{protein4}_cleaned.fa"
    output:
        merged = results + "seqs/blast/merged/{protein4}.fa"
    shell:
        """
        cp {input.eukeuk} {output.merged}
        """

rule copy_euk_euk_domain:
    input:
        eukeuk = results + "domain/blast/euk-euk/{protein4}.tsv"
    output:
        merged = results + "domain/blast/merged/{protein4}.tsv"
    shell:
        """
        cp {input.eukeuk} {output.merged}
        """

rule copy_euk_euk_blast:
    input:
        eukeuk = results + "blast/blast/euk-euk/{protein4}.tab"
    output:
        merged = results + "blast/blast/merged/{protein4}.tab"
    shell:
        """
        cp {input.eukeuk} {output.merged}
        """