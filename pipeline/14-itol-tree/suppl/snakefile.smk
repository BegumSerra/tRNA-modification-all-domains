# Snakefile for generating iTOL taxonomy range files
#
# Usage:
#   snakemake -s snakefile.smk -c1              (runs all proteins in config.yaml)
#   snakemake -s snakefile.smk -c1 --configfile config.yaml
#
# config.yaml format:
#   proteins:
#     Elp3: k01_f01_e00_c01_m01_nogaps_r01
#     ADAT2: k01_f01_e01_c01_nogaps_r01

configfile: "config.yaml"

# ── config ────────────────────────────────────────────────────────────
PROTEINS          = config.get("proteins", {})           # single-protein trees
COMBINED_PROTEINS = config.get("combined_proteins", {})  # merged / combined trees

# Unified dict of ALL proteins (single + combined) for wildcard iteration
ALL_PROTEINS = {}
ALL_PROTEINS.update(PROTEINS)
ALL_PROTEINS.update(COMBINED_PROTEINS)

# ── paths ─────────────────────────────────────────────────────────────
BASE_DIR    = "/home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA"
OUTPUT_DIR  = BASE_DIR + "/output"
SCRIPT_DIR  = BASE_DIR + "/pipeline/14-itol-tree/suppl"
COLOR_TSV   = BASE_DIR + "/input/table/color.tsv"

GRADIENT_DIR = OUTPUT_DIR + "/itol/gradient"
LABELS_DIR   = OUTPUT_DIR + "/itol/templates/labels"

# ── collect all outputs across proteins ──────────────────────────────
ALL_OUTPUTS = []
for p in ALL_PROTEINS:
    ALL_OUTPUTS += [
        GRADIENT_DIR + "/" + p + "_archaea_gradient_itol.txt",
        GRADIENT_DIR + "/" + p + "_eukaryota_gradient_itol.txt",
        GRADIENT_DIR + "/" + p + "_bacteria_gradient_itol.txt",
        GRADIENT_DIR + "/" + p + "_archaea_group_text_right_itol.txt",
        GRADIENT_DIR + "/" + p + "_eukaryota_group_text_right_itol.txt",
        GRADIENT_DIR + "/" + p + "_bacteria_group_text_right_itol.txt",
        LABELS_DIR   + "/" + p + "_long_names_labels.txt",
    ]

# ── helper: look up treefile / phylum-table path for a protein wildcard ──
#   Combined trees live under  combined_tree/  subdirectories;
#   single-protein trees live directly under  iqtree/  and  phylum_table/.
def treefile_for(wildcards):
    p = wildcards.protein
    tid = ALL_PROTEINS[p]
    if p in COMBINED_PROTEINS:
        return OUTPUT_DIR + "/iqtree/combined_tree/" + p + "_" + tid + "_Cmod.treefile"
    return OUTPUT_DIR + "/iqtree/" + p + "_" + tid + "_Cmod.treefile"

def phylum_table_for(wildcards):
    p = wildcards.protein
    if p in COMBINED_PROTEINS:
        return OUTPUT_DIR + "/table/phylum_table/combined_tree/" + p + "_phylum.table.tsv"
    return OUTPUT_DIR + "/table/phylum_table/" + p + "_phylum.table.tsv"

# ── rules ─────────────────────────────────────────────────────────────
rule all:
    input:
        ALL_OUTPUTS

rule domain_gradient_ranges:
    input:
        treefile = treefile_for,
        phylum_table = phylum_table_for,
        color_tsv = COLOR_TSV
    output:
        arc  = GRADIENT_DIR + "/{protein}_archaea_gradient_itol.txt",
        euk  = GRADIENT_DIR + "/{protein}_eukaryota_gradient_itol.txt",
        bact = GRADIENT_DIR + "/{protein}_bacteria_gradient_itol.txt"
    shell:
        """
        python {SCRIPT_DIR}/all-domain-group-ranges-gradient.py \
            {input.treefile} {input.phylum_table} {input.color_tsv} \
            {output.arc} {output.euk} {output.bact}
        """

rule archaea_text_labels:
    input:
        gradient = GRADIENT_DIR + "/{protein}_archaea_gradient_itol.txt"
    output:
        GRADIENT_DIR + "/{protein}_archaea_group_text_right_itol.txt"
    shell:
        """
        python {SCRIPT_DIR}/all-domain-groups-labels.py {input.gradient} Archaea {output}
        """

rule eukaryota_text_labels:
    input:
        gradient = GRADIENT_DIR + "/{protein}_eukaryota_gradient_itol.txt"
    output:
        GRADIENT_DIR + "/{protein}_eukaryota_group_text_right_itol.txt"
    shell:
        """
        python {SCRIPT_DIR}/all-domain-groups-labels.py {input.gradient} Eukaryota {output}
        """

rule bacteria_text_labels:
    input:
        gradient = GRADIENT_DIR + "/{protein}_bacteria_gradient_itol.txt"
    output:
        GRADIENT_DIR + "/{protein}_bacteria_group_text_right_itol.txt"
    shell:
        """
        python {SCRIPT_DIR}/all-domain-groups-labels.py {input.gradient} Bacteria {output}
        """

rule long_name_labels:
    input:
        treefile = treefile_for,
        phylum_table = phylum_table_for
    output:
        LABELS_DIR + "/{protein}_long_names_labels.txt"
    shell:
        """
        python {SCRIPT_DIR}/all-domain-branch-long-names.py {input.treefile} {input.phylum_table} {output}
        """
