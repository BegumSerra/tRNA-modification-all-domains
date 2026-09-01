All model information printed to ../../output/iqtree/TusB.model.gz
CPU time for ModelFinder: 0.068 seconds (0h:0m:0s)
Wall-clock time for ModelFinder: 0.072 seconds (0h:0m:0s)
Generating 1000 samples for ultrafast bootstrap (seed: 239281)...
ERROR: It makes no sense to perform bootstrap with less than 4 sequences.

Building DAG of jobs...
MissingInputException in line 133 of /home/begum/PycharmProjects/PythonProject/Nov_2025_tRNA/pipeline/05-tree/tree.smk:
Missing input files for rule phylum_table:
    output: ../../output/table/phylum_table/MnmH_phylum.table.tsv
    wildcards: protein=MnmH
    affected files:
        ../../output/seqs/blast/merged/MnmH.fa

[Wed Jan  7 10:12:29 2026]
Error in rule tree:
    jobid: 238
    output: ../../output/iqtree/TusB.treefile
    shell:
        iqtree -s ../../output/alignments/trimmed_alignments/TusB_nogaps.trimal -m MFP -mset LG -bb 1000 --prefix ../../output/iqtree/TusB -redo
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2026-01-07T101228.465546.snakemake.log