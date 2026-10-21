# tRNA modification proteins across all domains of life

Snakemake-based bioinformatics pipeline for identifying, filtering, and
analyzing homologs of tRNA-modification enzymes across archaea, bacteria,
and eukaryotes. The pipeline builds curated protein datasets, runs homology
searches, constructs phylogenetic trees, and produces summary
visualizations (presence/absence tables, dot graphs, domain figures,
subcellular localization annotations)

## Repository Layout

```
pipeline/               Ordered pipeline steps (00 → 17), each a self-contained Snakemake module
```

Each numbered folder in `pipeline/` typically contains:
- One or more `*.smk` Snakemake workflow files
- Supporting Python scripts invoked by the rules
- A `config.yaml` / `*.yaml` file with run parameters
- `pixi.toml` / `pixi.lock` for reproducible environment management (via [pixi](https://pixi.sh))
- A `README.txt` with step-specific notes, where present

Run each step's Snakemake file from within its own directory (e.g. `cd pipeline/00-create-dataset && snakemake -s create-dataset.smk`), respecting the numeric order, since later steps consume outputs from earlier ones.

## Filtering Suffix Glossary

Filenames throughout `output/` encode which filtering stages a sequence set
has passed through:

| Suffix | Meaning |
|--------|---------|
| `k01`  | Filtered by KEGG orthology — sequences not matching the reference KEGG ID are removed |
| `f01`  | Filtered by domain architecture — sequences with missing/extra/aberrant domains relative to the reference profile are removed |
| `e01`  | Eukaryotic contamination removed via manual inspection of the phylogenetic tree |
| `c01`  | CD-HIT redundancy removal |
| `m01` / `m02` | Sequences removed after manual inspection of the phylogenetic tree (rounds 1 and 2) |
| `r01`  | Redundant eukaryotic paralogs removed — when multiple sequences from the same species cluster together, only the longest is kept |
| `a01`  | Sequences added back based on a curated ID list (KEGG-negative but below an e-value threshold, e.g. Tuc1 < 1e-25, TRM9B < 1e-20) |

## Pipeline Steps (00–17)

| Step | Folder | Purpose |
|------|--------|---------|
| 00 | `00-create-dataset` | Build the initial genome/protein dataset: download NCBI proteins and GTDB taxonomy, filter GTDB metadata, map taxa to genomes, and deduplicate proteins. |
| 01 | `01-order-color-dataset` | Assign a display order and color scheme to taxa; merge taxonomy metadata and build the protein-to-genome mapping used for downstream tables and trees. |
| 02 | `02-database` | Build local BLAST databases for archaea, bacteria, and eukaryota from the reference sequence sets. |
| 03 | `03-homolog-search` | BLAST each query protein against the archaeal, bacterial, and eukaryotic databases; extract hit sequences, annotate domains (InterPro/Pfam), and filter BLAST results by e-value and domain composition. |
| 04 | `04-merge` | Merge per-domain BLAST/annotation results across archaea, bacteria, and eukaryota into combined phylum-level tables. |
| 05 | `05-summary-stats` | Generate protein presence tables and summary statistics per KEGG reference, including highlighting reference KEGG hits. |
| 06 | `06-homolog-filter` | Core quality-control filtering: apply KEGG ID filtering (`k01`), domain-based filtering (`f01`), removal of eukaryotic contaminants (`e01`), CD-HIT deduplication (`c01`), and manual filtering rounds (`m01`, `m02`). Produces a filtering summary and verification report. |
| 07 | `07-tree` | Align filtered sequences, trim gapped alignment columns (trimAl), and build maximum-likelihood phylogenetic trees with IQ-TREE. |
| 08 | `08-combined-trees` | Combine/re-filter sequences across related protein families, regenerate alignments and trees, optionally annotating trees with domain information. |
| 09 | `09-filter-euk-copies` | Remove redundant/duplicate eukaryotic paralogs from trees, producing the final cleaned eukaryotic tree (`r01`-style filtering). |
| 10 | `10-special-filter-cases` | Handle protein-family-specific edge cases: add back specific archaeal/bacterial sequences (`a01`) based on curated ID lists and e-value thresholds. |
| 11 | `11-tree-Cmodels` | Collect and summarize the best-fit substitution models used for each protein's tree (from IQ-TREE ModelFinder), producing a models summary table and special-case tree configs. |
| 12 | `12-final-edit-tree` | Apply final manual edits (pruning/relabeling) to trees before downstream visualization. |
| 13 | `13-deeploc` | Extract eukaryotic sequences and run/parse DeepLoc subcellular localization predictions; filter and collapse localization calls (cytoplasm/nucleus/mitochondria/plastid) for iTOL annotation. |
| 14 | `14-itol-tree` | Prepare iTOL-compatible datasets, domain annotations, and label/template files for interactive tree visualization. |
| 15 | `15-PA-table` | Build a protein presence/absence (PA) table across taxa, order and color it by taxonomic group, and analyze protein co-occurrence patterns. |
| 16 | `16-dot-graph` | Generate dot-plot graphs (main and supplementary) visualizing protein presence/absence patterns, with FASTA comparison utilities. |
| 17 | `17-domain-figure` | Produce protein domain architecture figures for selected proteins (e.g. MnmE, TadA, TilS) based on per-protein configs. |

## Environment Management

Most step directories use [pixi](https://pixi.sh) (`pixi.toml`/`pixi.lock`) to
pin tool versions (BLAST, IQ-TREE, trimAl, CD-HIT, DeepLoc, etc.). Run
`pixi install` inside a step's directory before executing its Snakemake
workflow.

## Notes

See `notes/README.txt` for manual curation notes on ambiguous eukaryotic
paralog assignments encountered while filtering trees, and the
per-step `README.txt` files inside `pipeline/*/` for additional details
specific to that stage.
