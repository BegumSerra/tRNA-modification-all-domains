###The definitions of different filtering suffixes.###

--k01--
filtered sequences based on kegg. The IDs that don't match with the reference KEGG id will be filtered out

--f01--
filtered sequences based domain information if they dont match with the reference domain structures. See below "01-filterextract-ids.smk Filtering Overview" for more information.

--e01--
removing eukaryotic contaminations by manual detection from analysing the phylogenetic tree

--c01--
cdhit removal

--m01--
filtered sequences after by manual detection from analysing the phylogenetic tree

--m02--
filtered sequences after by manual detection from analysing the phylogenetic tree

--r01--
filtered eukaryotic sequences. When multiple eukaryotic sequences from the same species cluster together, keep only the longest one.

--a01--
Add archaeal sequences based on a01 ID list.
This adds sequences with kegg=False but e-value below specified threshold.
- Tuc1: e-value < 1e-25
- TRM9B: e-value < 1e-20


### 01-filterextract-ids.smk ###

F01 Filtering Overview
The f01 filtering is a quality control step that identifies and removes redundant or problematic protein sequences from phylogenetic trees based on domain architecture. It consists of two complementary filters:

Filter 1: filter_wo_domain (without domain)
Removes protein copies that lack domain information when a species has multiple copies:

For species with multiple protein copies in the tree
Identifies the most common domain profile across all proteins
If a species has at least one copy with the regular domain profile
Removes copies that have no domain information (missing domains)
Purpose: Eliminates incomplete sequences that lack domain annotations

Filter 2: filter_w_extra_domain (with extra domains)
Removes protein copies with aberrant domain architectures in multi-copy species:

For species with multiple copies, it removes proteins that:

Have extra domains: Contains all regular domains PLUS additional ones
Have completely different domains: No overlap with the regular domain profile
Keeps proteins that:
Match the regular domain profile exactly
Have multiple copies of regular domains (e.g., 2x the expected domain count)
Purpose: Removes fusion proteins, chimeras, or misannotated sequences with unusual domain architectures

Output
Creates merged list of all protein IDs to remove (*_merged.txt)
Generates 02-filter-from-tree.yaml listing proteins that had sequences filtered out
These filtered IDs are used in downstream steps to clean the phylogenetic trees
This filtering ensures your phylogenetic analysis uses high-quality, properly annotated sequences with expected domain architectures


