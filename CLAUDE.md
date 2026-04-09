# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PathFX is a bioinformatics tool that analyzes drug-protein interactions and phenotype associations. Given a drug name/ID and/or gene targets, it extracts a protein-protein interaction subnetwork and performs hypergeometric enrichment analysis to identify associated phenotypes. Optionally, phenotypes can be clustered by semantic similarity using UMLS.

Published: Wilson et al, PLoS Comp Bio, 2018.

## Running PathFX

Scripts must be run from the `scripts/` directory — they use relative paths to `../rscs/` and `../results/`.

```bash
cd scripts/

# Basic analysis (drug targets fetched from DrugBank automatically)
python phenotype_enrichment_pathway_Pfx050120.py -d "DrugName" -a "analysis_name"

# With explicit gene targets
python phenotype_enrichment_pathway_Pfx050120.py -d "DrugName" -a "analysis_name" -t "GENE1,GENE2"

# With phenotype clustering (requires UMLS setup — see install.md)
python phenotype_enrichment_pathway_Pfx050120.py -d "DrugName" -a "analysis_name" -c True

# Standalone phenotype clustering (no PathFX run needed)
python run_stand_alone_phen_clust.py          # without real-time UMLS
python run_stand_alone_phen_clust_ex2.py      # with real-time UMLS
```

## Environment Setup

```bash
conda env create -f environment.yml
conda activate pathfx   # or whichever name is set in environment.yml
```

Key dependencies: numpy, pandas, scipy, matplotlib, seaborn, networkx, fastcluster.

For phenotype clustering: UMLS 2020AA database + Perl modules (`UMLS::Interface`, `UMLS::Similarity`). See `install.md`.

## Architecture

### Data Flow

```
Drug name/ID → DrugBank lookup (pfxDB050620_dint.pkl)
    → Get targets → Hash-based neighborhood lookup (pfx041520_*.pkl)
    → Subnetwork extraction
    → Hypergeometric enrichment (get_network_associations_Pfx050120.py)
    → BH multiple testing correction + background p-value filter
    → Association table output
    → [Optional] UMLS Lin similarity (calc_lin_matrix_umls.py via Perl)
    → [Optional] Ward hierarchical clustering + dendrogram (plot_and_cluster_phenotypes.py)
```

### Key Scripts

| Script | Role |
|--------|------|
| `phenotype_enrichment_pathway_Pfx050120.py` | Main orchestrator — CLI entry point |
| `get_network_associations_Pfx050120.py` | Hypergeometric enrichment + BH correction |
| `calc_lin_matrix_umls.py` | Computes UMLS Lin semantic similarity matrix via Perl |
| `plot_and_cluster_phenotypes.py` | Ward clustering + dendrogram generation |
| `get_network_associations_v3.py` | Legacy version of the association script |

### Resource Files (`rscs/`)

All pre-computed data is stored as pickled dictionaries (~247MB total):

- **Drug data**: `pfxDB050620_dint.pkl` (DrugBank ID → targets), `pfxDB050620_dbid2name.pkl`
- **Network**: `pfx041520_0.82_spec_nbhd_hash.pkl` (hash ID → neighborhood), `pfx041520_0.8_node_to_hashID.pkl`
- **Gene–phenotype mappings**: `Pfx050120_merged_genes_to_cuis.pkl`, `Pfx050120_merged_unique_cuis2genes.pkl`
- **CUI lookups**: `Pfx050120_cui_to_phens.pkl`, `Pfx050120_all_phens_to_cuis.pkl`
- **Sources**: `Pfx050120_sourced_phens.pkl` (gene, CUI) → database sources
- **Statistics**: `Pfx050120_expected_pvalue_summary.pkl` (background p-value cutoffs)

Data version identifiers: network = `pfx041520`, drug DB = `pfxDB050620`, phenotype mappings = `Pfx050120`.

### Output Files

Results are written to `results/{analysis_name}/`:

- `{drug}_assoc_table_.txt` — ranked phenotypes with p-values, BH correction, gene lists
- `{drug}_assoc_database_sources_.txt` — per-gene phenotype source databases
- `{drug}_cui_list_.txt` — top CUIs (input for clustering)
- `{drug}_merged_neighborhood.txt` — combined PPI subnetwork edges
- `{drug}_network_nodeType.txt` — node classifications (drug/target/interactor/phenotype)
- (If clustered) `*_lin_matrix_out.txt`, `cluster_membership_1.7_*.txt`, `*_dendogram_*.png`

### Enrichment Algorithm

1. Extract all proteins in the drug's network neighborhood
2. Map proteins → UMLS CUI phenotypes via `merged_genes_to_cuis.pkl`
3. Hypergeometric test against full interactome background (`pfx041520_intome_size.pkl`)
4. Benjamini-Hochberg correction
5. Filter by expected p-value background (`expected_pvalue_summary.pkl`)

### Phenotype Clustering

Requires UMLS database. The pipeline calls `umls-similarity.pl` (Perl) to compute pairwise Lin semantic similarity between CUIs, then applies Ward linkage hierarchical clustering with a fixed cutoff of `link_cutoff = 1.7` to define clusters.

## Testing

Tests use Perl modules in `tests/lib/`:

```bash
perl tests/diff.pl   # Compare result files
```

`PathFX.pm` parses result files; `PathFX/Test.pm` compares results statistically (mean/std dev of neighborhoods and association p-values).

No Python test runner is configured — example runs in `results/` serve as integration references (e.g., `FCGR2B_phenotype_analysis/`, `IL1R2_phenotype_analysis/`).

## Important Notes

- No `setup.py` / `pyproject.toml` — not a Python package; run scripts directly.
- Multiple versioned script variants exist (e.g., `_v3`, `_Pfx050120`, `_SO`). The `_Pfx050120` versions are current.
- The `_SO` suffix variants use a "stand-alone" mode that doesn't require a live UMLS connection.
- UMLS config goes in `rscs/itfc_confg.txt` (not included in repo — see `install.md`).
