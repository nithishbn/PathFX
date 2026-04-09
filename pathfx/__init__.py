"""
PathFX as an importable library.

Usage:
    from pathfx import PathFX

    pfx = PathFX('/path/to/PathFX')
    sig_assoc, node_list = pfx.run_network(
        targets=['TP53', 'TOP2A'],
        drug_name='DB00997__TP53',
        outdir='/path/to/results/crispr_networks/DB00997__TP53',
    )

The PathFX instance is designed to be created once and reused across many
run_network() calls — all pickle data is loaded lazily on the first call
and cached for subsequent calls.
"""

import pickle
from collections import defaultdict
from pathlib import Path

from scipy.stats import hypergeom


class PathFX:
    def __init__(self, pathfx_dir):
        """
        Parameters
        ----------
        pathfx_dir : str or Path
            Root directory of the PathFX repository.
        """
        self.pathfx_dir = Path(pathfx_dir)
        self.rscs_dir = self.pathfx_dir / 'rscs'
        self.scripts_dir = self.pathfx_dir / 'scripts'
        self._data = None

    # ── Data loading ─────────────────────────────────────────────────────────

    def _load(self):
        if self._data is not None:
            return
        r = self.rscs_dir
        print('[pathfx] loading resource files...')
        self._data = {
            'node_to_hash':  pickle.load(open(r / 'pfx041520_0.8_node_to_hashID.pkl', 'rb')),
            'nbhd_hash':     pickle.load(open(r / 'pfx041520_0.82_spec_nbhd_hash.pkl', 'rb')),
            'unique_nodes':  set(pickle.load(open(r / 'pfx041520_unique_nodes.pkl', 'rb'))),
            'intome_size':   pickle.load(open(r / 'pfx041520_intome_size.pkl', 'rb')),
            'genes_to_cuis': pickle.load(open(r / 'Pfx050120_merged_genes_to_cuis.pkl', 'rb')),
            'cui_to_phens':  pickle.load(open(r / 'Pfx050120_cui_to_phens.pkl', 'rb')),
            'cui_to_genes':  pickle.load(open(r / 'Pfx050120_merged_unique_cuis2genes.pkl', 'rb')),
            'sourced_phens': pickle.load(open(r / 'Pfx050120_sourced_phens.pkl', 'rb')),
            'expect_pvals':  pickle.load(open(
                self.pathfx_dir / 'results' / 'Pfx050120random_networks' /
                'Pfx050120_expected_pvalue_summary.pkl', 'rb'
            )),
        }
        print('[pathfx] ready.')

    # ── PPI network extraction ────────────────────────────────────────────────

    def _get_neighborhood(self, gene):
        """Return the precomputed path→score dict for one gene, or None."""
        d = self._data
        if gene not in d['node_to_hash']:
            return None
        hash_id = d['node_to_hash'][gene]
        rel_path = d['nbhd_hash']['spn' + hash_id]
        abs_path = (self.scripts_dir / rel_path).resolve()
        if not abs_path.exists():
            return None
        return pickle.load(open(abs_path, 'rb'))

    def _build_merged_neighborhood(self, targets):
        """
        Merge PPI path dicts for all targets.
        Returns (merged_dict, skipped_genes).
        When the same path appears in multiple targets, the highest score wins.
        """
        merged = {}
        skipped = []
        for t in targets:
            nbhd = self._get_neighborhood(t)
            if nbhd is None:
                skipped.append(t)
                continue
            for path, score in nbhd.items():
                if score > merged.get(path, 0):
                    merged[path] = score
        return merged, skipped

    @staticmethod
    def _paths_to_edges(path_dict):
        """
        Convert path→score dict to (nodeA, nodeB)→max_score edge dict.
        Path format: 'GENE1@GENE2@GENE3' — only the last edge is kept.
        Single nodes (no '@') are stored as (node, '').
        """
        edges = defaultdict(float)
        for path, score in path_dict.items():
            if '@' in path:
                parts = path.split('@')
                a, b = parts[-2], parts[-1]
            else:
                a, b = path, ''
            if score > edges[(a, b)]:
                edges[(a, b)] = score
        return dict(edges)

    @staticmethod
    def _node_list_from_edges(edges):
        nodes = set()
        for a, b in edges:
            nodes.add(a)
            if b:
                nodes.add(b)
        nodes.discard('')
        return list(nodes)

    # ── Phenotype enrichment ──────────────────────────────────────────────────

    def _get_assoc(self, node_list):
        genes_to_cuis = self._data['genes_to_cuis']
        assoc_count = defaultdict(int)
        assoc_genes = defaultdict(list)
        for n in node_list:
            if n in genes_to_cuis:
                for cui in genes_to_cuis[n]:
                    assoc_count[cui] += 1
                    assoc_genes[cui].append(n)
        return assoc_count, assoc_genes

    def _calc_hyp(self, node_list, Q=0.001):
        """Hypergeometric test + Benjamini-Hochberg correction, matching original logic."""
        d = self._data
        N = d['intome_size']
        cui_to_genes = d['cui_to_genes']
        cui_to_phens = d['cui_to_phens']
        n = len(node_list)

        assoc_count, assoc_genes = self._get_assoc(node_list)

        # Compute p-values
        assoc_analy = []
        for cui, k in assoc_count.items():
            K = len(cui_to_genes[cui])
            p = 1 - hypergeom.cdf(k, N, K, n)
            assoc_analy.append([cui, k, K, p])

        assoc_analy.sort(key=lambda x: (x[3], x[0]))
        m = len(assoc_analy)

        # BH correction and filtering — matches original exactly:
        # keep rows where p < BH threshold AND intome_count > 24,
        # stop early once p > BH (list is sorted by p)
        sig_assoc = []
        for i, (cui, k, K, p) in enumerate(assoc_analy):
            BH = (float(i + 1) / m) * Q
            if p < BH and K > 24:
                phen_term = cui_to_phens[cui][0]
                gene_str = ','.join(sorted(assoc_genes[cui]))
                sig_assoc.append([i + 1, phen_term, cui, k, K, p, BH, gene_str])
            elif p > BH:
                break

        return sig_assoc

    def _background_check(self, num_targets, sig_assoc):
        """Filter enriched phenotypes against pre-computed background p-values."""
        expect_pvals = self._data['expect_pvals']
        bg = expect_pvals.get(num_targets, {})
        return [
            row for row in sig_assoc
            if row[6] < bg.get(row[2], 1.0)  # BH < expected_pval for this CUI
        ]

    # ── File writing ──────────────────────────────────────────────────────────

    @staticmethod
    def _write_edges(edges, path):
        with open(path, 'w') as f:
            for (a, b), score in edges.items():
                f.write(f'{a}\t{b}\t{score}\n')

    @staticmethod
    def _write_assoc_table(sig_assoc, path):
        header = 'rank\tphenotype\tcui\tassoc in neigh\tassoc in intom\tprobability\tBenjamini-Hochberg\tgenes\n'
        with open(path, 'w') as f:
            f.write(header)
            for row in sig_assoc:
                f.write('\t'.join(str(x) for x in row) + '\n')

    def _write_sources(self, sig_assoc, path):
        sourced_phens = self._data['sourced_phens']
        with open(path, 'w') as f:
            f.write('Gene\tCUI\tSource Databases\n')
            for _rank, phen_term, cui, _k, _K, _p, _BH, gene_str in sig_assoc:
                for g in gene_str.split(','):
                    sources = ','.join(sorted(sourced_phens.get((g, cui), [])))
                    f.write(f'{g}\t{cui}\t{phen_term}\t{sources}\n')

    @staticmethod
    def _write_cui_list(sig_assoc, path, limit=50):
        cuis = [row[2] for row in (sig_assoc[:limit] if len(sig_assoc) >= 100 else sig_assoc)]
        with open(path, 'w') as f:
            f.write('\n'.join(cuis))

    def _write_viz_files(self, drug_name, targets, sig_assoc, edges, outdir):
        outdir = Path(outdir)

        # Full network: drug → targets, then all PPI edges, then phenotype → gene edges
        net_path = outdir / f'{drug_name}_merged_neighborhood__withDrugTargsAndPhens.txt'
        with open(net_path, 'w') as f:
            f.write('node1\tnode2\tedge_score\n')
            for t in targets:
                f.write(f'{drug_name}\t{t}\t1.0\n')
            for (a, b), score in edges.items():
                if b:
                    f.write(f'{a}\t{b}\t{score}\n')
            for _rank, phen_term, _cui, _k, _K, _p, _BH, gene_str in sig_assoc:
                for g in gene_str.split(','):
                    f.write(f'{phen_term}\t{g}\t1.0\n')

        # Node type annotations
        nt_path = outdir / f'{drug_name}_network_nodeType.txt'
        with open(nt_path, 'w') as f:
            f.write('node_name\tnode_type\n')
            f.write(f'{drug_name}\tdrug\n')
            for t in targets:
                f.write(f'{t}\tdrug_target\n')
            for _rank, phen_term, _cui, _k, _K, _p, _BH, _gene_str in sig_assoc:
                f.write(f'{phen_term}\tphenotype\n')

    # ── Public API ────────────────────────────────────────────────────────────

    def run_network(self, targets, drug_name, outdir):
        """
        Run PathFX for a list of protein targets and write results to outdir.

        Parameters
        ----------
        targets   : list of str — HUGO gene symbols (drug targets + any extras)
        drug_name : str — label for this run, used in output filenames
        outdir    : str or Path — directory to write results into (created if absent)

        Returns
        -------
        sig_assoc : list of enriched phenotype rows (empty list if none found)
                    Each row: [rank, phenotype, cui, k, K, p, BH, gene_str]
        node_list : list of str — all proteins in the merged PPI neighborhood
        """
        self._load()
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Filter to genes present in the interactome
        valid = [t for t in targets if t in self._data['unique_nodes']]
        if not valid:
            print(f'[pathfx] {drug_name}: no targets found in interactome, skipping')
            return [], []

        # Build merged PPI neighborhood
        path_dict, skipped = self._build_merged_neighborhood(valid)
        if skipped:
            print(f'[pathfx] {drug_name}: skipped (no precomputed neighborhood): {skipped}')

        edges = self._paths_to_edges(path_dict)
        node_list = self._node_list_from_edges(edges)

        # Write merged network file
        self._write_edges(edges, outdir / f'{drug_name}_merged_neighborhood_.txt')

        # Phenotype enrichment
        sig_assoc = self._calc_hyp(node_list)
        sig_assoc = self._background_check(len(valid), sig_assoc)

        # Write output files (matching original naming conventions)
        self._write_assoc_table(sig_assoc, outdir / f'{drug_name}_merged_neighborhood__assoc_table_.txt')
        self._write_sources(sig_assoc, outdir / f'{drug_name}_merged_neighborhood__assoc_database_sources_.txt')
        self._write_cui_list(sig_assoc, outdir / f'{drug_name}_merged_neighborhood__cui_list_.txt')
        self._write_viz_files(drug_name, valid, sig_assoc, edges, outdir)

        return sig_assoc, node_list
