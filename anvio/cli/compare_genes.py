#!/usr/bin/env python
"""Compare all genes in a contigs database using k-mer/Jaccard similarity."""

import os
import sys
import sqlite3
import pickle
import multiprocess as multiprocessing
import hashlib
import time
import struct
from functools import partial

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.utils as utils

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import ContigsSuperclass


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['karkman']
__requires__ = ['contigs-db']
__provides__ = []
__description__ = ("Compares all genes in a contigs database using k-mer/Jaccard similarity "
                   "for genes and their flanking regions. Supports k-mers up to 13 via caching "
                   "and probabilistic pre-filtering using MinHash sketches.")


# Global variables for workers
gene_data_global = None
min_similarity_global = 0.0
kmer_size_global = 3


def init_worker(data, min_sim, k_size):
    """Initialize worker process with shared data."""
    global gene_data_global, min_similarity_global, kmer_size_global
    gene_data_global = data
    min_similarity_global = min_sim
    kmer_size_global = k_size


def get_kmers_packed(seq, k):
    """Extract k-mer set from a sequence using fast 2-bit packing."""
    if not seq or len(seq) < k:
        return frozenset()

    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    kmers = set()
    current_pack = 0
    mask = (1 << (2 * k)) - 1
    valid_len = 0

    for base in seq:
        if base in base_map:
            current_pack = ((current_pack << 2) | base_map[base]) & mask
            valid_len += 1
            if valid_len >= k:
                kmers.add(current_pack)
        else:
            valid_len = 0
            current_pack = 0

    return frozenset(kmers)


def jaccard_similarity_sets(set1, set2):
    """Compute Jaccard similarity between two sets."""
    inter = len(set1 & set2)
    union = len(set1 | set2)
    return inter / union if union > 0 else 0.0


def get_minhash_sketch(kmers, num_hashes=100, seed=1):
    """Compute a MinHash sketch using multiple hash functions.

    Uses hashlib.md5 with different seeds to generate independent hash functions.
    Returns a frozenset of the smallest num_hashes hash values.
    """
    if not kmers:
        return frozenset()

    hashes = []
    for seed_val in range(1, num_hashes + 1):
        h = hashlib.md5()
        h.update(struct.pack('I', seed_val))
        for kmer in kmers:
            h.update(struct.pack('I', kmer))
        hashes.append(h.digest())

    # Use the first num_hashes distinct hash values
    # Sort by the hash bytes to get a consistent ordering
    sorted_hashes = sorted(hashes)[:num_hashes]
    return frozenset(sorted_hashes)


def minhash_jaccard(sketch1, sketch2):
    """Estimate Jaccard similarity from two MinHash sketches."""
    inter = len(sketch1 & sketch2)
    union = len(sketch1 | sketch2)
    return inter / union if union > 0 else 0.0


def worker(pairs):
    """Worker function to compute similarities for a given list of gene pairs.

    Each pair is a tuple (gid1, gid2, ann1, ann2).
    Returns a list of result tuples.
    """
    results = []
    global gene_data_global, min_similarity_global, kmer_size_global

    for gid1, gid2, ann1, ann2 in pairs:
        data1 = gene_data_global[gid1]
        data2 = gene_data_global[gid2]

        passed_filter = True
        if min_similarity_global > 0:
            if minhash_jaccard(data1['sketch'], data2['sketch']) < min_similarity_global:
                passed_filter = False

        if passed_filter:
            gene_sim = jaccard_similarity_sets(data1['gene_kmers'], data2['gene_kmers'])
            up_sim = jaccard_similarity_sets(data1['upstream_kmers'], data2['upstream_kmers'])
            down_sim = jaccard_similarity_sets(data1['downstream_kmers'], data2['downstream_kmers'])
            comb_sim = jaccard_similarity_sets(data1['combined_kmers'], data2['combined_kmers'])
            results.append((gid1, gid2, gene_sim, up_sim, down_sim, comb_sim, ann1, ann2))
        else:
            results.append((gid1, gid2, None, None, None, None, ann1, ann2))

    return results


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    contigs_db_path = A('contigs_db')
    output_path = A('output_file')
    kmer_size = A('kmer_size') or 3
    flank_length = A('flank_length') or 500
    num_threads = A('num_threads') or 1
    cache_file = A('cache_file')
    annotation_source = A('compare_by_annotation_source')

    try:
        if not contigs_db_path:
            raise ConfigError("You must provide a contigs database.")
        c = ContigsSuperclass(args)

        # Validate annotation sources early
        available_sources = c.a_meta.get('gene_function_sources', [])
        if args.list_annotation_sources:
            run.warning('', 'ANNOTATION SOURCES FOUND', lc='yellow')
            for s in available_sources: run.info_single(s)
            sys.exit(0)

        if annotation_source:
            requested_sources = [s.strip() for s in annotation_source.split(',')]
            for s in requested_sources:
                if s not in available_sources:
                    raise ConfigError(f"Annotation source '{s}' not found in this database. "
                                      f"Available sources are: {', '.join(available_sources)}")

        if not output_path: raise ConfigError("Missing output file path.")

        # Require at least one comparison mode
        if not any([annotation_source, args.gene_caller_ids]):
            raise ConfigError("You must specify a comparison mode. Use --compare-by-annotation-source "
                              "to group genes by functional annotation, or --gene-caller-ids to "
                              "compare a specific set of genes. For a full all-against-all comparison, "
                              "provide a list of all gene caller IDs via --gene-caller-ids.")

        filesnpaths.is_output_file_writable(output_path)
        c.init_functions()

        if args.gene_caller_ids:
            with open(args.gene_caller_ids, 'r') as f:
                gene_ids = sorted([int(line.strip()) for line in f if line.strip()])
        else:
            gene_ids = sorted(list(c.genes_in_contigs_dict.keys()))

        if not gene_ids:
            raise ConfigError("No genes found in the contigs database or provided list.")

        # --- STEP 1: INITIAL GROUPING ---
        gene_groups = []
        gene_annotations = {}  # Map gene_id -> annotation string
        # Standard anvi'o strings for hypothetical proteins
        hypothetical_terms = ["hypothetical", "hypothetical protein", "conserved hypothetical",
                              "conserved hypotheticals", "Conserved hypothetical protein"]

        if annotation_source:
            sources = [s.strip() for s in annotation_source.split(',')]
            run.info("Grouping genes by sources", ', '.join(sources))
            groups = {}
            skipped_genes = 0
            for gid in gene_ids:
                fn = None
                for s in sources:
                    res = c.gene_function_calls_dict.get(gid, {}).get(s)
                    if res and res[1]:
                        fn = res[1]
                        break

                # Logic: Skip genes that have NO annotation hit in the requested source(s)
                if fn is None:
                    skipped_genes += 1
                    continue

                # Store the original annotation for this gene
                gene_annotations[gid] = fn

                # If the function is any kind of hypothetical, pool them all into one group
                if any(term.lower() in fn.lower() for term in hypothetical_terms):
                    fn = "Hypothetical"

                if fn not in groups:
                    groups[fn] = []
                groups[fn].append(gid)

            if skipped_genes > 0:
                run.info_single(f"Skipped {skipped_genes} genes without any annotation hits in the requested sources", level=1)

            gene_groups = [(n, sorted(g)) for n, g in groups.items() if len(g) > 1]
            total_possible_pairs = sum(len(g) * (len(g) - 1) // 2 for _, g in gene_groups)
            run.info("Functional categories identified", f"{len(gene_groups)} (covering {total_possible_pairs:,} total potential pairs)")
        else:
            # No annotation source: compare the provided gene IDs as one group, but still fetch annotations from all available sources
            gene_groups = [("All", gene_ids)]
            total_possible_pairs = len(gene_ids) * (len(gene_ids) - 1) // 2

            # Fetch annotations from ALL available sources and combine them
            available_sources = c.a_meta.get('gene_function_sources', [])
            if available_sources:
                run.info("Fetching annotations from all available sources", ', '.join(available_sources))
                for gid in gene_ids:
                    annotations = []
                    for source in available_sources:
                        res = c.gene_function_calls_dict.get(gid, {}).get(source)
                        if res and res[1]:
                            annotations.append(f"{source}: {res[1]}")
                    gene_annotations[gid] = " | ".join(annotations) if annotations else ""
            else:
                gene_annotations = {gid: "" for gid in gene_ids}

        # --- STEP 2: CACHE LOOKUP ---
        gene_data = {}
        cached_gids = set()
        cache_conn = None
        cache_cursor = None
        if cache_file:
            db_exists = os.path.exists(cache_file)
            cache_conn = sqlite3.connect(cache_file)
            cache_cursor = cache_conn.cursor()
            if db_exists:
                # Validate k-mer size matches
                cache_cursor.execute("SELECT value FROM meta WHERE key = 'kmer_size'")
                row = cache_cursor.fetchone()
                if row:
                    cached_kmer_size = int(row[0])
                    if cached_kmer_size != kmer_size:
                        raise ConfigError(f"Cache file k-mer size ({cached_kmer_size}) does not match "
                                          f"current k-mer size ({kmer_size}). Please use a different cache file "
                                          f"or re-run with --kmer-size {cached_kmer_size}.")
                cache_cursor.execute("SELECT gene_callers_id, data FROM kmers")
                for gid, blob in cache_cursor.fetchall():
                    if gid in gene_ids:
                        gene_data[gid] = pickle.loads(blob)
                        cached_gids.add(gid)
            else:
                cache_cursor.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
                cache_cursor.execute("CREATE TABLE kmers (gene_callers_id INTEGER PRIMARY KEY, data BLOB)")
                cache_cursor.execute("INSERT INTO meta VALUES ('kmer_size', ?)", (str(kmer_size),))
                cache_conn.commit()

        # --- STEP 3: SELECTIVE LOADING & TOKENIZATION ---
        gene_ids_to_load = [gid for gid in gene_ids if gid not in cached_gids]
        if gene_ids_to_load:
            progress.new('Selective loading and tokenization', progress_total_items=len(gene_ids_to_load))
            chunk_size = 5000
            new_data_to_cache = []
            old_quiet, old_verbose = anvio.QUIET, c.run.verbose
            anvio.QUIET, c.run.verbose = True, False
            try:
                for i in range(0, len(gene_ids_to_load), chunk_size):
                    batch = gene_ids_to_load[i:i + chunk_size]
                    _, seqs_raw = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=batch, flank_length=0)
                    _, full_raw = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=batch, flank_length=flank_length)
                    for gid in batch:
                        g, f = seqs_raw[gid]['sequence'], full_raw[gid]['sequence']
                        idx = f.find(g)
                        up, down = (f[:idx], f[idx + len(g):]) if idx != -1 else ("", "")
                        gene_kmers = get_kmers_packed(g, kmer_size)
                        upstream_kmers = get_kmers_packed(up, kmer_size)
                        downstream_kmers = get_kmers_packed(down, kmer_size)
                        combined_kmers = get_kmers_packed(up + down, kmer_size)
                        data = {
                            'gene_kmers': gene_kmers,
                            'upstream_kmers': upstream_kmers,
                            'downstream_kmers': downstream_kmers,
                            'combined_kmers': combined_kmers,
                            'sketch': get_minhash_sketch(gene_kmers)
                        }
                        gene_data[gid] = data
                        if cache_file: new_data_to_cache.append((gid, pickle.dumps(data)))
                        progress.increment()
            finally:
                anvio.QUIET, c.run.verbose = old_quiet, old_verbose
            progress.end()
            if cache_file and new_data_to_cache:
                cache_cursor.executemany("INSERT INTO kmers VALUES (?, ?)", new_data_to_cache)
                cache_conn.commit()
        if cache_file: cache_conn.close()

        # --- STEP 4: FINAL COMPARISON & WRITING ---
        clusters_graph = None
        if args.cluster_results:
            import networkx as nx
            clusters_graph = nx.Graph(); clusters_graph.add_nodes_from(gene_ids)

        with open(output_path, 'w') as outf:
            outf.write("gene_callers_id_1\tgene_callers_id_2\tgene_similarity\tupstream_similarity\tdownstream_similarity\tcombined_flank_similarity\tannotation_1\tannotation_2\n")
            progress.new('Computing final similarities', progress_total_items=total_possible_pairs)
            pool = multiprocessing.Pool(processes=num_threads, initializer=init_worker, initargs=(gene_data, args.min_similarity, kmer_size))
            count, start_t, last_u = 0, time.time(), time.time()

            def pair_gen():
                for name, group in gene_groups:
                    pairs = []
                    for i in range(len(group)):
                        for j in range(i + 1, len(group)):
                            ann1 = gene_annotations.get(group[i], "")
                            ann2 = gene_annotations.get(group[j], "")
                            pairs.append((group[i], group[j], ann1, ann2))
                    yield pairs

            it = pool.imap_unordered(worker, pair_gen())

            for results_list in it:
                for gid1, gid2, g_s, u_s, d_s, c_s, ann1, ann2 in results_list:
                    count += 1
                    if g_s is not None:
                        outf.write(f"{gid1}\t{gid2}\t{g_s:.6f}\t{u_s:.6f}\t{d_s:.6f}\t{c_s:.6f}\t{ann1}\t{ann2}\n")
                        if clusters_graph is not None and g_s >= args.clustering_similarity_threshold: clusters_graph.add_edge(gid1, gid2)
                    if count % 1000 == 0 or time.time() - last_u >= 5.0:
                        progress.increment(increment_to=count)
                        rate = count / (time.time() - start_t) if time.time() > start_t else 0
                        progress.update(f'Compared {count:,} pairs - {rate:6.0f} pairs/sec')
                        last_u = time.time()
            progress.end()
            pool.close(); pool.join()

        if clusters_graph is not None:
            clusters = list(nx.connected_components(clusters_graph))
            c_out = os.path.splitext(output_path)[0] + ".clusters.txt"
            with open(c_out, 'w') as f:
                f.write("gene_caller_id\tcluster_id\n")
                for i, cl in enumerate(clusters):
                    for g in sorted(list(cl)): f.write(f"{g}\t{i}\n")
            run.info("Clustering completed", c_out)
        run.info('Comparison completed', output_path)

    except Exception as e:
        progress.end()
        if anvio.DEBUG: raise
        run.warning(f"Error: {e}"); sys.exit(-1)

def get_args():
    from anvio.argparse import ArgumentParser
    p = ArgumentParser(description=__description__)
    p.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    p.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': False}))
    p.add_argument(*anvio.A('kmer-size'), **anvio.K('kmer-size', {'type': int, 'default': 3}))
    p.add_argument(*anvio.A('flank-length'), **anvio.K('flank-length', {'type': int, 'default': 500}))
    p.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads', {'type': int, 'default': 1}))
    p.add_argument('--gene-caller-ids')
    p.add_argument('--cache-file')
    p.add_argument('--min-similarity', type=float, default=0.0)
    p.add_argument('--compare-by-annotation-source')
    p.add_argument('--list-annotation-sources', action='store_true')
    p.add_argument('--cluster-results', action='store_true')
    p.add_argument('--clustering-similarity-threshold', type=float, default=0.98)
    return p.get_args(p)

if __name__ == '__main__':
    main()