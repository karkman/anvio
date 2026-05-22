#!/usr/bin/env python
"""Compare all genes in a contigs database using k-mer/Jaccard similarity."""

import os
import sys
import sqlite3
import pickle
import multiprocess as multiprocessing
import hashlib
import time

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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
                   "for genes and their flanking regions. Supports k-mers up to 13 via caching.")


# Global variables for workers
gene_data_global = None
gene_ids_global = None
min_similarity_global = 0.0


def init_worker(data, ids, min_sim):
    """Initialize worker process with shared data."""
    global gene_data_global, gene_ids_global, min_similarity_global
    gene_data_global = data
    gene_ids_global = ids
    min_similarity_global = min_sim


def jaccard_similarity_sets(set1, set2):
    """Compute Jaccard similarity between two sets."""
    inter = len(set1 & set2)
    union = len(set1 | set2)
    return inter / union if union > 0 else 0.0


def worker(i_idx, group_gene_ids=None):
    """Worker function to compute similarities for one gene against all subsequent genes in a group."""
    target_gene_ids = group_gene_ids if group_gene_ids is not None else gene_ids_global
    gid1 = target_gene_ids[i_idx]
    data1 = gene_data_global[gid1]
    results = []

    global min_similarity_global

    for j in range(i_idx + 1, len(target_gene_ids)):
        gid2 = target_gene_ids[j]
        data2 = gene_data_global[gid2]

        # MinHash Filter
        passed_filter = True
        if min_similarity_global > 0:
            if minhash_jaccard(data1['sketch'], data2['sketch']) < min_similarity_global:
                passed_filter = False

        # Only compute similarities if we passed the filter
        if passed_filter:
            # Compute similarities using pre-computed sets
            gene_sim = jaccard_similarity_sets(data1['gene_kmers'], data2['gene_kmers'])
            up_sim = jaccard_similarity_sets(data1['upstream_kmers'], data2['upstream_kmers'])
            down_sim = jaccard_similarity_sets(data1['downstream_kmers'], data2['downstream_kmers'])
            comb_sim = jaccard_similarity_sets(data1['combined_kmers'], data2['combined_kmers'])
            results.append((gid1, gid2, gene_sim, up_sim, down_sim, comb_sim))
        else:
            # Still return a result for filtering tracking, but with None similarities
            results.append((gid1, gid2, None, None, None, None))

    return results


import hashlib

def get_kmers(seq, k):
    """Extract k-mer set from a sequence."""
    if len(seq) < k:
        return set()
    return {int(hashlib.md5(seq[i:i+k].encode()).hexdigest(), 16)
            for i in range(len(seq) - k + 1)}


def get_minhash_sketch(kmers, num_hashes=100):
    """Create a MinHash sketch from a set of k-mers."""
    hashes = []
    for kmer in kmers:
        h = int(hashlib.md5(str(kmer).encode()).hexdigest(), 16)
        hashes.append(h)
    hashes.sort()
    return set(hashes[:num_hashes])


def minhash_jaccard(sketch1, sketch2):
    """Estimate Jaccard similarity using MinHash sketches."""
    inter = len(sketch1 & sketch2)
    union = len(sketch1 | sketch2)
    return inter / union if union > 0 else 0.0


def main():
    """Main function for anvi-compare-genes."""
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
    run_all_against_all = A('run_all_against_all')

    min_gene_length = A('min_gene_length') or 0

    try:
        if not contigs_db_path:
            raise ConfigError("You must provide a contigs database.")

        # Use ContigsSuperclass for high-level access
        c = ContigsSuperclass(args)

        if args.list_annotation_sources:
            annotation_sources = c.a_meta.get('gene_function_sources')
            if not annotation_sources:
                raise ConfigError("This contigs database does not have any functional annotations :/")

            run.warning('', 'FUNCTIONAL ANNOTATION SOURCE%s FOUND' % ('S' if len(annotation_sources) > 1 else ''), lc='yellow')
            for annotation_source_name in annotation_sources:
                if annotation_source_name == annotation_sources[-1]:
                    run.info_single('%s' % annotation_source_name, nl_after=1)
                else:
                    run.info_single('%s' % annotation_source_name)
            sys.exit(0)

        if not output_path:
            raise ConfigError("You must provide an output file path (unless you are using --list-annotation-sources).")

        # Validate that a comparison mode is selected
        if not any([annotation_source, args.gene_caller_ids, run_all_against_all]):
            raise ConfigError("You must specify a comparison mode. Please use either '--compare-by-annotation-source', "
                              "'--gene-caller-ids', or '--run-all-against-all'.")

        if run_all_against_all and not (annotation_source or args.gene_caller_ids):
            run.warning("You have selected '--run-all-against-all'. For large metagenomes, this may be extremely "
                        "slow and memory-heavy. Consider using functional grouping instead.")

        filesnpaths.is_output_file_writable(output_path)
        c.init_functions()

        # Get gene IDs
        if args.gene_caller_ids:
            if not os.path.exists(args.gene_caller_ids):
                raise ConfigError(f"Gene caller IDs file not found: {args.gene_caller_ids}")
            with open(args.gene_caller_ids, 'r') as f:
                gene_ids = sorted([int(line.strip()) for line in f if line.strip()])
        else:
            gene_ids = sorted(list(c.genes_in_contigs_dict.keys()))

        if not gene_ids:
            raise ConfigError("No genes found in the contigs database or provided list.")

        # Filter by length
        if min_gene_length > 0:
            gene_ids = [gid for gid in gene_ids if (c.genes_in_contigs_dict[gid]['stop'] - c.genes_in_contigs_dict[gid]['start']) >= min_gene_length]
            if not gene_ids:
                raise ConfigError(f"No genes left after filtering with --min-gene-length {min_gene_length}.")
            run.info("Genes remaining after length filter", len(gene_ids))

        # Initialize clustering graph if requested
        clusters_graph = None
        if args.cluster_results:
            import networkx as nx
            clusters_graph = nx.Graph()
            clusters_graph.add_nodes_from(gene_ids)

        # Group genes
        gene_groups = [] # List of tuples: (group_name, list_of_ids)
        if annotation_source:
            # Support comma-separated annotation sources
            annotation_sources = [s.strip() for s in annotation_source.split(',')]
            run.info("Grouping genes by sources", ', '.join(annotation_sources))
            groups = {}
            skipped_genes = 0
            for gid in gene_ids:
                fn_name = None
                # Check each annotation source in order until we find one with an annotation
                for src in annotation_sources:
                    if gid in c.gene_function_calls_dict:
                        f = c.gene_function_calls_dict[gid].get(src)
                        if f and f[1]:
                            fn_name = f[1]
                            break  # Found an annotation, use this source

                # Skip genes that don't have any of the requested annotations
                if fn_name is None:
                    skipped_genes += 1
                    continue

                if fn_name not in groups:
                    groups[fn_name] = []
                groups[fn_name].append(gid)

            if skipped_genes > 0:
                run.info_single(f"Skipped {skipped_genes} genes without any of the requested annotations")

            gene_groups = [(name, sorted(g)) for name, g in groups.items() if len(g) > 1]
            total_pairs = sum(len(g) * (len(g) - 1) // 2 for _, g in gene_groups)
            run.info("Total pairs to compare", f"{total_pairs:,} (grouped into {len(gene_groups)} functional categories)")
        else:
            gene_groups = [("All", gene_ids)]
            total_pairs = len(gene_ids) * (len(gene_ids) - 1) // 2

        # Handle Caching
        gene_data = {}
        cached_gids = set()
        if cache_file:
            run.info("Cache file", cache_file)
            db_exists = os.path.exists(cache_file)
            conn = sqlite3.connect(cache_file)
            cursor = conn.cursor()
            if db_exists:
                cursor.execute("SELECT value FROM meta WHERE key='kmer_size'")
                row = cursor.fetchone()
                if row and int(row[0]) != kmer_size:
                    conn.close()
                    raise ConfigError(f"Cache file {cache_file} was created with k={row[0]}, but you requested k={kmer_size}.")
                
                cursor.execute("SELECT gene_callers_id, data FROM kmers")
                for gid, blob in cursor.fetchall():
                    if gid in gene_ids:
                        gene_data[gid] = pickle.loads(blob)
                        cached_gids.add(gid)
            else:
                cursor.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
                cursor.execute("CREATE TABLE kmers (gene_callers_id INTEGER PRIMARY KEY, data BLOB)")
                cursor.execute("INSERT INTO meta VALUES ('kmer_size', ?)", (str(kmer_size),))
                conn.commit()

        gids_to_compute = [gid for gid in gene_ids if gid not in cached_gids]

        # Show progress for k-mer computation
        if gids_to_compute:
            progress.new('Computing k-mers for genes', progress_total_items=len(gids_to_compute))

            _, gene_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gids_to_compute,
                                                                flank_length=0)
            _, full_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gids_to_compute,
                                                                flank_length=flank_length)

            new_data_to_cache = []
            for i, gid in enumerate(gids_to_compute):
                g_seq = gene_seqs[gid]['sequence']
                f_seq = full_seqs[gid]['sequence']

                if flank_length > 0:
                    gene_len = len(g_seq)
                    idx = f_seq.find(g_seq)
                    if idx == -1:
                        run.warning(f"Gene sequence not found in flanking sequence for gene {gid}. Using empty flanking regions.")
                        upstream_seq = ""
                        downstream_seq = ""
                    else:
                        upstream_seq = f_seq[:idx]
                        downstream_seq = f_seq[idx + gene_len:]
                else:
                    upstream_seq = ""
                    downstream_seq = ""

                kmers = get_kmers(g_seq, kmer_size)
                data = {
                    'gene_kmers': kmers,
                    'upstream_kmers': get_kmers(upstream_seq, kmer_size),
                    'downstream_kmers': get_kmers(downstream_seq, kmer_size),
                    'combined_kmers': get_kmers(upstream_seq + downstream_seq, kmer_size),
                    'sketch': get_minhash_sketch(kmers)
                }
                gene_data[gid] = data
                if cache_file:
                    new_data_to_cache.append((gid, pickle.dumps(data)))

                # Update progress
                progress.increment()
                if (i + 1) % 500 == 0 or (i + 1) == len(gids_to_compute):
                    progress.update(f'Processed {i + 1}/{len(gids_to_compute)} genes')

            progress.end()

            if cache_file and new_data_to_cache:
                cursor.executemany("INSERT INTO kmers VALUES (?, ?)", new_data_to_cache)
                conn.commit()
        
        if cache_file:
            conn.close()

        # Write results
        with open(output_path, 'w') as outf:
            outf.write("gene_callers_id_1\tgene_callers_id_2\tgene_similarity\tupstream_similarity\tdownstream_similarity\tcombined_flank_similarity\tannotation_1\tannotation_2\n")

            progress.new('Computing gene similarities', progress_total_items=total_pairs)

            pool = multiprocessing.Pool(processes=num_threads, initializer=init_worker, initargs=(gene_data, gene_ids, args.min_similarity))

            count = 0
            start_time = time.time()
            last_update_time = start_time
            UPDATE_INTERVAL_SECONDS = 5.0  # Update every 5 seconds

            for name, group in gene_groups:
                num_genes_in_group = len(group)
                from functools import partial
                worker_with_group = partial(worker, group_gene_ids=group)

                for results in pool.imap_unordered(worker_with_group, range(num_genes_in_group)):
                    for gid1, gid2, g_sim, u_sim, d_sim, c_sim in results:
                        count += 1

                        # Only write out valid comparisons (where similarity is not None)
                        if g_sim is not None:
                            if annotation_source:
                                # For annotation source mode, both genes in a pair have the same annotation (the group name)
                                outf.write(f"{gid1}\t{gid2}\t{g_sim:.6f}\t{u_sim:.6f}\t{d_sim:.6f}\t{c_sim:.6f}\t{name}\t{name}\n")
                            else:
                                # Get annotations
                                ann1_list = []
                                if gid1 in c.gene_function_calls_dict:
                                    for s, f in c.gene_function_calls_dict[gid1].items():
                                        if f and f[1]:
                                            ann1_list.append(f"{s}:{f[1]}")
                                ann1 = "; ".join(ann1_list)

                                ann2_list = []
                                if gid2 in c.gene_function_calls_dict:
                                    for s, f in c.gene_function_calls_dict[gid2].items():
                                        if f and f[1]:
                                            ann2_list.append(f"{s}:{f[1]}")
                                ann2 = "; ".join(ann2_list)

                                outf.write(f"{gid1}\t{gid2}\t{g_sim:.6f}\t{u_sim:.6f}\t{d_sim:.6f}\t{c_sim:.6f}\t{ann1}\t{ann2}\n")

                            # Add to clustering graph
                            if clusters_graph is not None and g_sim >= args.clustering_similarity_threshold:
                                clusters_graph.add_edge(gid1, gid2)

                        # Update progress based on time interval or every 10000 pairs (reduced frequency to minimize flickering)
                        current_time = time.time()
                        if count % 10000 == 0 or (current_time - last_update_time) >= UPDATE_INTERVAL_SECONDS:
                            elapsed = current_time - start_time
                            rate = count / elapsed if elapsed > 0 else 0
                            percent = (count / total_pairs) * 100 if total_pairs > 0 else 0
                            progress.increment(increment_to=count)
                            progress.update(
                                f'Compared {count:>12,} pairs ({percent:5.1f}%) - '
                                f'{rate:>6.0f} pairs/sec'
                            )
                            last_update_time = current_time

            progress.end()
            pool.close()
            pool.join()

        # Perform clustering if requested
        if clusters_graph is not None:
            import networkx as nx
            clusters = list(nx.connected_components(clusters_graph))
            cluster_output = os.path.splitext(output_path)[0] + ".clusters.txt"
            with open(cluster_output, 'w') as f:
                f.write("gene_caller_id\tcluster_id\n")
                for i, cluster in enumerate(clusters):
                    for gid in sorted(list(cluster)):
                        f.write(f"{gid}\t{i}\n")
            run.info("Clustering completed", f"Results written to {cluster_output}")

        run.info('Comparison completed', f'Results written to {output_path}')

    except ConfigError as e:
        progress.end()
        run.warning(str(e))
        sys.exit(-1)
    except FilesNPathsError as e:
        progress.end()
        run.warning(str(e))
        sys.exit(-1)
    except Exception as e:
        progress.end()
        if anvio.DEBUG:
            raise
        else:
            run.warning(f"An unexpected error occurred: {e}")
            sys.exit(-1)


def get_args():
    """Define arguments for the program."""
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': False, 'help': 'Output file path.'}))
    parser.add_argument(*anvio.A('kmer-size'), **anvio.K('kmer-size', {'type': int, 'default': 3}))
    parser.add_argument(*anvio.A('flank-length'), **anvio.K('flank-length', {'type': int, 'default': 500}))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads', {'type': int, 'default': 1}))
    parser.add_argument('--min-gene-length', type=int, default=0, help='Minimum length of a gene call to be included in the comparison.')
    parser.add_argument('--gene-caller-ids', help='Path to a file containing a list of gene caller IDs (one per line).')
    parser.add_argument('--cache-file', help='Optional SQLite file to cache k-mer sets.')
    parser.add_argument('--min-similarity', type=float, default=0.0, help='MinHash Jaccard similarity threshold for filtering gene pairs. Recommended: 0.1-0.3.')
    parser.add_argument('--compare-by-annotation-source', help='Optionally compare only genes that have the same annotation in this source (e.g. COG20_FUNCTION).')
    parser.add_argument('--list-annotation-sources', action='store_true', default=False, help='List available functional annotation sources and quit.')
    parser.add_argument('--run-all-against-all', action='store_true', default=False, help='Explicitly run an all-against-all comparison (NOT RECOMMENDED for large metagenomes).')
    parser.add_argument('--cluster-results', action='store_true', default=False, help='Cluster highly similar genes based on connected components.')
    parser.add_argument('--clustering-similarity-threshold', type=float, default=0.98, help='Jaccard similarity threshold for clustering. Default: 0.98.')
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
