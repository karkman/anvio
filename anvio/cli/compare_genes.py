#!/usr/bin/env python
"""Compare all genes in a contigs database using k-mer/Jaccard similarity."""

import os
import sys
import sqlite3
import pickle
import multiprocess as multiprocessing

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

def worker(i_idx):
    """Worker function to compute similarities for one gene against all subsequent genes."""
    gid1 = gene_ids_global[i_idx]
    data1 = gene_data_global[gid1]
    results = []

    # Get threshold from shared state if possible, or assume 0
    # For simplicity, pass it as part of init_worker if needed.
    # Here we'll assume a global threshold 'min_similarity_global'
    global min_similarity_global

    for j in range(i_idx + 1, len(gene_ids_global)):
        gid2 = gene_ids_global[j]
        data2 = gene_data_global[gid2]

        # MinHash Filter
        if min_similarity_global > 0:
            if minhash_jaccard(data1['sketch'], data2['sketch']) < min_similarity_global:
                continue

        # Compute similarities using pre-computed sets
        gene_sim = jaccard_similarity_sets(data1['gene_kmers'], data2['gene_kmers'])
        up_sim = jaccard_similarity_sets(data1['upstream_kmers'], data2['upstream_kmers'])
        down_sim = jaccard_similarity_sets(data1['downstream_kmers'], data2['downstream_kmers'])
        comb_sim = jaccard_similarity_sets(data1['combined_kmers'], data2['combined_kmers'])

        results.append((gid1, gid2, gene_sim, up_sim, down_sim, comb_sim))

    return results



import hashlib

def get_kmers(seq, k):
    """Extract k-mer set from a sequence."""
    if len(seq) < k:
        return set()
    return {hash(seq[i:i+k]) for i in range(len(seq) - k + 1)}


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

    try:
        if not contigs_db_path:
            raise ConfigError("You must provide a contigs database.")
        if not output_path:
            raise ConfigError("You must provide an output file path.")

        filesnpaths.is_output_file_writable(output_path)

        # Use ContigsSuperclass for high-level access
        c = ContigsSuperclass(args)
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
                    gene_data[gid] = pickle.loads(blob)
                    cached_gids.add(gid)
            else:
                cursor.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
                cursor.execute("CREATE TABLE kmers (gene_callers_id INTEGER PRIMARY KEY, data BLOB)")
                cursor.execute("INSERT INTO meta VALUES ('kmer_size', ?)", (str(kmer_size),))
                conn.commit()

        gids_to_compute = [gid for gid in gene_ids if gid not in cached_gids]
        
        if gids_to_compute:
            _, gene_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gids_to_compute,
                                                                flank_length=0)
            _, full_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gids_to_compute,
                                                                flank_length=flank_length)

            new_data_to_cache = []
            for gid in gids_to_compute:
                g_seq = gene_seqs[gid]['sequence']
                f_seq = full_seqs[gid]['sequence']

                if flank_length > 0:
                    gene_len = len(g_seq)
                    idx = f_seq.find(g_seq)
                    upstream_seq = f_seq[:idx] if idx != -1 else ""
                    downstream_seq = f_seq[idx + gene_len:] if idx != -1 else ""
                else:
                    upstream_seq = ""
                    downstream_seq = ""

                data = {
                    'gene_kmers': get_kmers(g_seq, kmer_size),
                    'upstream_kmers': get_kmers(upstream_seq, kmer_size),
                    'downstream_kmers': get_kmers(downstream_seq, kmer_size),
                    'combined_kmers': get_kmers(upstream_seq + downstream_seq, kmer_size),
                    'sketch': get_minhash_sketch(get_kmers(g_seq, kmer_size))
                }
                gene_data[gid] = data
                if cache_file:
                    new_data_to_cache.append((gid, pickle.dumps(data)))
            
            if cache_file and new_data_to_cache:
                cursor.executemany("INSERT INTO kmers VALUES (?, ?)", new_data_to_cache)
                conn.commit()
        
        if cache_file:
            conn.close()

        # Write results
        with open(output_path, 'w') as outf:
            outf.write("gene_callers_id_1\tgene_callers_id_2\tgene_similarity\tupstream_similarity\tdownstream_similarity\tcombined_flank_similarity\tannotations_1\tannotations_2\n")

            total_genes = len(gene_ids)
            total_pairs = total_genes * (total_genes - 1) // 2
            progress.new('Computing gene similarities', progress_total_items=total_pairs)

            pool = multiprocessing.Pool(processes=num_threads, initializer=init_worker, initargs=(gene_data, gene_ids, args.min_similarity))
            
            count = 0
            # We use imap_unordered for better memory efficiency and responsiveness
            for results in pool.imap_unordered(worker, range(total_genes)):
                for gid1, gid2, g_sim, u_sim, d_sim, c_sim in results:
                    # Get annotations
                    ann1 = "; ".join([f"{s}:{f[0]}" for s, f in c.gene_function_calls_dict.get(gid1, {}).items() if f])
                    ann2 = "; ".join([f"{s}:{f[0]}" for s, f in c.gene_function_calls_dict.get(gid2, {}).items() if f])

                    outf.write(f"{gid1}\t{gid2}\t{g_sim:.6f}\t{u_sim:.6f}\t{d_sim:.6f}\t{c_sim:.6f}\t{ann1}\t{ann2}\n")

                    count += 1
                    progress.increment()
                    if count % 1000 == 0:
                        progress.update(f'Compared {count} pairs')

            progress.end()
            pool.close()
            pool.join()


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
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True, 'help': 'Output file path.'}))
    parser.add_argument(*anvio.A('kmer-size'), **anvio.K('kmer-size', {'type': int, 'default': 3}))
    parser.add_argument(*anvio.A('flank-length'), **anvio.K('flank-length', {'type': int, 'default': 500}))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads', {'type': int, 'default': 1}))
    parser.add_argument('--gene-caller-ids', help='Path to a file containing a list of gene caller IDs (one per line).')
    parser.add_argument('--cache-file', help='Optional SQLite file to cache k-mer sets.')
    parser.add_argument('--min-similarity', type=float, default=0.0, help='MinHash Jaccard similarity threshold for filtering gene pairs. Recommended: 0.1-0.3.')
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
