#!/usr/bin/env python
"""Compare all genes in a contigs database using k-mer/Jaccard similarity."""

import sys

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
                   "for genes and their flanking regions.")


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

    try:
        if not contigs_db_path:
            raise ConfigError("You must provide a contigs database.")
        if not output_path:
            raise ConfigError("You must provide an output file path.")

        filesnpaths.is_output_file_writable(output_path)

        # Use ContigsSuperclass for high-level access
        c = ContigsSuperclass(args)

        # Get gene IDs
        gene_ids = sorted(list(c.genes_in_contigs_dict.keys()))
        if not gene_ids:
            raise ConfigError("No genes found in the contigs database.")

        # Extract sequences using proper anvi'o methods
        # 1. Get gene sequences (flank_length=0)
        progress.new('Initializing gene sequences')
        _, gene_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gene_ids,
                                                            flank_length=0)

        # 2. Get sequences with flanks
        progress.new('Initializing flanking sequences')
        _, full_seqs = c.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gene_ids,
                                                            flank_length=flank_length)

        # Prepare data structures
        gene_data = {}
        for gid in gene_ids:
            g_seq = gene_seqs[gid]['sequence']
            f_seq = full_seqs[gid]['sequence']

            # If flank_length > 0, we can extract upstream/downstream
            if flank_length > 0:
                gene_len = len(g_seq)

                # The gene starts at some index 'idx' in f_seq.
                idx = f_seq.find(g_seq)
                if idx == -1:
                    # This should never happen if methods are consistent
                    upstream_seq = ""
                    downstream_seq = ""
                else:
                    upstream_seq = f_seq[:idx]
                    downstream_seq = f_seq[idx + gene_len:]
            else:
                upstream_seq = ""
                downstream_seq = ""

            gene_data[gid] = {
                'gene_seq': g_seq,
                'upstream_seq': upstream_seq,
                'downstream_seq': downstream_seq
            }

        # Function to compute Jaccard similarity
        def jaccard_similarity(seq1, seq2, k):
            if len(seq1) < k or len(seq2) < k:
                return 0.0
            kmers1 = {seq1[i:i+k] for i in range(len(seq1) - k + 1)}
            kmers2 = {seq2[i:i+k] for i in range(len(seq2) - k + 1)}
            inter = len(kmers1 & kmers2)
            union = len(kmers1 | kmers2)
            return inter / union if union > 0 else 0.0

        # Write results
        with open(output_path, 'w') as outf:
            outf.write("gene_callers_id_1\tgene_callers_id_2\tgene_similarity\tupstream_similarity\tdownstream_similarity\tcombined_flank_similarity\n")

            total_pairs = len(gene_ids) * (len(gene_ids) - 1) // 2
            progress.new('Computing gene similarities', progress_total_items=total_pairs)

            count = 0
            for i in range(len(gene_ids)):
                gid1 = gene_ids[i]
                data1 = gene_data[gid1]
                for j in range(i+1, len(gene_ids)):
                    gid2 = gene_ids[j]
                    data2 = gene_data[gid2]

                    # Compute similarities
                    gene_sim = jaccard_similarity(data1['gene_seq'], data2['gene_seq'], kmer_size)
                    up_sim = jaccard_similarity(data1['upstream_seq'], data2['upstream_seq'], kmer_size)
                    down_sim = jaccard_similarity(data1['downstream_seq'], data2['downstream_seq'], kmer_size)

                    # Combined flanks
                    combined1 = data1['upstream_seq'] + data1['downstream_seq']
                    combined2 = data2['upstream_seq'] + data2['downstream_seq']
                    comb_sim = jaccard_similarity(combined1, combined2, kmer_size)

                    outf.write(f"{gid1}\t{gid2}\t{gene_sim:.6f}\t{up_sim:.6f}\t{down_sim:.6f}\t{comb_sim:.6f}\n")

                    count += 1
                    if count % 1000 == 0:
                        progress.increment(increment_to=count)

            progress.increment(increment_to=total_pairs)
            progress.end()

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
    return parser.get_args(parser)


if __name__ == '__main__':
    main()
