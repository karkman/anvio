# anvi-compare-genes

Compares all genes in a contigs database using k-mer/Jaccard similarity for genes and their flanking regions.

Artifacts required: %(contigs-db)s
Artifacts produced: None

## Description

This program iterates through every pair of genes in a given %(contigs-db)s and computes the Jaccard similarity based on k-mer compositions of their sequences. It reports similarities for the gene sequences themselves, their upstream flanking regions, and their downstream flanking regions.

The program is useful for identifying genes with similar sequence contexts, which can be an indicator of horizontal gene transfer or conserved genomic neighborhoods.

### Overcoming the SQLite k-mer limit

Standard anvi'o programs (like `anvi-gen-contigs-database`) use a "wide" table format where every k-mer is a column. Because SQLite has a limit of 2,000 columns, these programs are typically limited to $k=5$.

**anvi-compare-genes** bypasses this limit by using an optimized "side-car" caching strategy. It does not create new columns in your %(contigs-db)s. Instead, it computes k-mers in memory and can optionally store them in a separate sparse SQLite database via the `--cache-file` flag. This allows you to use k-mer sizes up to **$k=13$** (and beyond) without hitting database limitations.

## Usage

{{ codestart }}
anvi-compare-genes -c %(contigs-db)s \
                   -o results.txt \
                   --kmer-size 4 \
                   --flank-length 100
{{ codestop }}

### Using a Cache File

For large k-mer sizes (e.g., $k=13$) or large metagenomes, k-mer tokenization can be slow. You can use a cache file to store these pre-computed sets for future runs:

{{ codestart }}
anvi-compare-genes -c %(contigs-db)s \
                   -o results.txt \
                   --kmer-size 13 \
                   --cache-file my_kmers.cache
{{ codestop }}

The output is a TAB-delimited file containing the following columns:

1. `gene_callers_id_1`: The ID of the first gene in the comparison.
2. `gene_callers_id_2`: The ID of the second gene in the comparison.
3. `gene_similarity`: Jaccard similarity of the gene sequences.
4. `upstream_similarity`: Jaccard similarity of the upstream flanking regions (5' end).
5. `downstream_similarity`: Jaccard similarity of the downstream flanking regions (3' end).
6. `combined_flank_similarity`: Jaccard similarity of the concatenated upstream and downstream flanking regions.

## Notes

- For genes on the reverse strand, the program correctly identifies the 5' (upstream) and 3' (downstream) ends and reverse-complements them before comparison.
- If a gene is too close to the start or end of a contig, the flanking regions will be truncated accordingly.
- The complexity of this program is O(N^2) where N is the number of genes. For large databases, this may take a significant amount of time. Use `--num-threads` to speed up the process.
