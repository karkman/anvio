# anvi-compare-genes

Compares all genes in a contigs database using k-mer/Jaccard similarity for genes and their flanking regions.

Artifacts required: %(contigs-db)s
Artifacts produced: None

## Description

This program iterates through pairs of genes in a given %(contigs-db)s and computes the Jaccard similarity based on k-mer compositions of their sequences. It reports similarities for the gene sequences themselves, their upstream flanking regions, and their downstream flanking regions.

The program is useful for identifying genes with similar sequence contexts, which can be an indicator of horizontal gene transfer or conserved genomic neighborhoods.

### Overcoming the SQLite k-mer limit

Standard anvi'o programs (like `anvi-gen-contigs-database`) use a "wide" table format where every k-mer is a column. Because SQLite has a limit of 2,000 columns, these programs are typically limited to $k=5$.

**anvi-compare-genes** bypasses this limit by using an optimized "side-car" caching strategy. It does not create new columns in your %(contigs-db)s. Instead, it computes k-mers in memory and can optionally store them in a separate sparse SQLite database via the `--cache-file` flag. This allows you to use k-mer sizes up to **$k=13$** (and beyond) without hitting database limitations.

### Comparison modes

You must specify exactly one comparison mode:

- `--compare-by-annotation-source`: Group genes by functional annotation and compare within each group. This is the recommended mode for most analyses.
- `--gene-caller-ids`: Compare a specific list of genes (one gene caller ID per line in a text file). For an all-against-all comparison of the entire database, provide a file with all gene IDs.

## Usage

{{ codestart }}
anvi-compare-genes -c %(contigs-db)s \
                   -o results.txt \
                   --compare-by-annotation-source COG_FUNCTION \
                   --kmer-size 4 \
                   --flank-length 100
{{ codestop }}

### Using a Cache File

For large k-mer sizes (e.g., $k=13$) or large metagenomes, k-mer tokenization can be slow. You can use a cache file to store these pre-computed sets for future runs:

{{ codestart }}
anvi-compare-genes -c %(contigs-db)s \
                   -o results.txt \
                   --compare-by-annotation-source COG_FUNCTION \
                   --kmer-size 13 \
                   --cache-file my_kmers.cache
{{ codestop }}

The cache file validates the k-mer size on subsequent runs and will raise an error if the k-mer size has changed.

## Output Format

The output is a TAB-delimited file containing the following columns:

1. `gene_callers_id_1`: The ID of the first gene in the comparison.
2. `gene_callers_id_2`: The ID of the second gene in the comparison.
3. `gene_similarity`: Jaccard similarity of the gene sequences.
4. `upstream_similarity`: Jaccard similarity of the upstream flanking regions (5' end).
5. `downstream_similarity`: Jaccard similarity of the downstream flanking regions (3' end).
6. `combined_flank_similarity`: Jaccard similarity of the concatenated upstream and downstream flanking regions.
7. `annotation_1`: Functional annotations for the first gene.
8. `annotation_2`: Functional annotations for the second gene.

### Annotation columns

- When using `--compare-by-annotation-source`, only the specified source(s) are shown.
- When NOT using `--compare-by-annotation-source`, annotations from **all available sources** in the database are fetched and displayed in the format: `SOURCE1: annotation | SOURCE2: annotation | ...`
- If a gene has no annotations, the column will be empty.

## Notes

- For genes on the reverse strand, the program correctly identifies the 5' (upstream) and 3' (downstream) ends and reverse-complements them before comparison.
- If a gene is too close to the start or end of a contig, the flanking regions will be truncated accordingly.
- The complexity of this program is O(N^2) where N is the number of genes. For large databases, this may take a significant amount of time. Use `--num-threads` to speed up the process.
- When using `--min-similarity`, pairs are pre-filtered using MinHash sketches before full comparison, which can significantly reduce runtime for large datasets.
