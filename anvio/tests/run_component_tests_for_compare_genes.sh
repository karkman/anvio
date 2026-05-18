#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

files_dir=$(pwd)/$files

INFO "Setting up the compare-genes test directory"
mkdir -p $output_dir/compare_genes
cd $output_dir/compare_genes

INFO "Generating a new contigs database"
anvi-gen-contigs-database -f $files_dir/contigs.fa \
                          -o CONTIGS.db \
                          --project-name "Test for compare-genes" \
                          --no-progress \
                          $thread_controller

INFO "Running anvi-compare-genes"
anvi-compare-genes -c CONTIGS.db \
                   -o results.txt \
                   --kmer-size 4 \
                   --flank-length 100 \
                   $thread_controller

INFO "Checking output"
if [ ! -f results.txt ]; then
    echo "ERROR: Output file results.txt was not created"
    exit 1
fi

if [ ! -s results.txt ]; then
    echo "ERROR: Output file results.txt is empty"
    exit 1
fi

# Check header (tabs)
header=$(head -n 1 results.txt)
expected_header="gene_callers_id_1	gene_callers_id_2	gene_similarity	upstream_similarity	downstream_similarity	combined_flank_similarity	annotations_1	annotations_2"
if [ "$header" != "$expected_header" ]; then
    echo "ERROR: Header mismatch"
    echo "Got:      $header"
    echo "Expected: $expected_header"
    exit 1
fi

# Check if there are at least some comparisons (CONTIGS.db from contigs.fa should have many genes)
num_lines=$(wc -l < results.txt)
if [ $num_lines -lt 2 ]; then
    echo "ERROR: Not enough results in output file"
    exit 1
fi

INFO "SUCCESS: anvi-compare-genes test passed"
