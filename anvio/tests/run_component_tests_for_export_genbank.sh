#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

files_dir=$(pwd)/$files

INFO "Setting up the export-genbank test directory"
mkdir -p $output_dir/export_genbank
cd $output_dir/export_genbank

# Create a mock database first (reusing the previous logic)
INFO "Generating a mock GenBank file"
cat > mock.gbk <<EOF
LOCUS       MOCK_GENOME              1000 bp    DNA     linear   BCT 18-MAY-2026
DEFINITION  Mock genome for testing.
ACCESSION   MOCK_01
VERSION     MOCK_01
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Mock bacterium"
     CDS             100..400
                     /locus_tag="MOCK_001"
                     /product="Healthy coding protein"
                     /translation="MTVKVGINGFGRIGR"
     tRNA            450..520
                     /locus_tag="MOCK_002"
                     /product="tRNA-Ala"
     rRNA            600..800
                     /locus_tag="MOCK_003"
                     /product="16S ribosomal RNA"
ORIGIN
        1 gacctacgtc gttttgcgcc atttatacgc cacatcgcgt gggtttgcgc ttggtaatgg
       61 ggcttagatg ctcctatgtg gtcctctccg gggtgtgttc taaggctaca aggatgtatt
      121 ggatattccc gtgcagcgac cacacggcgg ggtgatctag tggaaatttc tgatctcgat
      181 gtcgaaagct agtcttcctg tgggagccgt cattatactg tgcagctcta tcgtatgatt
      241 gctcagctta cttttgaatt cgcaactttc tgctagccgg ctccgccgat cgatgtctca
      301 atactgcaga tactaactcg cgatcgacac tcgggtgggt tttgcgatcc gagaagtgaa
      361 atttcgaagc tgcggggacc atctcggtag atctagtatt ctgaagggag gatggaactt
      421 agcttgggga ctgctacccg cggagggggg tctggataga ccacttggtc gtatactacg
      481 ggtttgccaa gagccgtcga atcagagaat ctgcataccg ctgtctaact tcctgaggag
      541 aacccctaaa tctcgtttac tttataactt aaagctaaat tattctgacg tgagggatcc
      601 acactatcag cagcgctgtt gtcttacctt taggggttta catcatccag gacgcaattt
      661 agtgtttgga gcctttttcc aagtacacca cctcggatag tcggtatgtg tacattagcc
      721 ttcgccttgt caccatcggc cacccacgtt gttagggttt caagtgggag cccagtcatc
      781 ggttttggta aatcgtacga cggtagtgcc atgcagtcat acacagagcc acaattaacc
      841 tcgacgagtt tggcacccac cgagtggagc tccgtagcac gcgaggttct tacagcaact
      901 gtactttaac ctgtgtccgg ttcgtgatac ggacagagga gtgtagcacc agagactgag
      961 ccggttgcaa attagtgctc aatgcgtgca ctgtttcata
//
EOF

INFO "Converting GenBank to anvi'o artifacts"
anvi-script-process-genbank -i mock.gbk -O MOCK

INFO "Generating contigs database"
anvi-gen-contigs-database -f MOCK-contigs.fa \
                          --external-gene-calls MOCK-external-gene-calls.txt \
                          -o MOCK.db \
                          --project-name "Test for export-genbank" \
                          --no-progress

INFO "Importing functions"
anvi-import-functions -c MOCK.db -i MOCK-external-functions.txt

INFO "Running anvi-export-genbank"
anvi-export-genbank -c MOCK.db -o exported.gbk

INFO "Checking exported GenBank file"
if [ ! -f exported.gbk ]; then
    echo "ERROR: Output file exported.gbk was not created"
    exit 1
fi

# Check for features in exported file
grep -q "CDS             100..400" exported.gbk || { echo "ERROR: CDS feature missing or coordinates wrong"; exit 1; }
grep -q "tRNA            450..520" exported.gbk || { echo "ERROR: tRNA feature missing or coordinates wrong"; exit 1; }
grep -q "rRNA            600..800" exported.gbk || { echo "ERROR: rRNA feature missing or coordinates wrong"; exit 1; }
grep -q "/product=\"Healthy coding protein\"" exported.gbk || { echo "ERROR: product annotation missing"; exit 1; }

INFO "Running anvi-export-locus with --export-genbank"
# Anchor on the tRNA (ID 1 usually)
anvi-export-locus -c MOCK.db --gene-caller-ids 1 -n 1,1 -O LOCUS --export-genbank

INFO "Checking exported locus GenBank file"
if [ ! -f LOCUS_0001.gbk ]; then
    echo "ERROR: Output file LOCUS_0001.gbk was not created"
    exit 1
fi

# The locus should contain CDS (100-400), tRNA (450-520), and rRNA (600-800)
# coordinates will be adjusted to the locus sequence.
grep -q "tRNA" LOCUS_0001.gbk || { echo "ERROR: tRNA missing from locus GenBank"; exit 1; }

INFO "SUCCESS: anvi-export-genbank and anvi-export-locus --export-genbank tests passed"
