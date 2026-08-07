This program lets you export the contents of your %(contigs-db)s as a GenBank file.

It will include all contig sequences and their associated gene calls (CDS, tRNA, and rRNA), as well as any functional annotations that have been imported into the database.

### Basic usage

You can export your entire database to a single GenBank file like this:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o my_genome.gbk
{{ codestop }}

### Exporting specific contigs or genes

If you only want to export specific parts of your database, you can use the `--contigs-of-interest` or `--gene-caller-ids` flags:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o subset.gbk \
                    --contigs-of-interest my_contigs.txt
{{ codestop }}

### Controlling functional annotations

By default, anvi'o will include all functional annotations in the GenBank file (under `note` and `db_xref` qualifiers). If you only want to include specific sources, you can use the `--annotation-sources` flag:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o annotated.gbk \
                    --annotation-sources COG20_FUNCTION,KOfam
{{ codestop }}
