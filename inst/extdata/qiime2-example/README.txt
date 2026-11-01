Simulated data in the layout QIIME 2 writes. Not real data.

feature-table.tsv   the output of `biom convert --to-tsv` (50 samples)
taxonomy.tsv        the output of `qiime tools export` for a taxonomy artifact
sample-metadata.tsv a QIIME 2 metadata file with the #q2:types row (48 samples with metadata; 2 are controls not in the table)

Built by data-raw/make-qiime2-example.R. Try:  ?read_qiime2
