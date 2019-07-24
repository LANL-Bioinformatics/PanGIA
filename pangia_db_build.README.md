### Update database

##### Adding sequences of existing strains

PanGIA database is a standard multi-sequences FASTA file with a specific format in FASTA header. As long as the header of new sequences are modified to the rule below, you can append original database (FASTA). Note that the TAXIDs have to already be included in `taxonomy.tsv` or `taxonomy.custom.tsv`.

```
>[ACC#]|[SEQ_LENGTH]|[TAXID]|[FLAGS]
```

FLAGS is combined with [A|B|E|V|H] and [c|p|g].

 - 'A' for 'Archaea'
 - 'B' for 'Bacteria'
 - 'V' for 'Virus'
 - 'H' for 'Host'
 - 'c' for 'chromosom'
 - 'p' for 'plasmid'
 - 'g' for 'Phage'

For example, JIDR01000002.1, is a plasmid belongs to Francisella tularensis strain (taxid: 1341661) and the length of the sequence is 3194bp:

```
>JIDR01000002.1|3194|1341661|Bp
```

##### Adding sequences of new strains

If your update involved in generating new custom strains, you can use `pangia_db_build.py` to modify FASTA headers. A mapping table, `accession2taxid.tsv`, is necessary to facilitate converting acc# to taxonomy ID. The file can be downloaded at LANL's ftp server `ftp://ftp.lanl.gov/public/genome/PanGIA/v2/accession2taxid.tsv`.

Persume general taxonomy info file `taxonomy.tsv`, custom taxonomy info file `taxonomy.custom.tsv` and access# to taxid mapping table `accession2taxid.tsv` are stored in `database/` directory. Other files look like the setup below:

 - pangia_db_build.py
 - taxonomy.py
 - database / taxonomy.tsv
 - database / taxonomy.custom.tsv
 - database / accession2taxid.tsv
 - database / NCBI_genomes_111216_p_GRCh38.fa
 - sequences.fa

Let's say we want to add `sequences.fa` to current database `NCBI_genomes_111216_p_GRCh38.fa`. There are 3 major steps:

 - Modify FASTA header using `pangia_db_build.py`.

```sh
cat sequences.fa | pangia_db_build.py \
    --kingdom B \
    -dp ./database \
    --custax sequences.custom_taxonomy.tsv
    --update > sequences.pangia.fa
```
This step will generate two output, 1) modified FASTA file `sequences.pangia.fa` and 2) new custom taxonomy file `sequences.custom_taxonomy.tsv`.

 - Append new data to old database

```sh
cat sequences.pangia.fa >> database/NCBI_genomes_111216_p_GRCh38.fa
cat sequences.custom_taxonomy.tsv >> database/taxonomy.custom.tsv
```

 - Indexing new database

```
bwa index database/NCBI_genomes_111216_p_GRCh38.fa
```

### Build taxonomy lineage and acc# mapping table

Two files need to be generated `accession2taxid.tsv` and `taxonomy.tsv`. `accession2taxid.tsv`

Retrieve taxonomy from NCBI ftp:

```sh
cd ../
nohup rsync -auvzh --delete rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy taxonomy &
```

Produce an acc# to taxonomy ID mapping table for all nucleotide GenBank acc#. Note that the version number (.#) is revmoed from acc# to increase compatibility.

```sh
mkdir database
cd database
zcat ../taxonomy/accession2taxid/*nucl*.gz | cut -f 2,3 | sed 's/\.[[:digit:]]*//' | LC_ALL=C sort > accession2taxid.nucl.tsv
ln -s accession2taxid.nucl.tsv accession2taxid.tsv
```

Produce taoxnomy lineage table `taxonomy.tsv`:

```sh
tar -xzf ../taxonomy/taxdump.tar.gz
../scripts/extractTaxonomy.pl ../taxonomy > taxonomy.tsv
```

Then, we initial a `taxonomy.custom.tsv` file to store our custom taxonomies. Since there are no high-level ranks for viruses, the following text is added to `taxonomy.custom.tsv` to increase the resolutions.

```sh
131567	1	1	subroot	cellular organisms
35237	2	10239	phylum	dsDNA viruses, no RNA stage
35325	2	10239	phylum	dsRNA viruses
12333	2	10239	phylum	unclassified phages
12429	2	10239	phylum	unclassified viruses
12877	2	10239	phylum	Satellites
29258	2	10239	phylum	ssDNA viruses
35268	2	10239	phylum	Retro-transcribing viruses
186616	2	10239	phylum	environmental samples
439488	2	10239	phylum	ssRNA viruses
451344	2	10239	phylum	unclassified archaeal viruses
552364	2	10239	phylum	unclassified virophages
686617	2	10239	phylum	unassigned viruses
1425366	2	10239	phylum	Virus-associated RNAs
1714266	2	10239	phylum	Virus families not assigned to an order
```
