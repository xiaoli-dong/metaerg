
# MetaErg
MetaErg is a stand-alone and fully automated metagenome and metaproteome annotation pipeline. MetaErg bundles essential annotation tasks such as feature prediction, functional annotation with Hidden Markov Model (HMM) searches as well as blast and diamond searches. It estimates and visualizes quantitative taxonomic and pathway compositions of multiple metagenomic and proteomics samples using sequencing coverage and proteomics spectral counts, respectively. For visualization, MetaErg provides a HTML interface, bringing all annotation results together, and producing sortable and searchable tables, collapsible trees, and other graphic representations enabling intuitive navigation of complex data.

MetaErg is developed using perl, html, and javascript and it requires Perl 5.6.0 or higher. It has been tested on Linux platform. 


# Required perl modules
MetaErg requires the following perl modules, which are not included in the perl core modules, to be installed:
```
* Archive::Extract;
* Bio::Perl;
* Bio::DB::EUtilities
* DBD::SQLite
* DBI;
* File::Copy::Recursive
* LWP::Protocol::https
* SWISS::Entry;
* SWISS::KW;
```
# Third-party software
MetaErg depends on a list of third-party programs to do the gene predication and function annotations. Make sure that all MetaErg's dependencies to be instaled and are in your system's path. MetaErg depends on:

* [ARAGORN](http://mbio-serv2.mbioekol.lu.se/ARAGORN):  a program to detect tRNA genes and tmRNA genes in nucleotide sequences 
* [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): The Basic Local Alignment Search Tool (BLAST) finds regions of	local similarity between sequences.
* [DIAMOND](https://github.com/bbuchfink/diamond): a new high-throughput program for aligning DNA reads or protein sequences against a protein reference database 
* [Hmmer](http://hmmer.org): HMMER is used for searching sequence databases for sequence homologs, and for making sequence alignments.
* [MinCED](https://github.com/ctSkennerton/minced): a program to find Clustered Regularly Interspaced Short Palindromic Repeats (CRISPRs) in full genomes or environmental datasets such as assembled contigs from metagenomes.
* [MinPath](http://omics.informatics.indiana.edu/MinPath): a parsimony approach for biological pathway reconstructions using protein family predictions, achieving a more conservative, yet more faithful, estimation of the biological pathways for a query dataset.
* [Prodigal](https://github.com/hyattpd/Prodigal): Fast, reliable protein-coding gene prediction for prokaryotic genomes.
* [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp): The program predicts the presence of signal peptides and the location of their cleavage sites in proteins from Archaea, Gram-positive Bacteria, Gram-negative Bacteria and Eukarya.
* [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm): a method for prediction transmembrane helices based on a hidden Markov model


# Installation
```
#This command will install metaerg to your home directory
git clone https://github.com/xiaoli-dong/metaerg.git $HOME/metaerg
```

MetaErg require external database to assign the taxonomic, functinal, and pathway annotations to the predicted genes. There are two ways to obtain the metaerg databases: 

The exteral data can be download and unarchived:
```
#download and unarchive the metaerg database to your home directory. the database base will be sitting in $HOME/db 
wget http://ebg.ucalgary.ca/metaerg/db.tar.gz -P $HOME
tar -xvzf $HOME/db.tar.gz
```

The external data can also be build by metaerg supplied script:
```
#build metaerg datbase to your home direcoty. the database will be in $HOME/db
$HOME/metaerg/bin/setup_db.pl -o $HOME -v 132
```
MetaErg databases were built based on the following public available databases:

* [casgene.hmm](https://www.nature.com/articles/nature21059)
* [FOAM](https://cbb.pnnl.gov/portal/software/FOAM.html)
* [metabolic-hmms](https://github.com/banfieldlab/metabolic-hmms)
* [Pfam](http://pfam.xfam.org)
* [SwissProt](https://www.uniprot.org/)
* [SILVA](https://www.arb-silva.de/download/archive/)
* [TIGRFAMS](http://tigrfams.jcvi.org/cgi-bin/index.cgi)
* [GTDBTK](https://github.com/Ecogenomics/GTDBTk)
* [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

# Running MetaErg

```
# Checking MetaErg's command line options
> perl $HOME/metaerg/bin/metaerg.pl --help
```
```
# Running MetaErg with default parameters
perl ../bin/metaerg.pl test.fasta
```
\# Running MetaErg with signal peptide and cleavage site and transmembrane helics predication enabled. MetaErg will be slow down with those two feature enabled

> perl ../bin/metaerg.pl --sp --tm test.fasta

\# Running MetaErg with all the output to a user defined directory. By default, MetaErg will generate a output directory in the format of "locusttag_timestamp"

> perl ../bin/metaerg.pl --sp --tm --outdir metaerg_out --prefix test --locustag metaerg_test test.fasta

\# Overwrite previous results with --force option

> perl  ../bin/metaerg.pl  --force test.fasta

\# Running MetaErg with user provided coverage profile. With this supplied coverage profile, MetaErg will generate quantitative taxonomic, functional, and pathway compositions of multiple metagenomic samples.  The coverage profile must be in a well defined format in order for MetaErg  program to parse it correctly. For the required coverage profile format, you can refer to "cyano.depth.txt" file included in the "example" directory of the MetaErg installation. The example file can be generated using "**jgi_summarize_bam_contig_depths**" command from [MetaBat](https://bitbucket.org/berkeleylab/metabat) program after you finish the short reads to contig mapping process.

>perl ../bin/metaerg.pl --outdir metaerg_cyano --prefix cyano --locustag cyano --sp --tm --depth cyano.depth.txt cyano.fna

\# Running MetaErg with with user provided protein expression level profile. With protein expression level profile provided, MetaErg will generate functional, pathway profiles based on the proteins expressed in the metagenomic samples. For the required protein expression profile format, you can refer to "cyano.proteomics.txt" file included in the "example" of the MetaErg installation.

>perl ../bin/metaerg.pl --outdir cyano --prefix cyano --locustag cyano --sp --tm --depth cyano.depth.txt --plevel cyano.proteomics.txt cyano.fna

---
##Utility scripts

MetaErg also includes some utility scripts to filtering contigs, add bin ids to the coding sequecnes, generate input for VizBin program. Some of the examples are listed as below:

##### Filter out fasta format sequences shorter than a defined length ####

> perl ../bin/filterContigByLength.pl test.fasta 500

The above command filters out the sequences shorter than 500bp from the input test.fasta file

##### Generate input files for VisBin program to visualize binning results ##########

>perl ../bin/getVizBinInput.pl -d binning_dir

The "binning_dir" contains all the bin files. Each bin file is named in the format of "Bin.binid.fa" and binid is a number. Each bin file contains all the contigs binned together. The above command writes two files into the "binning_dir":  "binned_concat.fasta" and "binned_annotation.list".  In the VizBin application, the "binned_concat.fasta" will be uploaded to "File to Visualize" field and "binned_annotation.list" will be loaded to "Annotation file(optional)" field.

##### Extract MetaErg annotations for a subset of input contig sequences #####

Let's assume you are in the "example" directory of the MetaErg installation and "subset.fasta" file contains a subset of contig sequences from "cyano.fasta" file. The following two commands will generate the html reports for the contigs included in "subset.fasta" file

\# Step1, extracting the annotations belonging to all the contigs contained in "subset.fasta" in gff format :

>perl ../bin/fastaContig2Gff.pl -c subset.fasta -g cyano/data/master.gff  > subset.gff

\# Step 2, generating the annotation results and html reports to "metaerg_subset_output"

>perl ../bin/output_reports.pl  -g subset.gff -f subset.fasta -o metaerg_subset_output

##### Add bin ids to the MetaErg generated fiiles #####

Let's assume you are in "example" directory of the MetaErg installation and your binning results are in the "binning" directory.  "Bin.1.fa", "Bin.2.fa",  and "Bin.3.fa" files sitting in the "binning" directory contain all the fasta format contig sequences belonging to bin1, bin2, and bin3, respectivly.  

\# Add bin id to the front of the protein coding sequence id in the format of "binid_" 

>perl ../bin/add_binid2cds.pl -d binning -c cyano/data/cds.faa -g cyano/data/master.gff

\# Add bin ids to master.tsv file  as the first column

>perl ../bin/add_binid2master_dot_tsv.pl -d binning -t cyano/data/master.tsv

---
## MeteErg outputs
MetaErg writes all the output files into a user defined or MetaErg generated output directory. The following is an example of the MetaErg output file structure layout:

[[img src=metaerg_output_home.png alt=MetaErg_output_fist_level_layout]]

| Output        | Description|
|: --- |: --- |
| cyan.fna | Reformated and filtered fasta format input contig sequences |
| data | A directory contains all the MetaErg generated annotation summary files in differenct formats. Although the files have different suffix, they are all text files, whch can be opened in any text editor |
| html | A directory contains all the HTML pages for various type of  HTML reports |
| images | A directory contains all the image files for the html reports such as logo, banner|
| index.html | An interactive HTML report, which links all the MetaErg annotation results together |
| js | A directory contains all the required Javascript libraries for the interactive html reports |
| style.css | A HTML style sheet, which controls the look of the html reports |
| tmp | A dirctory contains all the MetaErg intermediate outputs. It is also useful when MetaErg fails in the middle of the run. With this directory in place, when you restart the job using the exact same parameters after MetaErg failling, MetaErg will start from the place it failed.  After MetaErg job finishes successfully, this directory can be deleted before you transfer the results to your local computers|

### data directory list

| File|Description|
    |---|---|
     |\*.tab.txt | Tab separated feature summary file|
    |\*.json | Feature summary files in json formats, which are used to populate trees, tables, and sunburst visualizations in html reports.|
     |\*.minpath| MinPath program generated pathway report files, on how to interpret the results, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt) |
    |\*.minpath.details|MinPath program generated pathway report files, on how to interpret the results, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt) |
    |master.gff| Total MetaErg annotation in  GFF3 format, It can be viewed directly in Artemis or IGV|
    | master.tsv| A tab-separated file including all the selected annotation fields |
    |master.tbl| A total MetaErg annotation file in feature table format | 
    | master.stats.txt | A statistic summary file relating to the annotated features found |
    | msampleName2shortName.txt | A tab-separated file has two lines. The first line is the tab-separated original sample names MetaErg gets from the coverage file provided by user and the second line is the tab-separated new shorter sample names generated by MetaErg |
    |psampleName2shortName.txt |A tab-separated file has two lines. The first line is tab-separated original sample names MetaErg gets from the protein expression level file provided by user and the second line is the tab-separated new shorter sample names generated by MetaErg |
    |cds.foam.\* | FOAM database searching output summary files in different formats|
    | cds.metabolic.\* | Protein coding genes searching against [metabolic database](https://github.com/banfieldlab/metabolic-hmms) searching result summaries in different formats|
    | cds.gene2ec.txt | Gene to ec number mapping file of the total metegenome samples|
    |cds.gene2ec.M\*.txt| Gene to ec number mapping file of each input metagenome sample|
    |cds.gene2ec.plevel.p\*.txt| Gene to ec number mapping file, which is generated based on only the expressed protein from each sample|
    |cds.gene2ec.minpath\*|MetaCyc pathway report based on the total metagenome |
    |cds.gene2ec.M\*.minpath\*|MetaCyc pathway report of each metagenome sample|
    |cds.gene2ec.plevel.P\*.minpath\*|MetaCyc pathway report from only the expressed proteins of an individual metagenome sample|
    |cds.gene2ko.txt | Gene to ko number mapping files of the total metagenome samples|
    |cds.gene2ko.M\*.txt|Gene to ko number mapping file of the total metagenome samples|
    |cds.gene2ko.plevel.p\*.txt|Gene to ko number mapping file,  which is generated based on only the expressed protein from each sample |
    |cds.gene2ko.minpath\*|KEGG pathway report of the total metagenome samples|
    |cds.gene2ko.M\*.minpath|KEGG pathway report of each metagenome sample|
    |cds.gene2ko.plevel.P\*.minpath|KEGG pathway report from only the expressed proteins of an individual metagenome sample|
    |taxon.cds.\*| Taxnomy assignment sumary reports of the total metagenome smaple in different formats based on the coding sequences|
    |taxon.\*SrRNA.\*|The rRNA gene taxonomy assignment summary reports in different formats |
    |cds.gene2sport.\*| Protein coding gene searching against SwissPort database result summaries in different format |
    |cds.gene2pfam.\*|Protein coding sequence searching against pfam database result summaries in different format|
    |cds.gene2tigrfam.\*|Protein coding sequence searching against tigrfam database result summaries in different format|
    |rRNA.\*| Predicated rRNA genes (including both SSU-rRNA and LSU-rRNA) and their associcated taxonomy summary reports|
    |CRISPR.\*| Predicated CRISPR and its associated taxonomy summary reports|
    |\*.ffn |FASTA format nucleotide sequence files for the MetaErg predicated features (tRNA, crispr, cds, rRNA)|
    |cds.faa | Predicated protein coding amino acid sequence in FASTA format |
