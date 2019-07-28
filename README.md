
# MetaErg
MetaErg is a stand-alone and fully automated metagenomic and metaproteomic data annotation pipeline. It bundles essential annotation tasks such as feature prediction, functional annotation with Hidden Markov Model (HMM) searches as well as blast and diamond searches. It estimates and visualizes quantitative taxonomic and pathway compositions of multiple metagenomic and proteomics samples using sequencing coverage and proteomics spectral counts, respectively. For visualization, MetaErg provides a HTML interface, bringing all annotation results together, and producing sortable and searchable tables, collapsible trees, and other graphic representations enabling intuitive navigation of complex data.

MetaErg analysis output demo page can be access at: https://xiaoli-dong.github.io/metaerg/

# Required perl modules
MetaErg was developed using perl, html, and javascript. It requires Perl 5.6.0 or higher and runs on Linux platform. Besides the perl core modules, it also requires the following perl modules to be installed:
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
MetaErg makes use of the following 3rd party dependencies and assumes these are on your system path:

* [ARAGORN](http://mbio-serv2.mbioekol.lu.se/ARAGORN):  a program to detect tRNA genes and tmRNA genes in nucleotide sequences 
* [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): The Basic Local Alignment Search Tool (BLAST) finds regions of	local similarity between sequences.
* [DIAMOND](https://github.com/bbuchfink/diamond): a new high-throughput program for aligning DNA reads or protein sequences against a protein reference database 
* [Hmmer](http://hmmer.org): HMMER is used for searching sequence databases for sequence homologs, and for making sequence alignments.
* [MinCED](https://github.com/ctSkennerton/minced): a program to find Clustered Regularly Interspaced Short Palindromic Repeats (CRISPRs) in full genomes or environmental datasets such as assembled contigs from metagenomes.
* [MinPath](http://omics.informatics.indiana.edu/MinPath): a parsimony approach for biological pathway reconstructions using protein family predictions, achieving a more conservative, yet more faithful, estimation of the biological pathways for a query dataset.
* [Prodigal](https://github.com/hyattpd/Prodigal): Fast, reliable protein-coding gene prediction for prokaryotic genomes.
* [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp): The program predicts the presence of signal peptides and the location of their cleavage sites in proteins from Archaea, Gram-positive Bacteria, Gram-negative Bacteria and Eukarya.
* [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm): a method for prediction transmembrane helices based on a hidden Markov model

# MetaErg reference DB
MetaErg requires external databases that need to be downloaded and unarchived
```
# Retrieve the prebuilt database
wget http://ebg.ucalgary.ca/metaerg/db.tar.gz -P $HOME
tar -xvzf $HOME/db.ar.tz
```
or built using the supplied script (see running with docker and installation sections). 

MetaErg databases were built based on the following publicly available databases:

* [casgene.hmm](https://www.nature.com/articles/nature21059)
* [FOAM](https://cbb.pnnl.gov/portal/software/FOAM.html)
* [metabolic-hmms](https://github.com/banfieldlab/metabolic-hmms)
* [Pfam](http://pfam.xfam.org)
* [SwissProt](https://www.uniprot.org/)
* [SILVA](https://www.arb-silva.de/download/archive/)
* [TIGRFAMS](http://tigrfams.jcvi.org/cgi-bin/index.cgi)
* [GTDBTK](https://github.com/Ecogenomics/GTDBTk)
* [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

# Running with docker
MetaErg docker image is host on the docker hub: https://hub.docker.com/r/xiaolidong/docker-metaerg. Due to licences permissions, this image does not contain [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp) and [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm). When running with docker image, "--sp --tm" options cannot be enabled.
```
# Get Docker image
docker pull xiaolidong/docker-metaerg

# Using the downloaded prebuilt database or build database using MetaErg supplied script. Building the database process will take a while to run:
docker run --shm-size 2g --rm -u $(id -u):$(id -g) -it -v my_local_dir:/data/ docker-metaerg setup_db.pl -o /data -v 132

#Running MetaErg with default options
docker run --shm-size 2g --rm -u $(id -u):$(id -g) -it -v my_local_dir:/data/ docker-metaerg metaerg.pl --dbdir /data/db /data/contig.fasta
```

# Installation
```
# Install metaerg to your home directory
git clone https://github.com/xiaoli-dong/metaerg.git $HOME/metaerg

# Using the downloaded prebuilt database
# or build database using MetaErg supplied script and the building process will take a while to run
$HOME/metaerg/bin/setup_db.pl -o $HOME -v 132
```

# Quick start:
The functionality provided by MetaErg can be accessed through the help menu:
```
>perl $HOME/metaerg/bin/metaerg.pl --help
```
Running MetaErg with the default parameters will output the final and intermediate results into a directory named metaerg.pl_ddmmyyyy directory
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db contig.fasta
```
Running MetaErg with user defined output directory and file names
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --outdir mydir --prefix mycontigs contig.fasta
```
With a user provided depth file, MetaErg can quantify the taxonomic, functional, and pathway compositions of multiple metagenomic samples. An example depth file was included in the "example" direcotry and it was genereated using "jgi_summarize_bam_contig_depths" from [MetaBat](https://bitbucket.org/berkeleylab/metabat) using BAM files, which are created by aligning the reads of each metagenomic sample separately to the contigs
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --depth demo.depth.txt demo.fasta
```
With a user provided protein expression level profile, MetaErg can also quantify functional, pathway profiles based on the active expressed protein genes from the metagenomic samples. An example protein expression profile was also included in the "example".
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --plevel demo.proteomics.txt demo.fna
```
# MetaErg utility scripts
MetaErg also includes some utility scripts to filtering contigs, add bin ids to the coding sequecnes, generate input for VizBin program. Some of the examples are listed as below:
```
#Filter out fasta format sequences shorter than a defined length
>perl $HOME/metaerg/bin/filterContigByLength.pl contig.fasta 500
```
The above command filters out the contigs <500bp from contig.fasta file
```
#Generate input files for VisBin program to visualize binning results
>perl $HOME/metaerg/bin/getVizBinInput.pl -d binning_dir
```
The "binning_dir" contains all the bin files. Each bin file is named in the format of "Bin.binid.fa" and binid is a number. Each bin file contains all the contigs binned together. The above command writes two files into the "binning_dir":  "binned_concat.fasta" and "binned_annotation.list".  In the VizBin application, the "binned_concat.fasta" will be uploaded to "File to Visualize" field and "binned_annotation.list" will be loaded to "Annotation file(optional)" field.
# Extract the subset of the MetaErg annotation results
Step1, extracting the annotations belonging to all the contigs contained in "subset.fasta" in gff format from the total dataset annotation :
```
#the total annotations are in mydir direcotry and the subset.fasta is a subset of the total input sequences to MetaErg annnotation
>perl $HOME/metaerg/bin/fastaContig2Gff.pl -c subset.fasta -g mydir/data/master.gff  > subset.gff
```
Step 2, generating the annotation results and html reports for the subset sequences
```
>perl $HOME/metaerg/bin/output_reports.pl  -g subset.gff -f subset.fasta -o mysubsetdir
```
# Add bin ids to the MetaErg generated files
Let's assume you are in "example" directory of the MetaErg installation and your binning results are in the "binning" directory.  "Bin.1.fa", "Bin.2.fa",  and "Bin.3.fa" files sitting in the "binning" directory contain all the fasta format contig sequences belonging to bin1, bin2, and bin3, respectivly.  
```
#Add bin id to the front of the protein coding sequence id in the format of "binid_" 
>perl $HOME/metaerg/bin/add_binid2cds.pl -d binning -c mydir/data/cds.faa -g mydir/data/master.gff
```
```
# Add bin ids to master.tsv file  as the first column
>perl $HOME/metaerg/bin/add_binid2master_dot_tsv.pl -d binning -t mydir/data/master.tsv
```
---
# MeteErg outputs
| Output        | Description|
|:--- |:--- |
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
|:---|:---|
|\*.tab.txt | Tab separated feature summary file|
|\*.json | Feature summary files in json formats, which are used to populate trees, tables, and sunburst visualizations in html reports.|
|\*.minpath| MinPath program generated pathway report files, on how to interpret the results, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt) |
|\*.minpath.details|MinPath program generated pathway report files, on how to interpret the results, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt) |
|master.gff| Total MetaErg annotation in  GFF3 format, It can be viewed directly in Artemis or IGV|
|master.tsv| A tab-separated file including all the selected annotation fields |
|master.tbl| A total MetaErg annotation file in feature table format | 
|master.stats.txt | A statistic summary file relating to the annotated features found |
|msampleName2shortName.txt | A tab-separated file has two lines. The first line is the tab-separated original sample names MetaErg gets from the coverage file provided by user and the second line is the tab-separated new shorter sample names generated by MetaErg |
|psampleName2shortName.txt |A tab-separated file has two lines. The first line is tab-separated original sample names MetaErg gets from the protein expression level file provided by user and the second line is the tab-separated new shorter sample names generated by MetaErg |
|cds.foam.\* | FOAM database searching output summary files in different formats|
|cds.metabolic.\* | Protein coding genes searching against [metabolic database](https://github.com/banfieldlab/metabolic-hmms) searching result summaries in different formats|
|cds.gene2ec.txt | Gene to ec number mapping file of the total metegenome samples|
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
