
# MetaErg
MetaErg is a stand-alone and fully automated metagenomic and metaproteomic data annotation pipeline. It bundles essential annotation tasks such as feature prediction, functional annotation wit
h Hidden Markov Model (HMM) searches as well as blast and diamond searches. It estimates and visualizes quantitative taxonomic and pathway compositions of multiple metagenomic and proteomics s
amples using sequencing coverage and proteomics spectral counts, respectively. For visualization, MetaErg provides a HTML interface, bringing all annotation results together in sortable and se
archable tables, collapsible trees, and other graphic representations, enabling intuitive navigation of complex data.

A MetaErg analysis output demo page can be found at: https://xiaoli-dong.github.io/metaerg/

A MetaErg Docker application can be found here: https://hub.docker.com/r/xiaolidong/docker-metaerg

Please cite the following: [Dong X and Strous M (2019) An Integrated Pipeline for Annotation and Visualization of Metagenomic Contigs. Front. Genet. 10:999. doi: 10.3389/fgene.2019.00999](http
s://www.frontiersin.org/articles/10.3389/fgene.2019.00999/full)

# Required perl modules
If you do not use Docker, you will need to install perl modules. MetaErg requires Perl 5.6.0 or higher and runs on Linux platforms. Besides the perl core modules, it also requires the followin
g perl modules to be installed:
```
* Archive::Extract;
* Bio::Perl;
* Bio::DB::EUtilities
* DBD::SQLite
* DBI;
* File::Copy::Recursive
* HTML::Entities
* LWP::Protocol::https
* SWISS::Entry;
* SWISS::KW;
```
# Third-party software
If you do not use Docker, you will need to install the following 3rd party dependencies and make sure they are on your system path:

* [ARAGORN](http://mbio-serv2.mbioekol.lu.se/ARAGORN):  a program to detect tRNA genes and tmRNA genes in nucleotide sequences
* [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): The Basic Local Alignment Search Tool (BLAST) finds regions of  local similarity between
 sequences.
* [DIAMOND](https://github.com/bbuchfink/diamond): a new high-throughput program for aligning DNA reads or protein sequences against a protein reference database
* [Hmmer](http://hmmer.org): HMMER is used for searching sequence databases for sequence homologs, and for making sequence alignments.
* [MinCED](https://github.com/ctSkennerton/minced): a program to find Clustered Regularly Interspaced Short Palindromic Repeats (CRISPRs) in full genomes or environmental datasets such as asse
mbled contigs from metagenomes.
* [MinPath](http://omics.informatics.indiana.edu/MinPath): a parsimony approach for biological pathway reconstructions using protein family predictions, achieving a more conservative, yet more
 faithful, estimation of the biological pathways for a query dataset.
* [Prodigal](https://github.com/hyattpd/Prodigal): Fast, reliable protein-coding gene prediction for prokaryotic genomes.
* [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp): The program predicts the presence of signal peptides and the location of their cleavage sites in proteins from Archaea, Gram-
positive Bacteria, Gram-negative Bacteria and Eukarya.
* [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm): a method for predicting transmembrane helices based on a hidden Markov model

# MetaErg reference DB
MetaErg requires external databases that need to be downloaded and unarchived
```
# Retrieve the prebuilt database
wget http://ebg.ucalgary.ca/metaerg/db.tar.gz -P $HOME
tar -xvzf $HOME/db.tar.gz
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
The MetaErg docker image is hosted on the docker hub: https://hub.docker.com/r/xiaolidong/docker-metaerg. Due to licencing permissions, this image does not contain [SignalP](http://www.cbs.dtu
.dk/cgi-bin/nph-sw_request?signalp) and [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm). When running with docker image, "--sp --tm" options cannot be enabled.
```
# Get the Docker image
docker pull xiaolidong/docker-metaerg

# Databases and Docker
With Docker, you can either use the downloaded prebuilt database or build a database with the command below. Building the database process will take a while to run:
docker run --shm-size 2g --rm -u $(id -u):$(id -g) -it -v my_local_dir:/data/ xiaolidong/docker-metaerg setup_db.pl -o /data -v 132

#Running MetaErg with default options
docker run --shm-size 2g --rm -u $(id -u):$(id -g) -it -v my_local_dir:/data/ xiaolidong/docker-metaerg metaerg.pl --dbdir /data/db --outdir /data/my_metaerg_output /data/contig.fasta
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
Running MetaErg with the default parameters will output the final and intermediate results into a directory named metaerg.pl_ddmmyyyy. Without --dbdir option, MetaErg will look for the databas
e directory "db" inside the metaerg directory
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db contig.fasta
```
Running MetaErg with user defined output directory and file names:
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --outdir mydir --prefix mycontigs contig.fasta
```
With a user provided "depth file", MetaErg can quantify the taxonomic, functional, and pathway compositions of multiple metagenomic samples. An example "depth file" is included in the "example
" directory. The depth file was generated with the script "jgi_summarize_bam_contig_depths" from [MetaBat](https://bitbucket.org/berkeleylab/metabat) using BAM files. BAM files are created by
aligning the reads of each metagenomic sample to the assembled contigs, using a program such as BBMap or bwa.
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --depth demo.depth.txt demo.fasta
```
<!--
With a user provided protein expression level profile, MetaErg can also quantify functional, pathway profiles based on actively expressed proteins. An example protein expression profile is als
o included in the "example" folder.
```
>perl $HOME/metaerg/bin/metaerg.pl --dbdir $home/db --plevel demo.proteomics.txt demo.fna
```
-->
# Utility scripts
MetaErg includes some utility perl scripts and they can be used to filter contigs by length, add bin identifiers to predicted coding sequences:
```
#Filter out contig sequences shorter than 500bp
>perl $HOME/metaerg/bin/filterContigByLength.pl contig.fasta 500
```
Let's assume you have annotated all the contigs of your metagenome. MetaErg can extract the subset of annotation results and produce html summary pages for an individual bin as follows:
```
#Extract annotation results for individual bins
#Step1, extracting the gff-format annotations for the contigs included in "mybin.nt.fasta" from the total metaerg dataset annotation:
>perl $HOME/metaerg/bin/fastaContig2Gff.pl -c mybin.nt.fasta -g mydir/data/all.gff  > mybin.gff

#Step 2, generating the html reports for the extracted contig subset
>perl $HOME/metaerg/bin/output_reports.pl  -g mybin.gff -f subset.fasta -o mybindir
```
Let's assume mybindir contains many nucleotide fasta files, one for each bin: Bin.1.fa", "Bin.2.fa", "Bin.3.fa"... files. The following commands will:

```
#Add bin id to the fasta format of the protein coding sequence and protein coding sequence id will be in the format of "binid_geneid"
>perl $HOME/metaerg/bin/add_binid2cds.pl -d binning -c mydir/data/cds.faa -g mydir/data/all.gff

#Add bin ids to master.tsv file, as the first column
>perl $HOME/metaerg/bin/add_binid2master_dot_tsv.pl -d binning -t mydir/data/master.tsv
```
# MetaErg output directory layout
MetaErg writes the output files into a user defined or MetaErg autogenerated output directory. An example of the MetaErg output directory layout is as following:

| Output        | Description|
|:--- |:--- |
| \*.fna | Reformated and filtered fasta format input contig sequences |
| data | A directory containing all the MetaErg generated annotation summary files in different formats. Although the files have different suffixes, they are all text files and can be opened i
n any text editor |
| html | A directory containing all the HTML pages for various type of  HTML reports and visualizations |
| images | A directory containing all the image files for the html reports such as logo, banner|
| index.html | An interactive HTML summary report page, which links all the MetaErg annotations and visualizations together |
| js | A directory containing all the required Javascript libraries for the interactive html reports |
| style.css | A HTML style sheet, which controls the look of the html reports |
| tmp | A dirctory containing all the MetaErg intermediate outputs. It is useful when MetaErg fails in the middle of the run. With this directory in place, when you restart metaerg using the e
xact same parameters, MetaErg will start from the place it failed.  After a MetaErg job finishes successfully, this directory can be deleted before you transfer the results to your local compu
ters|

# MetaErg data directory file descriptions

| File|Description|
|:---|:---|
|all.gff| Total MetaErg annotation in  GFF3 format, It contains all the gene prediction, annotation information of the analysis|
|master.gff.txt| Total MetaErg annotation in  GFF3 format, It's a simplified version of the all.gff file and can be viewed directly in Artemis or IGV|
|master.tsv.txt| A tab-separated file including all the selected annotation fields |
|master.tbl.txt| A total MetaErg annotation file in feature table format |
|master.stats.txt | A statistic summary file relating to the annotated features found |
|16SrRNA.ffn|FASTA format nucleotide sequence file of the 16S rRNA genes|
|18SrRNA.ffn|FASTA format nucleotide sequence file of the 18S rRNA genes|
|23SrRNA.ffn|FASTA format nucleotide sequence file of the 23S rRNA genes|
|28SrRNA.ffn|FASTA format nucleotide sequence file of the 28S rRNA genes|
|5SrRNA.ffn|FASTA format nucleotide sequence file of the 5S rRNA genes|
|tRNA.ffn|FASTA format nucleotide sequence file of the tRNA genes|
|cds.faa|FASTA format amino acid sequence file of the protein coding genes|
|cds.ffn|FASTA format nucleotide sequence file of the protein coding genes|
|crispr.ffn|FASTA format nucleotide sequence file of the CRISPRs|
|crispr.tab.txt|Tab separated CRISPR summary file|
|rRNA.tab.txt|Tab separated rRNA gene summary file|
|tRNA.tab.txt|Tab separated tRNA gene summary file|
|cds.gene2casgene.tab.txt|Protein coding gene annotation summary with casgene.hmm database searching|
|cds.gene2ec.tab.txt|Protein coding gene annotation summary with enzyme assignments|
|cds.gene2genomedb.tab.txt|Protein coding gene annotation summary with GenomeDB database searching|
|cds.gene2go.tab.txt|Protein coding gene annotation summary with GO term assignments|
|cds.gene2ko.tab.txt|Protein coding gene annotation summary with KEGG orthology assignments|
|cds.gene2metabolic.tab.txt|Protein coding gene annotation summary with metabolic-hmms database searching|
|cds.gene2pfam.tab.txt|Protein coding gene annotation summary with pfam database searching|
|cds.gene2sprot.tab.txt|Protein coding gene annotation summary with swiss-prot database searching|
|cds.gene2tigrfam.tab.txt|Protein coding gene annotation summary with tigrfams database searching|
|cds.gene2kegg.tab.txt|Protein coding gene annotation summary with KEGG pathway assignments|
|cds.gene2metacyc.tab.txt|Protein coding gene annotation summary with MetaCyc pathway assignments|
|cds.gene2ec.mapping.txt|Protein coding gene to ec number mapping file, which was used for MetaCyc pathway prediction|
|cds.gene2ko.mapping.txt|Protein coding gene to ec number mapping file, which was used for KEGG pathway prediction|
|ec.cds.profile.tab.txt|Protein coding gene funciton profile based on enzyme assignments|
|go.cds.profile.tab.txt|Protein coding gene funciton profile based on GO term assignments|
|ko.cds.profile.tab.txt|Protein coding gene funciton profile based on KO assignments|
|metabolic.cds.profile.tab.txt|Protein coding gene funciton profile based on the metabolic marker gene annotations|
|pfam.cds.profile.tab.txt|Protein coding gene funciton profile based on pfam annotations|
|taxon.cds.profile.tab.txt|Taxonomic profile based on protein coding gene taxonomic classification|
|taxon.lsu.profile.tab.txt|Taxonomic profile based on LSU rNRA gene taxonomic classification|
|taxon.ssu.profile.tab.txt|Taxonomic profile based on SSU rNRA gene taxonomic classification|
|tigrfam.cds.profile.tab.txt|Protein coding gene funciton profile based on tigrfam annotations|
|kegg.cds.profile.tab.txt|KEGG pathway profile predicated using MinPath program|
|metacyc.cds.profile.tab.txt|MetaCyc pathway profile predicated using MinPath program|
|cds.gene2ec.minpath|MinPath program output for predicting MetaCyc pathways, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt)|
|cds.gene2ec.minpath.details|MinPath program output for predicting MetaCyc pathways, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt)|
|cds.gene2ko.minpath|MinPath program output for predicting KEGG pathways, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt)|
|cds.gene2ko.minpath.details|MinPath program output for predicting KEGG pathways, please refer to [how to read MinPath report](http://omics.informatics.indiana.edu/MinPath/readme.txt)|
|\*.json | Feature summary files in json formats, which are used to populate trees, tables, and sunburst visualizations in html reports|
