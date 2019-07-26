#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use File::Fetch;
use lib "$FindBin::Bin";
use Archive::Extract;
use SWISS::Entry;
use SWISS::KW;
use XML::Simple;
use XML::Parser;
use File::Copy;
use Getopt::Long;
use Time::Piece;
use File::Basename;
use Cwd 'abs_path';
my ($outdir, $s_version);
$s_version = "132";
&GetOptions(
    "o=s" =>\$outdir,
    "v=s" =>\$s_version
    );

($outdir) ||
    die "usage: $0 OPTIONS
where options are:\n -o  <data output direcoty><-v silva database release version, for example 132\n";

my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";
my $sqlfile = "$bin/../metaerg.sql";
my $DBDIR = "$outdir/db";
my $tmp_dir = "$outdir/tmp";
my $diamond_dir = "$DBDIR/diamond";
my $protein_hmm_dir = "$DBDIR/hmm";
my $rna_hmm_dir = "$DBDIR/hmm/rna";
my $sqlite_dir = "$DBDIR/sqlite3";
my $blast_dir = "$DBDIR/blast";

msg("construct db directories");
runcmd("mkdir -p \Q$outdir\E") if (! -d $outdir);
runcmd("mkdir -p \Q$DBDIR\E") if (! -d $DBDIR);
runcmd("mkdir -p \Q$tmp_dir\E") if (! -d $tmp_dir);
runcmd("mkdir -p \Q$diamond_dir\E") if(! -d $diamond_dir);
runcmd("mkdir -p \Q$protein_hmm_dir\E") if (! -d $protein_hmm_dir);
runcmd("mkdir -p \Q$sqlite_dir\E") if (! -d $sqlite_dir);
runcmd("mkdir -p \Q$blast_dir\E") if (! -d $blast_dir);

#build_rRNAFinder_hmmdb();
build_rRNAFinder_txondb();
build_uniprot_sprot_db();
build_pfam_db();
build_tigrfam_db();
build_go_db();
build_foam_db();
build_genomedb();
build_metabolic_hmmdb();
build_casgene_hmmdb();
build_sqlite_db($sqlfile);
sub build_uniprot_sprot_db{

    my @files = ("uniprot_sprot.fasta.gz", "uniprot_sprot.dat.gz", "reldate.txt");

    foreach my $file (@files){
	msg("Start to fetch $file");
	if(! -e "$tmp_dir/$file"){
	    # build a File::Fetch object
	    my $ff = File::Fetch->new(uri => "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$file");

	    #fetch the uri to local directory
	    my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	}
    }

    if(! -e "$diamond_dir/uniprot_sprot.dmnd"){
	my $cmd = "diamond makedb --tmpdir /dev/shm --in $tmp_dir/uniprot_sprot.fasta.gz -d $diamond_dir/uniprot_sprot";
	msg("Start running:$cmd");
	runcmd($cmd);
    }

    #unzip, extract info to tabular format
    if(! -e "$tmp_dir/uniprot_sprot.dat"){
	msg("Uncompressing $tmp_dir/uniprot_sprot.dat.gz");
	gunzip "$tmp_dir/uniprot_sprot.dat.gz" => "$tmp_dir/uniprot_sprot.dat" or die "gunzip failed: $GunzipError\n";
    }
    #$cmd = "$^X $bin/uniprot_dat2tbl.pl $tmp_dir/uniprot_sprot.dat > $sqlite3_dir/sprot.sqldb.txt";
    #mgs("Running $cmd");
    #runcmd($cmd);
    msg("Start to parse $tmp_dir/uniprot_sprot.dat");
    if(! -e "$tmp_dir/sprot.sqldb.tsv"){
	parse_uniprot_sport_dat("$tmp_dir/uniprot_sprot.dat", "$tmp_dir/sprot.sqldb.tsv");
    }
}

sub parse_uniprot_sport_dat{
    my ($dat, $output) = @_;

    open(OUT, ">$output") or die "Could not open $output to write, $!\n";
    open(DAT, $dat) or die "Could not open $dat to read, $!\n";
    $/ = "\n//\n";

    my $HYPO = 'hypothetical protein';

    local $/ = "\n//\n";
    my $in=0;
    my $out=0;

    while (<DAT>)
    {
	# Read the entry
	my $entry = SWISS::Entry->fromText($_);
	$in++;

	my $gn = $entry->GNs->getFirst() || '';
	# print "gn=$gn\n";

	my $ec = '';
	my $desc = 'hypothetical protein';

	if (1) {
	    for my $de ($entry->DEs->elements) {
		my @inces = $entry->DEs->Includes->elements;
		foreach my $ince (@inces){
		    foreach my $ide ($ince->elements){
			if ($ide->type eq 'EC') {
			    my $iec = $ide->text;
			    $iec =~ s/^\D*//;
			    $ec .="$iec;";
			}
			elsif ($ide->type eq 'Full' && $ide->category eq 'RecName' ) {
			    $desc .= ";" . $ide->text;
			}
		    }
		}

		if ($de->type eq 'EC') {
		    $ec = $de->text;
		    $ec =~ s/^\D*//;
		}
		elsif ($de->type eq 'Full' and $de->category eq 'RecName') {
		    $desc = $de->text;
		    if ($desc =~ m/^UPF\d|^Uncharacterized protein|^ORF|^Protein /) {
			$desc = $HYPO;
		    }
		}
		last if $desc and $ec;  # we have some data now, exit out
	    }
	}
	$desc ||= $HYPO;

	# skip hypthetical proteins, unless user has overridden this with --hypo
	#next if !$hypo and $desc eq $HYPO;

	#DRs
	my $biocyc = '';
	my $kegg = '';
	my $unipathway = '';
	my $eupathdb = '';
	my $go = '';
	my $reactome = '';
	my $KO = '';
	#my $cazy = '';


	if (1) {
	    for my $de ($entry->DRs->elements) {
		#print STDERR Dumper($de);
		#print STDERR join("\t", @$de), "\n";

		if($de->[0] eq "BioCyc"){
		    $biocyc .= $de->[1] .";";
		}
		elsif($de->[0] eq "KEGG"){
		    $kegg .= $de->[1] .";";
		}
		elsif($de->[0] eq "KO"){
		    $KO .= $de->[1] .";";
		}

		elsif($de->[0] eq "UniPathway"){
		    $unipathway .=  $de->[1] .";";
		}

		elsif($de->[0] eq "GO"){
		    $go .= $de->[1] .";";
		}

	    }
	}


	$ec ||= "";
	$gn ||= "";
	$desc ||= "";


	my $id = "sp|" . $entry->AC . "|" . $entry->ID;
	my $oc = join(";", @{$entry->OCs->list});
	my $os = ${$entry->OSs->list}[0]->text;
	#print STDERR Dumper($entry);
	$desc =~ s/;$//;
	$KO =~ s/;$//;
	$kegg =~ s/;$//;
	$go =~ s/;$//;
	$ec =~ s/;$//;
	print  OUT "$id\t$KO\t$kegg\t$ec\t$desc\t$go\t$os\t$oc\n";
	$out++;

    }

}


sub build_pfam_db{

    my @files = ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz", "Pfam.version.gz", "userman.txt");

    foreach my $file (@files){

	if(! -e "$tmp_dir/$file"){

	    #build a File::Fetch object
	    my $ff = File::Fetch->new(uri => "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/$file");

	    #fetch the uri to local directory
	    msg("Start fetching $file");
	    my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	}
    }

    if(! -e "$tmp_dir/Pfam-A.hmm.dat"){
	msg("Uncompressing $tmp_dir/Pfam-A.hmm.dat.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/Pfam-A.hmm.dat.gz");
	$ae->extract(to => $tmp_dir);
    }

    if(! -e "$protein_hmm_dir/Pfam-A.hmm"){
	msg("Uncompressing $tmp_dir/Pfam-A.hmm.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/Pfam-A.hmm.gz");
	$ae->extract(to => $protein_hmm_dir);
    }
    if(! -e "$tmp_dir/pfam2go"){
	my $ff = File::Fetch->new(uri => "http://geneontology.org/external2go/pfam2go");
	msg("Fetching pfam2go");
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }
    if(! -e "$tmp_dir/Pfam-A.hmm.sql.db.tsv"){
	open(PFAM, "$tmp_dir/Pfam-A.hmm.dat") or die "Could not open $tmp_dir/Pfam-A.hmm.dat to read, $!\n";
	open(SQL, ">$tmp_dir/Pfam-A.hmm.sql.db.tsv")  or die "Could not open $tmp_dir/Pfam-A.hmm.sql.db.tsv to write, $!\n";
	
	
	$/ = "\n# STOCKHOLM";
	
	while(<PFAM>){
	    if(/ID\s+(\S+).*?AC\s+(\w+).*?\n.*DE\s+(\S.*?)\n.*?TP\s+(\S.*?)\n.*?ML\s+(\d+)/s){
		print SQL "$2\t$1\t$3\t$4\t$5\n";
	    }
	    
	}
	close(SQL);
	close(PFAM);
	
	$/ = "\n";
    }
    
    if(! -e "$tmp_dir/pfam2go.sqldb.tsv"){
    open(SQL, ">$tmp_dir/pfam2go.sqldb.tsv")  or die "Could not open $tmp_dir/pfam2go.sqldb.tsv to write, $!\n";
    open(INPUT, "$tmp_dir/pfam2go") or die "Could not open $tmp_dir/pfam2go to read, $!\n";

    while(<INPUT>){
	chomp;
	next if /^!/;
	if(/Pfam:(\S+).*?;\s+(GO\S+)$/){
	    print SQL "$1\t$2\n";
	}
    }
    close(INPUT);
    close(SQL);
    }
}

sub build_tigrfam_db{


    my @files = ("TIGRFAMs_15.0_HMM.LIB.gz", "TIGRFAMS_GO_LINK", "TIGRFAMs_15.0_INFO.tar.gz", "TIGR_ROLE_NAMES", "TIGRFAMS_ROLE_LINK", "RELEASE_NOTE_15.0");

    foreach my $file (@files){

	if(! -e "$tmp_dir/$file"){

	    #build a File::Fetch object
	    my $ff = File::Fetch->new(uri => "ftp://ftp.jcvi.org/pub/data/TIGRFAMs/$file");

	    #fetch the uri to local directory
	    msg("Fetching $file");
	    my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	    #my $cmd = "wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/$file -O $tmp_dir";
	    #system($cmd) == 0 or err("Could not run command:$cmd");
	}
    }

    #extract archive files
    if(! -e "$protein_hmm_dir/TIGRFAMs.hmm"){
	msg("Uncompressing $tmp_dir/TIGRFAMs_15.0_HMM.LIB.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/TIGRFAMs_15.0_HMM.LIB.gz");
	$ae->extract(to => "$protein_hmm_dir/TIGRFAMs.hmm");
    }
    if(! -e "$tmp_dir/TIGRFAMs_15.0_INFO"){
	if(! -d "$tmp_dir/TIGRFAMs_15.0_INFO"){
	    msg("Creating $tmp_dir/TIGRFAMs_15.0_INFO folder");
	    runcmd("mkdir -p $tmp_dir/TIGRFAMs_15.0_INFO");
	}
	msg("Uncompressing $tmp_dir/TIGRFAMs_15.0_INFO.tar.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/TIGRFAMs_15.0_INFO.tar.gz");
	$ae->extract(to => "$tmp_dir/TIGRFAMs_15.0_INFO");
    }
    if(! -e  "$tmp_dir/TIGRFAMs.hmm.sqldb.tsv"){
	tigrfam_id2info_table("$protein_hmm_dir/TIGRFAMs.hmm", "$tmp_dir/TIGRFAMS_ROLE_LINK" , "$tmp_dir/TIGR_ROLE_NAMES", "$tmp_dir/TIGRFAMs.hmm.sqldb.tsv");
    }
    if(! -e "$tmp_dir/TIGRFAMS2GO.sqldb.tsv"){
	tigrfam_id2go_table("$tmp_dir/TIGRFAMS_GO_LINK", "$tmp_dir/TIGRFAMS2GO.sqldb.tsv");
    }
}

sub tigrfam_id2info_table{

    my ($hmm, $roleL, $roleN, $out) = @_;
    print STDERR "$hmm\n$roleL\n$roleN\n$out\n";
    my %fams = ();
    open(HMM, $hmm) or die "Could not open $hmm to read, $!\n";


    open(RLINK, $roleL);
    open(RNAMES, $roleN);
    open(OUT, ">$out") or die "Could not open $out to write, $!\n";

    #TIGR00001       158
    while(<RLINK>){
	chomp;
	my($ac, $roleid) = split(/\t/, $_);
	$fams{$ac}->{role_id} = $roleid;
    }
    close(RLINK);

    my %roleid2role = ();
    #role_id:        100     mainrole:       Central intermediary metabolism
    while(<RNAMES>){
	s/\://g;
	chomp;
	my @l = split(/\t/, $_);
	shift @l;
	$roleid2role{$l[0]}->{$l[1]} = $l[2];
    }
    close(RNAMES);


    $/ = "\n\//";

    while(<HMM>){
	if(/ACC\s+(\S.*?)\n.*LENG\s+(\d+)/s){
	    $fams{$1}->{LENG} = $2;
	}

    }

    my $dir = "$tmp_dir/TIGRFAMs_15.0_INFO";
    opendir(INPUTDIR,$dir) || die "Cannot opendir $dir: $!\n";

    my @dir_data=grep {!/^\.+$/} readdir(INPUTDIR);


    foreach my $name (@dir_data) {
	$/ = undef;
	open(NAME,"$dir/$name") or die "Cannot open file $dir/$name: $!\n";;
	my $data = <NAME>;

	my($ac) = $data =~ /\nAC\s+(TIGR\S+)/;

	if($data =~ /\nDE\s+(\S.*?)\n/){$fams{$ac}->{DE} = $1} else{$fams{$ac}->{DE} = "";}
	if($data =~ /\nIT\s+(\S.*?)\n/){$fams{$ac}->{IT} = $1} else{$fams{$ac}->{IT} = "";}
	if($data =~ /\nGS\s+(\S.*?)\n/){$fams{$ac}->{GS} = $1} else{$fams{$ac}->{GS} = "";}
	if($data =~ /\nEC\s+(\S.*?)\n/) {$fams{$ac}->{EC} = $1} else{$fams{$ac}->{EC} = "";}
	if($data =~ /^ID\s+(\S+?)\n/) {$fams{$ac}->{ID} = $1} else{$fams{$ac}->{ID} = "";}
	close(NAME);

    }

    close(INPUTDIR);

    foreach (sort keys %fams){
	my $mainrole = "";
	my $sub1role = "";
	if(exists $fams{$_}->{role_id}){
	    my $roleid = $fams{$_}->{role_id};
	    $mainrole = (exists $roleid2role{$roleid}->{"mainrole"}) ? $roleid2role{$roleid}->{"mainrole"} : "";
	    $sub1role = (exists $roleid2role{$roleid}->{"sub1role"}) ? $roleid2role{$roleid}->{"sub1role"} : "";
	}
	print OUT "$_\t";
	print OUT $fams{$_}->{ID}, "\t";
	print OUT $fams{$_}->{DE}, "\t";
	print OUT $fams{$_}->{IT}, "\t";
	print OUT $fams{$_}->{LENG},"\t";
	print OUT $fams{$_}->{GS},"\t";
	print OUT $fams{$_}->{EC},"\t";
	print OUT "$mainrole\t";
	print OUT "$sub1role\n";
    }
    close(OUT);
}

sub tigrfam_id2go_table{

    my ($golink, $out) = @_;
    open(GO, $golink) or  die "Could not open $golink to read, $!\n";
    open(OUT, ">$out") or die "Could not open $out to write, $!\n";
    while(<GO>){
	chomp;
	my @line =split(/\t/, $_);
	print OUT "$line[0]\t$line[1]\n";
    }
    close(GO);
    close(OUT);
}

sub build_go_db{

    my $file = "go_daily-termdb.rdf-xml.gz";
    if(! -e "$tmp_dir/$file"){

	#build a File::Fetch object
	my $ff = File::Fetch->new(uri => "http://archive.geneontology.org/latest-termdb/$file");
	### fetch the uri to local directory###
	msg("Fetching $file");
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }


    if(! -e "$tmp_dir/go_daily-termdb.rdf-xml"){
	msg("Uncompressing $tmp_dir/$file");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/$file");
	$ae->extract(to => "$tmp_dir");
    }
    go_table("$tmp_dir/go_daily-termdb.rdf-xml", "$tmp_dir/go.sqldb.tsv");

}

sub go_table{

    my ($xml_input, $out) = @_;
    open(OUT, ">$out") or die "Could not open $out to write, $!\n";

    # create object
    my $xml = new XML::Simple;
    # read XML file
    my $data = $xml->XMLin($xml_input);
    #print Dumper($data);

    #print "--go:accession\tgo:name\n";
    foreach my $term ( @{$data->{"rdf:RDF"}->{"go:term"}}){
	my $go_acc = $term->{"go:accession"};
	my $go_name =  $term->{"go:name"};
	next if $go_name eq "all";
	print OUT "$go_acc\t$go_name\n";
    }
    close(OUT);
}

sub build_foam_db{
    my @files = ("FOAM-hmm_rel1a.hmm.gz", "", "FOAM-onto_rel1.uniq.tsv");
    if(! -e ("$protein_hmm_dir/FOAM-hmm_rel1a.hmm" || "$tmp_dir/FOAM-hmm_rel1a.hmm.gz")){
	# build a File::Fetch object
	my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/FOAM-hmm_rel1a.hmm.gz");
	msg("Fetching FOAM-hmm_rel1a.hmm.gz");
	#fetch the uri to local directory
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	msg("Uncompressing $tmp_dir/FOAM-hmm_rel1a.hmm.gz");
        gunzip "$tmp_dir/FOAM-hmm_rel1a.hmm.gz" => "$protein_hmm_dir/FOAM-hmm_rel1a.hmm" or die "gunzip failed: $GunzipError\n";
	my $cmd = "$^X $bin/split_hmm.pl $protein_hmm_dir/FOAM-hmm_rel1a.hmm 5 $protein_hmm_dir FOAM";
	msg("Split FOAM database to multiple smaller database");
	runcmd($cmd);
    }
    if(! -e "$tmp_dir/FOAM-onto_rel1.uniq.tsv"){
	# build a File::Fetch object
	my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/FOAM-onto_rel1.uniq.tsv");
	msg("Fetching FOAM-onto_rel1.uniq.tsv");
	#fetch the uri to local directory
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }
}

sub build_genomedb{

    if(! -e "$diamond_dir/genomedb.dmnd"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/genomedb.faa.gz");
        msg("Fetching http://ebg.ucalgary.ca/metaerg/genomedb.faa.gz");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	my $cmd .= "diamond makedb --tmpdir /dev/shm --in $tmp_dir/genomedb.faa.gz -d $diamond_dir/genomedb";
	msg("Start running:$cmd");
	runcmd($cmd);
    }
    if(! -e "$tmp_dir/genomedb.tax.tsv"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/genomedb.tax.tsv");
        msg("Fetching http://ebg.ucalgary.ca/metaerg/genomedb.tax.tsv");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }

}
sub build_metabolic_hmmdb{

    if(! -e "$protein_hmm_dir/metabolic.hmm"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/metabolic.hmm.gz");
        msg("Fetching http://ebg.ucalgary.ca/metaerg/metabolic.hmm.gz");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	msg("Uncompressing $tmp_dir/metabolic.hmm.gz");
        gunzip "$tmp_dir/metabolic.hmm.gz" => "$protein_hmm_dir/metabolic.hmm" or die "gunzip failed: $GunzipError\n";

    }
}
sub build_casgene_hmmdb{

    if(! -e "$protein_hmm_dir/casgenes.hmm"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "http://ebg.ucalgary.ca/metaerg/casgenes.hmm.gz");
        msg("Fetching http://ebg.ucalgary.ca/metaerg/casgenes.hmm.gz");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	msg("Uncompressing $tmp_dir/casgenes.hmm.gz");
        gunzip "$tmp_dir/casgenes.hmm.gz" => "$protein_hmm_dir/casgenes.hmm" or die "gunzip failed: $GunzipError\n";
    }
}


sub build_sqlite_db{

    my ($sqlfile) = @_;

    use Cwd qw(getcwd);
    my $workdir = getcwd();

    msg("Start building sqldb $sqlfile");

    if(! -e "$sqlite_dir/metaerg.db"){

	my @table_input_file_names = ("pfam2go.sqldb.tsv", "Pfam-A.hmm.sql.db.tsv", "TIGRFAMS2GO.sqldb.tsv", "TIGRFAMs.hmm.sqldb.tsv", "go.sqldb.tsv", "sprot.sqldb.tsv", "FOAM-onto_rel1.uniq.tsv", "genomedb.tax.tsv");
	foreach (@table_input_file_names){
	    err("$tmp_dir/$_ does not exist") if ! -e  "$tmp_dir/$_";
	}
	msg("*******$sqlfile");
	my $path = abs_path($sqlfile);
	msg("*******$path");
	my $sqlfilename = basename($sqlfile);
	msg("Start copy sql table file sqlfile=$sqlfile, path=$path, filename=$sqlfilename");
	copy($sqlfile, $tmp_dir) or die "Copy failed: $!\n";
	chdir $tmp_dir;

   my $cdir = getcwd();
   msg("current dir=$cdir");
	my $cmd = "cat $sqlfilename | sqlite3 ../db/sqlite3/metaerg.db";
	runcmd($cmd);
	chdir $workdir;

    }


}

sub build_rRNAFinder_hmmdb{

    #annotated seed alignments in STOCKHOLM format
    my $rfam = "Rfam.seed";
    if(! -e "$tmp_dir/$rfam.gz"){
	my $ff = File::Fetch->new(uri => "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/$rfam.gz");

### fetch the uri to local directory###
	msg("Fetching $rfam.gz");
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }

    if(! -e "$tmp_dir/$rfam"){
	msg("Uncompressing $tmp_dir/$rfam.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/$rfam.gz");
	$ae->extract(to => $tmp_dir);
    }
    my @rna_types = ("5S_rRNA","5_8S_rRNA","SSU_rRNA_bacteria","SSU_rRNA_archaea","SSU_rRNA_eukarya","LSU_rRNA_archaea","LSU_rRNA_bacteria","LSU_rRNA_eukarya");
    my %rna_align_filehandlers = ();


    my %queries = (
	"5S_rRNA" => "$tmp_dir/5S_rRNA.rfam.align.fasta",
	"5_8S_rRNA" =>"$tmp_dir/5_8S_rRNA.rfam.align.fasta",
	"SSU_rRNA_bacteria" =>"$tmp_dir/SSU_rRNA_bacteria.rfam.align.fasta",
	"SSU_rRNA_archaea" =>"$tmp_dir/SSU_rRNA_archaea.rfam.align.fasta",
	"SSU_rRNA_eukarya" =>"$tmp_dir/SSU_rRNA_eukarya.rfam.align.fasta",
	"LSU_rRNA_archaea" =>"$tmp_dir/LSU_rRNA_archaea.rfam.align.fasta",
	"LSU_rRNA_bacteria" =>"$tmp_dir/LSU_rRNA_bacteria.rfam.align.fasta",
	"LSU_rRNA_eukarya" =>"$tmp_dir/LSU_rRNA_eukarya.rfam.align.fasta"
	);
    foreach (keys %queries){
	open $rna_align_filehandlers{$_}, ">", $queries{$_} or die "cannot open $queries{$_} to write, $!\n";
    }

    open(STOCK, "$tmp_dir/$rfam") or die "cannot open $rfam for reading: $!";

    my $i = 0;
    $/ = "\n\//";
    while(<STOCK>){
	chomp;
	if(/\#=GF\s+ID\s+(\S+)/s){
	    my $id = $1;
	    if(exists $queries{$id}){
		my @lines = split(/\n/, $_);
		foreach my $line (@lines){
		    if($line =~ /^\#/ || $line !~ /\S/){
			next;
		    }
		    else{
			#seqid alignment
			my @a = split(/\s+/, $line);
			#msg("id=$id");
			next if $a[0] eq "AY544260.1/501-381";
			my $handler = $rna_align_filehandlers{$id};
			print $handler ">$a[0]\n$a[1]\n";
		    }

		}
	    }

	}
    }

    $/ = "\n";
    close(STOCK);
    #foreach (keys %rna_align_filehandlers){
#	close($_);
 #   }

    fasta2domain("5S_rRNA.rfam.align.fasta", $tmp_dir);
    fasta2domain("5_8S_rRNA.rfam.align.fasta", $tmp_dir);

    my %hmms = (
	"arc_5SrRNA" => "$tmp_dir/arc_5S_rRNA.rfam.align.fasta",
	"euk_5SrRNA" => "$tmp_dir/euk_5S_rRNA.rfam.align.fasta",
	"bac_5SrRNA" => "$tmp_dir/bac_5S_rRNA.rfam.align.fasta",
	"euk_5_8SrRNA" => "$tmp_dir/euk_5_8S_rRNA.rfam.align.fasta",
	"arc_16SrRNA" => "$tmp_dir/SSU_rRNA_archaea.rfam.align.fasta",
	"bac_16SrRNA" => "$tmp_dir/SSU_rRNA_bacteria.rfam.align.fasta",
	"euk_18SrRNA" => "$tmp_dir/SSU_rRNA_eukarya.rfam.align.fasta",
	"arc_23SrRNA" => "$tmp_dir/LSU_rRNA_archaea.rfam.align.fasta",
	"bac_23SrRNA" => "$tmp_dir/LSU_rRNA_bacteria.rfam.align.fasta",
	"euk_28SrRNA" => "$tmp_dir/LSU_rRNA_eukarya.rfam.align.fasta"
	);

    my $cmd = "";
    foreach my $name (keys %hmms){
	my $alignFile = $hmms{$name};
	$cmd .= "hmmbuild --rna -n $name $tmp_dir/$name.hmm $alignFile;"
    }

    msg("Start running:$cmd");
    runcmd($cmd);

    $cmd = "cat $tmp_dir/arc_16SrRNA.hmm $tmp_dir/arc_23SrRNA.hmm $tmp_dir/arc_5SrRNA.hmm > $rna_hmm_dir/arc.hmm;";
    $cmd .= "cat $tmp_dir/bac_16SrRNA.hmm $tmp_dir/bac_23SrRNA.hmm $tmp_dir/bac_5SrRNA.hmm > $rna_hmm_dir/bac.hmm;";
    $cmd .= "cat $tmp_dir/euk_18SrRNA.hmm $tmp_dir/euk_28SrRNA.hmm $tmp_dir/euk_5_8SrRNA.hmm $tmp_dir/euk_5SrRNA.hmm > $rna_hmm_dir/euk.hmm;";

    msg("Start running:$cmd");
    runcmd($cmd);

    #my @tmp_files = glob "$tmp_dir/*";
    #unlink glob "'$tmp_dir/*'";
#rmdir $tmp_dir or warn "couldn't rmdir $tmp_dir: $!\n";
}

sub fasta2domain{
    my ($fasta, $tmp_dir) = @_;
    use Bio::DB::EUtilities;
    use Bio::SeqIO;

    open (FASTA, "$tmp_dir/$fasta") or die "Could not open $tmp_dir/$fasta to read, $!\n";
    $/ = "\n>";
    open (BAC, ">$tmp_dir/bac\_$fasta") or die "Could not open $tmp_dir\/bac\_$fasta for writting, $!\n";
    open (ARC, ">$tmp_dir\/arc\_$fasta") or die "Could not open $tmp_dir\/arc\_$fasta for writting, $!\n";;
    open (EUK, ">$tmp_dir\/euk\_$fasta") or die "Could not open $tmp_dir\/euk\_$fasta for writting, $!\n";;

    while(<FASTA>){
	chomp;
	if(my ($seqid,$other, $align) =  /^>?(\S+?)(\/.*?)\n(.*)/s){
	    
	    sleep(1);
	    #can get too many request problems. By default, the request from NCBI cannot be >3 times/second.
	    #with the api_key, it can be increased to 10 times/second perl email. 
	    
	    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						   -db      => 'nucleotide',
						   -rettype => 'gb',
						   -email   => 'xiaolid@gmail.com',
						   -api_key => 'f0f0c0d0dedaa42991657439b7c577b0ce08',
						   -id      => $seqid);
	    my $file = "$tmp_dir/myseqs.gb";
	    
	    # dump HTTP::Response content to a file (not retained in memory)
	    $factory->get_Response(-file => $file);
	    my $seqin = Bio::SeqIO->new(-file   => $file,
					-format => 'genbank');

	    while (my $seq = $seqin->next_seq) {
		my @classification = $seq->species->classification;
		my $lineage = join('\t', @classification);
		if($lineage =~ /Bacteria/i){
		    print  BAC ">$seqid$other\n$align\n";
		}
		elsif($lineage =~ /Archaea/i){
		    print  ARC ">$seqid$other\n$align\n";
		}
		elsif($lineage =~ /Eukaryota/i){
		    print EUK ">$seqid$other\n$align\n";
		}
	    }

	}

    }

    $/ = "\n";
    close(FASTA);
    close(BAC);
    close(ARC);
    close(EUK);
}

sub build_rRNAFinder_txondb{

    my $ssu = "SILVA\_$s_version\_SSURef_Nr99_tax_silva_trunc.fasta";
    my $lsu = "SILVA\_$s_version\_LSURef_tax_silva_trunc.fasta";

    my @files = ("$ssu.gz", "$lsu.gz", "README.txt");

    foreach my $file (@files){
	if(! -e "$tmp_dir/$file"){

            ### build a File::Fetch object ###
	    my $ff = File::Fetch->new(uri => "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/$file");

            ### fetch the uri to local directory###
	    msg("Fetching $file");
	    my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
	}
    }

    if(! -e "$tmp_dir/$ssu"){

	msg("Uncompressing $tmp_dir/$ssu.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/$ssu.gz");
	$ae->extract(to => $tmp_dir );
    }
    if(! -e "$tmp_dir/$lsu"){
	msg("Uncompressing $tmp_dir/$lsu.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/$lsu.gz");
	$ae->extract(to => $tmp_dir);
    }

    open(SSU,"$tmp_dir/$ssu") or die "Could not open $tmp_dir/$ssu to read, $!\n";
    open(SSU_DNA, ">$blast_dir/silva_SSURef_Nr99.fasta") || die "Could not open $blast_dir/silva_SSURef_Nr99.fasta to write, $!\n";
    msg("get ssu DNA file from $tmp_dir/$ssu to $blast_dir/silva_SSURef_Nr99.fasta for blast search database making");
    while(<SSU>){
	if (/^>(\S+?)\s+(\S+.*)$/){

	    print SSU_DNA ">$1 [$2]\n";
	}
	else{
	    s/U/T/g;
	    print SSU_DNA $_;
	}
    }
    close(SSU);
    close(SSU_DNA);

    open(LSU,"$tmp_dir/$lsu") or die "Could not open $tmp_dir/$lsu to read, $!\n";
    open(LSU_DNA, ">$blast_dir/silva_LSURef.fasta")  || die "Could not open $blast_dir/silva_LSURef.fasta to write, $!\n";
    msg("get lsu DNA file from $tmp_dir/$lsu to $blast_dir/silva_LSURef.fasta for blast search database making");
    while(<LSU>){
	if (/^>(\S+?)\s+(\S+.*)$/){

	    print LSU_DNA ">$1 [$2]\n";
	}

	else{
	    s/U/T/g;
	    print LSU_DNA $_;
	}
    }
    close(LSU);
    close(LSU_DNA);


    my $cmd = "makeblastdb -input_type fasta -dbtype nucl  -in $blast_dir/silva_SSURef_Nr99.fasta;";
    $cmd .= "makeblastdb -input_type fasta -dbtype nucl  -in $blast_dir/silva_LSURef.fasta";

    msg("Start running:$cmd");
    runcmd($cmd);



}


sub err {
    my($txt) = @_;
    msg($txt);
    exit(2);
}

sub msg {
    my $t = localtime;
    #my $line = "[".$t->hms."] @_\n";
    my $line = "[".$t->cdate."] @_\n";
    print STDERR $line;
}

sub runcmd {
    msg("Running:", @_);
    my $cmd = join("\n", @_);
    system(@_)==0 or err("Could not run command:$cmd, $!\n");
}
