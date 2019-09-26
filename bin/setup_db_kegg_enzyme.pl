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
my $txt_dir = "$bin/../txt";
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

build_enzyme_table();
build_kegg_table();

sub build_enzyme_table{
    
    if(! -e "$tmp_dir/enzyme.dat"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat");
        msg("Fetching ftp://ftp.expasy.org/databases/enzyme/enzyme.dat");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }
    if(! -e "$tmp_dir/enzclass.txt"){
        # build a File::Fetch object
        my $ff = File::Fetch->new(uri => "ftp://ftp.expasy.org/databases/enzyme/enzclass.txt");
        msg("Fetching ftp://ftp.expasy.org/databases/enzyme/enzclass.txt");
        #fetch the uri to local directory
        my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }
    open(OUT, ">$tmp_dir/enzyme.sqldb.tsv") || die "Could not open $tmp_dir/enzyme.sqldb.tsv to write, $!\n";
    open(CLASS, "$tmp_dir/enzclass.txt") || die "Could not open $tmp_dir/enzclass.txt to read, $!\n";
    while(<CLASS>){
	chomp;
	if(/(^\d+?\.\s*[0-9\-]+\.\s*[0-9\-]+\.\s*[0-9\-]+?)\s+(\S.*)$/){
	    my $id = $1;
	    my $name = $2;
	    $id =~ s/\s+//g;
	    $name =~ s/^\s+//g;
	    $name =~ s/\s+$//g;
	    $name =~ s/\"/\"\"/g;
	    print OUT "$id\t$name\n";
	}
	
    }
    close(CLASS);
    
    $/= "\n\/\/";
    open(ENZYME, "$tmp_dir/enzyme.dat") || die "Could not open $tmp_dir/enzyme.dat to read, $!\n";
    
    while(<ENZYME>){
	
	if(/\nID\s+?(\S+?)\nDE\s+(\S.*?)\n/s){
	    my $id = $1;
	    my $name = $2;
	    $name =~ s/^\s+//g;
	    $name =~ s/\s+$//g;
	    $name =~ s/\"/\"\"/g;
	    print OUT "$id\t$name\n";
	}
	
    }
    close(ENZYME);
    close(OUT);
    $/= "\n";
}

sub build_kegg_table{
    
    open(KEGG, "$txt_dir/ko00001.keg") || die "Could not open $txt_dir/ko00001.keg to read, $!\n";
    open(OUT, ">$tmp_dir/ko.sqldb.tsv") || die "Could not open $tmp_dir/ko.sqldb.tsv to write, $!\n";
    
    my %kos = ();
    
    while(<KEGG>){
	chomp;
	if(/^D\s+(\S+?)\s+?(\S.*)$/){
	    my $koid = $1;
	    my $name = "unknown";
	    my $de = "";
	    
	    my @items = split(";", $2);
	    if(@items == 2){
		$name = $items[0];
		$de = $items[1];
	    }
	    else{
		$de = $2;
	    }
	    $name =~ s/^\s+//g;
	    $name =~ s/\s+$//g;
	    $name =~ s/\"/\"\"/g;
	    $de =~ s/\s+$//g; 
	    $de =~ s/^\s+//g;
	    $de =~ s/\"/\"\"/g;
	    
	    if(not exists $kos{$koid}){
		$kos{$koid}->{name} = $name;
		$kos{$koid}->{de} = $de;
	    }
	    
	}
	
    }
    close(KEGG);

    foreach my $id (sort keys %kos){
	print OUT "$id\t$kos{$id}->{name}\t$kos{$id}->{de}\n";
    }
    
    close(OUT);
    
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
