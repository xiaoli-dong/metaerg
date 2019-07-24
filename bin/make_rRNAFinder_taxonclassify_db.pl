#!/usr/bin/env perl
use File::Fetch;
use Archive::Extract;
use FindBin qw($Bin);

@ARGV == 1 or die "Usage: $0 <silva database version number, for example 132>\n",
      "download lsu and ssu database from silva\n";

my $version = $ARGV[0];

my $bin = "$FindBin::Bin";
my $tmpdir = "$bin/../db/tmp";
my $hmm_dir = "$bin/../db/hmm";
my $blastdb_dir = "$bin/../db/blast";


my $fp = find_exe("makeblastdb");
if(!$fp){
    err("Cannot find makeblastdb tool");
}

print STDERR "hmmdb_dir=$hmm_dir\n";
print STDERR "blastdb_dir=$blastdb_dir\n";


my $ssu = "SILVA\_$version\_SSURef_Nr99_tax_silva_trunc.fasta";
my $lsu = "SILVA\_$version\_LSURef_tax_silva_trunc.fasta";

@files = ("$ssu.gz", "$lsu.gz", "README.txt");

foreach $file (@files){

    print STDERR "$tmpdir/$file\n";
    if(! -e "$tmpdir/$file"){
### build a File::Fetch object ###
	my $ff = File::Fetch->new(uri => "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/$file");
### fetch the uri to local directory###
	print STDERR "Fetching $file\n";
	my $where = $ff->fetch(to => $tmpdir) or die $ff->error;
    }
}

if(! -e "$tmpdir/$ssu"){
    print STDERR "Uncompressing $tmpdir/$ssu.gz\n";
    my $ae = Archive::Extract->new(archive =>"$tmpdir/$ssu.gz");
    $ae->extract(to => $tmpdir );
}
if(! -e "$tmpdir/$lsu"){
    print STDERR "Uncompressing $tmpdir/$lsu.gz\n";
    my $ae = Archive::Extract->new(archive =>"$tmpdir/$lsu.gz");
    $ae->extract(to => $tmpdir);
}
open(SSU,"$tmpdir/$ssu") or die "Could not open $tmpdir/$ssu to read, $!\n";
open(SSU_DNA, ">$blastdb_dir/silva_SSURef_Nr99.fasta") || die "Could not open $blastdb_dir/silva_SSURef_Nr99.fasta to write, $!\n";
print STDERR "preparing $tmpdir/$lsu to $blastdb_dir/silva_SSURef_Nr99.fasta for blast search database making\n";
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

open(LSU,"$tmpdir/$lsu") or die "Could not open $tmpdir/$lsu to read, $!\n";
open(LSU_DNA, ">$blastdb_dir/silva_LSURef.fasta")  || die "Could not open $blastdb_dir/silva_LSURef.fasta to write, $!\n";
print STDERR "preparing $tmpdir/$ssu to $blastdb_dir/silva_LSURef.fasta for blast search database making\n";
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


my $cmd = "makeblastdb -input_type fasta -dbtype nucl  -in $blastdb_dir/silva_SSURef_Nr99.fasta;";
$cmd .= "makeblastdb -input_type fasta -dbtype nucl  -in $blastdb_dir/silva_LSURef.fasta";

print STDERR "Running $cmd\n";
system($cmd)==0 or err("Could not run command:$cmd");

sub find_exe {
    my($bin) = shift;
    for my $dir (File::Spec->path) {
	my $exe = File::Spec->catfile($dir, $bin);

	if(-x $exe){

	    return $exe;
	}
    }
    return;
}


