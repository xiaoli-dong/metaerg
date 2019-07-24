#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Tools::GFF;
my ($fasta, $bfile);

&GetOptions(
    "f=s" =>\$fasta,
    "b=s" =>\$bfile
    );

($fasta && $bfile) ||
    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "Split input fasta into two files: good bin file contains all the sequences belong to user supplied list;\notherwise, the sequences will be written to the bad bin file\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -f <protein coding genes with binids as part of sequence id. e.g. bin64_LCM1|00019>\n".
    "  -b <a good bin list file: each line in the file represent one bin. e.g. /gpfs/ebg_projects/CCSA/sodaLake/DLM2_metaspades/bin3/DLM2.Bin.10.fa>\n";

my %goodbins = ();

open (BFILE, $bfile)  || die "Could not open $bfile to read, $!\n";

while(<BFILE>){
    chomp;
    my ($sname, $binid) = $_ =~ /(\w+)\.Bin\.(\d+)/;
    print STDERR "$_\t$binid\tbin$binid\_$sname\n";
    $goodbins{"bin$binid\_$sname"} = 1;
}
close(BFILE);

my ($prefix) = $fasta =~ /^(\S+)\.\w+$/;
open (GOOD, ">$prefix.good.faa") or die "Could not open $prefix.good.faa to read, $!\n";
open (BAD, ">$prefix.unbinned.faa") or die "Could not open $prefix.unbinned.faa to read, $!\n";

open (FASTA, $fasta) or die "Could not open $fasta to read, $!\n";

$/ = "\n>";
while(<FASTA>){
    chomp;
    if(/^>?(\S+?)(\|.*)/s){
	my $locus = $1;
	my $other = $2;
	print STDERR "locus=$locus\n";
	my ($binid, $sname) = split(/_/, $locus);
	
	if(exists $goodbins{$locus}){
	    print GOOD ">$locus$other\n";
	}
	else{
	    print BAD ">unbinned\_$sname$other\n";
	}
	
    }
}
$/ = "\n";

close(FASTA);
close(GOOD);
close(BAD);
