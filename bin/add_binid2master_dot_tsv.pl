#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Tools::GFF;
my ($bindir, $tsv);

&GetOptions(
    "d=s" =>\$bindir,
    "t=s" =>\$tsv
    );

($bindir && $tsv) ||
    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  add binid to the master.tsv file as the first column\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -d <directory name containing all the metabat generated fasta format bin files in the format of \"Bin.*.fa\">\n".
    "  -t <master.tsv file generated from metaerg, which is located in data directory>\n";

my %contig2bin = ();

opendir(BINDIR, "$bindir") || die "Cannot opendir $bindir: $!\n";
my @dirdata=grep {/^Bin\.\d+\.fa/} readdir(BINDIR);
foreach my $fasta (@dirdata) {
    my ($binid) = $fasta =~ /Bin\.(\d+)\.fa/;
    print STDERR "found binning file: $fasta\n";
    open (FASTA, "$bindir/$fasta")  || die "Could not open $bindir/$fasta file to read, $!\n";
    $/ = "\n>";
    while(<FASTA>){
	chomp;
	if(my ($seqid,$seq) =  /^>?(\S+).*?\n(.*)/s){
	    #print ">$binid\_$seqid\n$seq\n";
	    $contig2bin{$seqid} = $binid;
	}
    }
    $/ = "\n";
    close(FASTA);
}
close(BINDIR);

open (TSV, "$tsv") or die "Could not open $tsv file to read, $!\n";
while(<TSV>){
    if(/^#/){
	s/^#//;
	print "#binid\t$_";
    }
    elsif(/^(\S+)\s/){
	my $cid = $1;
	my $binid = "unbinned";
	$binid = $contig2bin{$cid} if exists $contig2bin{$cid};
	print "$binid\t$_\n";
    }
}
close(TSV);

