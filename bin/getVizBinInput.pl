#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my ($bindir);

&GetOptions(
    "d=s" =>\$bindir
    );

($bindir) ||
    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  add bin ids to the binned contigs and also produce an annotaton file, which has all the bin ids for all the contigs in order\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -d <a directory having all the interested bin files. the files are named in the format of Bin.*.fa>\n";

opendir(BINDIR, "$bindir") || die "Cannot opendir $bindir: $!\n";
open(BINNED, ">$bindir/binned_concat.fasta") or  die "Cannot open $bindir/binned_concat.fasta to write: $!\n";
open(ANNOT, ">$bindir/binned_annotation.list") or  die "Cannot open $bindir/binned_annotation.list to write: $!\n";
print ANNOT "label\n";

my @dirdata=grep {/^Bin\.\d+\.fa/} readdir(BINDIR);
foreach my $fasta (@dirdata) {
    my ($binid) = $fasta =~ /Bin\.(\d+)\.fa/;
  
    print STDERR "found binning file: $fasta\n";
    open (FASTA, "$bindir/$fasta")  || die "Could not open $bindir/$fasta file to read, $!\n";
    $/ = "\n>";
    while(<FASTA>){
	chomp;
	if(my ($seqid,$seq) =  /^>?(\S+).*?\n(.*)/s){
	    print BINNED ">$seqid $binid\n$seq\n";
	    print ANNOT "$binid\n";
	}
    }
    $/ = "\n";
    close(FASTA);
}
close(BINDIR);
close(BINNED);
close(ANNOT);

