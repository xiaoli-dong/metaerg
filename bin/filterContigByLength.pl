#!/usr/bin/env perl
use strict;
use warnings;
$#ARGV == -1

    and die "Usage: $0 <fasta file> <minSeqLength>\n";

my $fasta = $ARGV[0];
my $minlen = $ARGV[1];

open(FASTA, "<$fasta") or die "cannot open $fasta for reading: $!";

my $i = 0;

$/ = "\n>";
while(<FASTA>){
    chomp;
    
    if(my ($seq_name,$seq) =  /^>?(\S+.*?)\n(.*)/s){
	
	$seq =~ tr/ \r\n\t//d;
	
	if(length($seq) >=$minlen ){
	    
	    print  ">$seq_name\n$seq\n";
	}
    }
}

$/ = "\n";

close(FASTA);


