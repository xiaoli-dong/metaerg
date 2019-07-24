#!/usr/bin/env perl

use strict;
use warnings;
$#ARGV == -1

    and die "Usage: $0 <img_metadata.tsv: located under/export/data/programs/CheckM/data/img>  <concatenated tree generated from checkM> \n";


open(JGI, "$ARGV[0]") or die "Could not open $ARGV[0] to read, $!\n";
open(TREE, "$ARGV[1]") or die "Could not open $ARGV[1] to read, $!\n";

my %id2lineage = ();

while(<JGI>){
      
    chomp;
    next if /^taxon_oid/;
    my @l = split(/\t/, $_);
    my $lineage = join("_", @l[1,3,4,5,6,7,8,9]);
    my $id = $l[0];
    $id2lineage{$id} = $lineage;
    #print STDERR $lineage;
}

close(JGI);

while(<TREE>){
    my $l = $_;
    foreach my $id (keys %id2lineage){
	my $value = $id2lineage{$id};
	#print "$id=$value\n";
	if(/$id/){
	    $l =~ s/$id/$id\_$value/g;
	}
    }
    print $l;
}
