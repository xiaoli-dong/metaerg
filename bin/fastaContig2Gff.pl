#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use HTML::Entities;
use Bio::Tools::GFF;

my ($fasta_contig, $gff);

&GetOptions(
    "c=s" =>\$fasta_contig,
    "g=s" =>\$gff
    );

($fasta_contig && $gff) ||
    die "Name:\n".
    "by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  extract subset of the gff file for the input contigs\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -c <fasta format contig file>\n".
    "  -g <gff format metagenome annotation file generated using meta-annotator>\n";

open (QUERY, "$fasta_contig");
my %contigs = ();

while(<QUERY>){
    chomp;
    if(/^>(\S+)/){
	$contigs{$1} = 1;
	#print "$1\n";
    }
}
close(QUERY);

open(GFF, "$gff") or die "cannot open $gff for reading: $!";
my %genes = ();
while(<GFF>){
    if(/##gff-version/){
	print $_;
	next;
    }
    chomp;
    ###sequence-region NODE_38_length_181617_cov_73.6676 1 181617
    if(/^##sequence-region\s+(\S+)/){
	my $id = $1;
	if(exists $contigs{$id}){
	    print "$_\n";
	}
    }
    elsif(/^(\S+?)\s+.*?ID=(\S+?);/){
	my $contigid = $1;
	my $fid = $2;
	$genes{$fid} = 1;
	if(exists $contigs{$contigid}){
            print "$_\n";
        }
    }
    if(/^##FASTA/){
	print "$_\n";
	last;
    }
}

$/ = "\n>";
while(<GFF>){
    chomp;
    if(!/>?(\S+)\n(.+)/s){
        die "Could not read FastA record #$.: $_\n";
    }
    my $id = $1;
    my $seq = $2;
    #my ($contigid) = $id =~ /(^\S+?)\_cds/;
    if(exists $genes{$id}){
	print ">$id\n$seq\n";
    }
}

close(GFF);
