#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($seqFile, $type, $query);

$type = "ids";
&GetOptions(
    "s=s" =>\$seqFile,
    "t=s" =>\$type,
    "q=s" =>\$query
    );

($seqFile) ||
    die "usage: $0 OPTIONS
where options are:\n -s  <input sequence file>\n -t <query type:ids|file>\n -q <list of seqids seperated by :[id1:id2...] or file with ids in each line>\n";

my %names = ();

if($type eq "ids"){
    %names = map { $_ => 1 } split(/:/, $query);

}
elsif($type eq "file"){

    open(NAMES, "<$query") or die $!;
    while(my $entry = <NAMES>){
$entry =~ s/^>//;
	if(my ($seqname) = $entry =~ /^(\S+)/){

	    $names{$seqname} = 1;
	    print STDERR "$seqname\n";
	}


    }
    close(NAMES);
}

open(SEQ, "<$seqFile") or die "cannot open $seqFile for reading: $!";

my $i = 0;

$/ = "\n>";
while(<SEQ>){
    chomp;

    if(my ($seq_name,$other, $seq) =  /^>?(\S+)(.*?)\n(.*)/s){

	if(exists $names{$seq_name}){
	    print  ">$seq_name$other\n$seq\n";
	}

    }

}

$/ = "\n";

close(SEQ);
