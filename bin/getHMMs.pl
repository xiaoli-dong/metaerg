#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($hmmFile, $type, $query);

$type = "ACC";
&GetOptions(
    "s=s" =>\$hmmFile,
    "t=s" =>\$type,
    "q=s" =>\$query
    );

($hmmFile) ||
    die "usage: $0 OPTIONS
where options are:\n -s  <input hmm file>\n -t <query type:ACC|file>\n -q <list of hmm names seperated by :[name1:name2...] or file with names in each line>\n";

my %names = ();

if($type eq "ACC"){
    %names = map { $_ => 1 } split(/:/, $query);
    
}
elsif($type eq "file"){
    
    open(NAMES, "<$query") or die $!;
    while(my $entry = <NAMES>){
	$entry =~ s/^>//;
	if(my ($hmmname) = $entry =~ /^(\S+)/){

	    $names{$hmmname} = 1;
	    print STDERR "$hmmname\n";
	}
	
	
    }
    close(NAMES);
}

open (HMMFILE, "$hmmFile") or die "Could not open $hmmFile to read, $!\n";

$/="//\n";
while(<HMMFILE>){
    if(/\nACC\s+(\S+)/s){ 
	my $name = $1;
	$name =~ s/\..*$//;
	if(exists $names{$name}){
	    print "$_";
	}
    }
    
}

close(HMMFILE);

