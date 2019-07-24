#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec;

my @tools = ("aragorn", "rRNAFinder.pl", "rna2taxon.pl", "minced", "prodigal", "signalp", "tmhmm", "diamond", "hmmsearch", "sqlite3", "MinPath.py");

check_all_tools();



sub check_all_tools {
    for my $toolname (sort @tools) {
	check_tool($toolname);
    }
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# get (and check) versions of tools

sub check_tool {
    my($toolname) = @_;

    my $fp = find_exe($toolname);
    if(!$fp){
	print STDERR "Looking for '$toolname' ******** Can't find required '$toolname' in your \$PATH\n";

    }
    else{
	print STDERR "Looking for '$toolname' - found $fp\n";

    }
}
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
