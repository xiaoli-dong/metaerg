#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my ($tmpdir, $fasta_aa, $cpus);
$cpus = 8;
$tmpdir = ".";

&GetOptions(
    "d=s" =>\$tmpdir,
    "f=s" =>\$fasta_aa,
    "c=i" =>\$cpus
    );

($fasta_aa) ||
    die "usage: $0 OPTIONS
where options are:\n -d  <tmp directoy for the intermediate files>\n -f <fasta aa file>\n";


print STDERR "****************** =$tmpdir\n";
my $inputdir = "$tmpdir/phobius_input";
my $outputdir = "$tmpdir/phobius_output";


print STDERR "$inputdir, $outputdir\n";
if (-d $inputdir) {
    print STDERR "Folder '$inputdir' already exists, !$\n";
    exit;
}
elsif (-d $outputdir) {
    print STDERR "Folder '$outputdir' already exists, !$\n";
    exit;
}
else {
    print STDERR "Creating new tmp folder for run phobius parallel: $tmpdir";
    &run_cmds(1, "mkdir -p \Q$inputdir\E");
    &run_cmds(1, "mkdir -p \Q$outputdir\E");
}

split_fasta($fasta_aa, $cpus, $inputdir);

my @cmds = ();

opendir(INPUTDIR, "$inputdir") || die "Cannot opendir $inputdir: $!";

my @dir_data=grep {!/^\.+$/} readdir(INPUTDIR);

foreach my $name (@dir_data) {

    print STDERR "file $name found\n";
    my $query = "$inputdir/$name";
    my $output = "$outputdir/$name.phobius.out";
    my $cmd = "phobius.pl -short $query  > $output";
    push(@cmds, $cmd);
}
close(INPUTDIR);

&run_cmds($cpus, @cmds);

my $cmd = "cat $outputdir/*.phobius.out > $fasta_aa.phobius;";
&run_cmds(1, $cmd);

sub split_fasta{
    
    my($infile, $num_files, $output_dir) = @_;
    
    chomp ($infile);
    open (INFILE, "<$infile");
    my $total_seqs = 0;
    while(<INFILE>){
	chomp;
	$total_seqs++ if /^>.*?/;
    }
    print STDERR "$total_seqs sequences we have\n";
    close(INFILE);
    
    my $remainder = $total_seqs % $num_files;
    my $seq_per_file = ($total_seqs - $remainder)/$num_files;
    
    open (INFILE, "<$infile");
    $/ = "\n>";
    my $count = 0;
    my $file_count = 1;
    my $output_file = "$output_dir/$file_count.fasta";
    open (OUTPUT, ">$output_file") || die "cannot open the file $output_file\n";
        
#split input file into multiple files
    
    while(my $entry = <INFILE>){
	
	chomp($entry);
	
	if(my ($name, $seq) = $entry =~ /^>?(\S+).*?\n(.*)/s){
	    print OUTPUT ">$name\n$seq\n";
	    $count++;
	}
	
	if($file_count != $num_files && $count == $seq_per_file){
	    close(OUTPUT);
	    $count = 0;
	    
	    if($file_count < $num_files){
		$file_count++;
		$output_file = "$output_dir/$file_count.fasta";
		open (OUTPUT, ">$output_file") || die "cannot open the file $output_file\n";
	    }
	}
	
    }
    
    close(INFILE);
    close(OUTPUT);
    
    
    
}
sub run_cmds {

        my ($max_cmds, @cmds) = @_;

        my ($num_children, $pid);

        return unless @cmds; # get out of the function if there is nothing in @cmds

        for($num_children = 1; $num_children < $max_cmds && @cmds; $num_children++){
        # initialize the number of child processes at 1, and increment it by one
        #while it is less than $max_cmds

                my $cmd = shift (@cmds);
                if($pid = fork) {
                        # do nothing if parent
                } elsif(defined $pid) { # $pid is zero here if defined
                        system $cmd;
                        exit;
                } else {
                        #weird fork error
                        die "Can't fork: $!\n";
                }
        }

        while(@cmds) {
                undef $pid;
                FORK: {
                        my $cmd = shift (@cmds);
                        if($pid = fork) {
                                # parent here
                                $num_children++;
                                wait;
                                $num_children--;
                                next;

                        } elsif(defined $pid) { # $pid is zero here if defined
                                system $cmd;
                                exit;

                        } else {
                                #weird fork error
                                die "Can't fork: $!\n";
                        }
                }
        }
        wait while $num_children--;
}
